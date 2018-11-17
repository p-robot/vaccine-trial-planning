#include <ctime>
#include <cmath>
#include <memory>
#include <numeric>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/assert.hpp>

#include "configuration.h"
#include "population.h"
#include "population_statistics.h"
#include "event.h"
#include "event_queue.h"
#include "random_number_generator.h"
#include "scribe.h"
#include "serotype.h"
#include "simulation.h"
#include "trial_arm.h"
#include "vaccination_program.h"

namespace fs = boost::filesystem;

Simulation::Simulation(int seed, std::string config_file, std::string output_folder) 
  : seed_                    { seed }
  , calendar_ptr_            { std::make_shared<Calendar>() }
  , rng_ptr_                 { std::make_shared<RandomNumberGenerator>(seed_) }
  , population_              { Population {calendar_ptr_, rng_ptr_, config_file} }
  , event_queue_             { }
  , start_time_              { std::time(nullptr) }
  , scribe_                  { output_folder }
  , trial_scribe_            { (fs::path(output_folder) / "trial_logs").string() }
  , trial_prevalence_scribe_ { (fs::path(output_folder) / "trial_prevalence").string() }
  , trial_exposure_scribe_   { (fs::path(output_folder) / "trial_exposure").string() }

{
  // define other members from configuration file
  configure(config_file);

  // copy configuration files and save the seed used
  fs::path source = (fs::path {config_file}).parent_path();
  scribe_.copy_folder(source.string());  
  scribe_.note_value(seed_, "seed.txt");
}

void Simulation::configure(std::string config_file) {
  boost::property_tree::ptree ptree = load_config(config_file, "simulation");
  num_initial_hosts_      = ptree.get<int>("num_initial_hosts");
  num_past_colonizations_ = ptree.get<int>("num_past_colonizations");
  num_days_demographic_   = 365 * ptree.get<int>("num_years_burn_in.demographic");
  num_days_epidemiologic_ = 365 * ptree.get<int>("num_years_burn_in.epidemiologic");
  num_days_simulation_    = 365 * ptree.get<int>("num_years_simulation");
  num_days_total_         = num_days_demographic_ + num_days_epidemiologic_ + num_days_simulation_;
  
  auto vac_programs_array = ptree.get_child("vaccinations");
  for (auto& v : vac_programs_array) {
    auto vaccine = v.second.get<std::string>("vaccine");
    auto start_day = 365 * v.second.get<int>("start_year") + v.second.get<int>("start_day");
    auto end_day = start_day + 365 * v.second.get<int>("duration") - 1;
    auto coverage = as_vector<double>(v.second.get_child("coverage"));
    auto schedule = as_vector<int>(v.second.get_child("schedule"));
    vaccination_programs_.emplace_back(vaccine, start_day, end_day, schedule, coverage);
  }
  
  auto trial_arms_array = ptree.get_child("trial_arms");
  for (auto arm : trial_arms_array) {
    auto name = arm.second.get<std::string>("name");
    auto num_subjects = arm.second.get<int>("num_subjects");
    auto start_day = 365 * arm.second.get<int>("start_year");
    auto vaccine = arm.second.get<std::string>("vaccine");
    auto schedule = as_vector<int>(arm.second.get_child("schedule"));
    trial_arms_.emplace(std::make_pair(name, TrialArm(name, num_subjects, start_day, vaccine, schedule)));
  }
}

void Simulation::set_population_beta(double beta) {
  population_.set_beta(beta);
}

double Simulation::get_population_beta() const {
  return population_.get_beta();
}

std::vector<Simulation::Snapshot> Simulation::run() {
  auto start_time = std::time(nullptr);
  std::cout << std::setprecision(2) << std::fixed; // format std::cout so our progress updates look nice
  std::cout << "Starting simulation... " << std::endl;

  std::vector<Snapshot> snapshots;
  initialize_population();
  while (true) {
    int current_day = calendar_ptr_->get_current_day();
    
    // check end condition
    if (current_day > num_days_total_)
      break;

    // determine which part of the simulation we're in
    auto mode = get_mode(current_day);

    // note our progress at the beginning of every year
    if (current_day % 365 == 0) {
      save_statistics();
      std::cout << "\r" << 100 * current_day / (double) num_days_total_ << "%";
      std::cout.flush();

      // if in Simulation mode, save snapshots to return later
      if (mode == Mode::Simulation) {
        auto snapshot = take_snapshot();
        snapshots.emplace_back(snapshot);
      }
    }
    
    // seed colonizations if we just finished the demographic-only part of the simulation
    if (current_day == num_days_demographic_) {
      population_.seed_past_colonizations(num_past_colonizations_);
      initialize_colonizations();
    }
    
    // trial arms
    for (auto& kv : trial_arms_) {
      auto& arm = kv.second;
      int ds = current_day - arm.get_start_day();
      
      // start a trial arm if it's time
      if (ds == 0) {
        initialize_trial_arm(arm);
      }
      
      // log prevalences and cumulative exposure in each trial arm every 30 days
      if (ds >= 0 && ds % 30 == 0) {
        int dv = ds - arm.get_last_vaccination_day();
        
        // point prevalence
        int n = population_.stats->get_num_participants(arm.get_name());
        int m = population_.stats->get_num_colonized_participants(arm.get_name());
        std::vector<int> row_prev { ds, dv, n, m };
        trial_prevalence_scribe_.note_values(row_prev, arm.get_name() + ".csv");

        // distribution of number of past colonizations
        auto exposure_dist = population_.stats->get_num_past_colonizations_dist(arm.get_name(), 50);
        std::vector<int> row_exp { ds, dv };
        row_exp.insert( row_exp.end(), begin(exposure_dist), end(exposure_dist) );
        trial_exposure_scribe_.note_values(row_exp, arm.get_name() + ".csv");
      }
    }
    
    // run one day of the simulation
    run_one_day(mode, current_day);
    
    // verify 5 times throughout
    int verify_interval = num_days_simulation_ / 5;
    int days_into_simulation = current_day - (num_days_demographic_ + num_days_epidemiologic_);
    if (mode == Simulation::Mode::Simulation && days_into_simulation % verify_interval == 0) {
      population_.verify();
    }
    
    // move to tomorrow
    calendar_ptr_->increment_current_day();
  }
  
  auto num_secs = std::time(nullptr) - start_time;
  std::cout << "...finished. Took ";
  if (num_secs < 60) 
    std::cout << num_secs << " seconds." << std::endl;
  else
    std::cout << num_secs / 60 << " minutes and " << num_secs % 60 << " seconds." << std::endl;
  return snapshots;
}

void Simulation::initialize_population() {
  // add initial population and schedule life events
  auto ids = population_.initialize(num_initial_hosts_);
  for (int id : ids) {
    auto host_ptr = population_.find_host(id);
    schedule_life_events(id);
  }
  // flush event queue of events in the past
  while (event_queue_.size() > 0 && event_queue_.day_of_next_event() < calendar_ptr_->get_current_day()) {
    event_queue_.pop();
  }
}

void Simulation::initialize_colonizations(){
  // seed initial infections and add recovery events to the queue
  int current_day = calendar_ptr_->get_current_day();
  auto host_serotype_pairs = population_.seed_colonizations();
  for (auto& pair : host_serotype_pairs) {
    schedule_colonization(pair.first, current_day, pair.second);
  }
}

void Simulation::initialize_trial_arm(TrialArm arm) {
  for (int i = 0; i < arm.get_num_subjects(); i++) {
    int host_id = population_.birth_host(arm.get_name());
    schedule_life_events(host_id);
  }
}

Simulation::Mode Simulation::get_mode(int day) const {
  if (day < num_days_demographic_) {
    return Simulation::Mode::Demographic;
  } else if (day < num_days_demographic_ + num_days_epidemiologic_) {
    return Simulation::Mode::Epidemiologic;
  } else {
    return Simulation::Mode::Simulation;
  }
}

int Simulation::run_one_day(Simulation::Mode mode, int current_day) {
  if (mode == Simulation::Mode::Epidemiologic || mode == Simulation::Mode::Simulation) {
    // colonize and schedule recovery events 
    auto host_serotype_pairs = population_.determine_colonizations();
    for (auto& pair : host_serotype_pairs) {
      schedule_colonization(pair.first, current_day, pair.second);
    }
  }

  // process events for today
  int num_events_processed = 0;
  while (event_queue_.size() > 0 && event_queue_.day_of_next_event() <= current_day) {
    Event event = event_queue_.pop();
    // check that this event isn't in the past somehow
    if (event.get_day() < current_day) {
      throw std::runtime_error { "Simulation::run_one_day: Next event occurs before current day. This should not occur." };
    }
    // log the event if it's relevant for the trial (before we process it, in case of death events)
    if (is_relevant_for_trial(event)) {
      log_for_trial(event);
    }
    // now process the event
    process_event(event, current_day);
    num_events_processed++;
  }

  return num_events_processed;
}

void Simulation::schedule_birthday(int host_id, int birthday) {
  event_queue_.schedule_event(Event { host_id, birthday, Event::Type::Birthday, true });
}
void Simulation::schedule_colonization(int host_id, int day_of_colonization, int serotype) {
  event_queue_.schedule_event(Event { host_id, day_of_colonization, Event::Type::Colonization, serotype });
}
void Simulation::schedule_recovery(int host_id, int day_of_recovery, int serotype, int day_of_colonization) {
  if (population_.has_host(host_id)) {
    auto host_ptr = population_.find_host(host_id);
    int day_of_death = host_ptr->get_day_of_death();
    if (day_of_recovery < day_of_death) {
      event_queue_.schedule_event(Event { host_id, day_of_recovery, Event::Type::Recovery, serotype, day_of_colonization });
    }
  }
}
void Simulation::schedule_death(int host_id) {
  auto host_ptr = population_.find_host(host_id);
  event_queue_.schedule_event(Event { host_id, host_ptr->get_day_of_death(), Event::Type::Death });
}
void Simulation::schedule_vaccination(int host_id, int day_of_vaccination, std::string vaccine_name) {
  event_queue_.schedule_event(Event { host_id, day_of_vaccination, Event::Type::Vaccination, vaccine_name });
}

void Simulation::schedule_life_events(int host_id) {
  auto host_ptr = population_.find_host(host_id);
  auto day_of_birth = host_ptr->get_day_of_birth();
  
  if (day_of_birth != calendar_ptr_->get_current_day()) {
    std::runtime_error { "Simulation::schedule_life_events: Scheduling events should occur on the host's day of birth." };
  }
  
  // schedule birthdays
  for (int i = 1; i <= host_ptr->get_lifespan(); i++) {
    int birthday =  day_of_birth + i * 365;
    schedule_birthday(host_id, birthday);
  }

  // schedule death
  schedule_death(host_id);
  
  // schedule vaccinations
  if (host_ptr->is_in_trial()) {
    // vaccination based on trial arm protocol
    auto arm = trial_arms_.at(host_ptr->get_trial_arm_name());
    for (auto d: arm.get_schedule()) {
      schedule_vaccination(host_id, day_of_birth + d, arm.get_vaccine_name());
    }
  } else {
    // vaccination based on programs for general population
    for (auto& program : vaccination_programs_) {
      auto coverage = program.get_coverage(day_of_birth);
      if (rng_ptr_->draw_bernoulli(coverage)) {
        for (auto d : program.get_schedule()) {
          schedule_vaccination(host_id, day_of_birth + d, program.get_vaccine_name());
        }
      }
    }
  }
}

void Simulation::process_event(const Event& event, int current_day) {
  switch (event.get_type()) {
    case Event::Type::Birthday: {
      population_.update_host_age_group(event.get_host_id());
      break;
    }
    case Event::Type::Death: {
      bool in_trial = population_.find_host(event.get_host_id())->is_in_trial();
      population_.remove_host(event.get_host_id());
      
      // birth a replacement ONLY if the host was not in a trial
      if (!in_trial) {
        int host_id = population_.birth_host();
        schedule_life_events(host_id);
      }
      break;
    }
    case Event::Type::Recovery: {
      population_.recover_host(event.get_host_id());
      break;
    }
    case Event::Type::Colonization: {
      int day_of_recovery = population_.colonize_host(
        event.get_host_id(), event.get_serotype()
      );
      schedule_recovery(event.get_host_id(), day_of_recovery, event.get_serotype(), current_day);
      break;
    }
    case Event::Type::Vaccination: {
      population_.vaccinate_host(event.get_host_id(), event.get_vaccine_name());
      break;
    }
  }
}

Simulation::Snapshot Simulation::take_snapshot() const {
  //auto n = population_.get_num_hosts_by_age();
  //int num_hosts_under_5 = std::accumulate(n.begin(), n.begin() + 5, 0);
  //std::vector<int> num_colonized_under_5 = population_.sample_hosts(std::make_pair(0, 4));
  return {
    population_.stats->get_num_hosts_u5(),
    population_.stats->get_num_colonized_u5()
  };  
}

void Simulation::save_statistics() {
  scribe_.note_value(calendar_ptr_->get_current_day(), "sampling_days.csv");
  scribe_.note_value(event_queue_.size(),              "num_events_queued.csv");

  scribe_.note_value(population_.stats->get_num_hosts(), "num_hosts.csv");
  scribe_.note_value(population_.stats->get_num_colonized(), "num_colonized.csv");
  scribe_.note_value(population_.stats->get_num_colonized(std::make_pair(0, 4)), "num_colonized_under_5.csv");
  scribe_.note_values(population_.stats->get_num_hosts_by_age(), "num_hosts_by_age.csv");
  scribe_.note_value(population_.stats->get_num_hosts(std::make_pair(390, 419)), "num_hosts_13mo.csv" );

  for (int s = 0; s < Serotype::get_num_serotypes(); s++) {
    scribe_.note_values(population_.stats->get_num_colonized_by_age(s, false), "num_colonized_by_age_" + Serotype::get_name(s) + ".csv");
    scribe_.note_values(population_.stats->get_num_colonized_by_age(s, true), "num_colonized_by_age_ss_" + Serotype::get_name(s) + ".csv");
//    scribe_.note_values(population_.get_num_previously_colonized_by_age(s), "num_previously_colonized_by_age_" + Serotype::get_name(s) + ".csv");
  }

  for (auto program : vaccination_programs_) {
    auto vaccine_name = program.get_vaccine_name();
    int num_doses = program.get_schedule().size();
    scribe_.note_value(population_.stats->get_num_vaccinated(vaccine_name, num_doses, std::make_pair(390, 419)),
                       "num_vaccinated_13mo_" + vaccine_name + ".csv");
    scribe_.note_value(population_.stats->get_num_vaccinated(vaccine_name, num_doses, std::make_pair(0, 1824)),
                       "num_vaccinated_under_5_" + vaccine_name + ".csv");
    scribe_.note_value(population_.stats->get_num_vaccinated(vaccine_name, num_doses),
                       "num_vaccinated_all_ages_" + vaccine_name + ".csv");
  }

  scribe_.note_values(population_.stats->get_num_past_colonizations_dist(std::make_pair(0, 4), 20),
                      "num_past_colonizations_dist_under_5.csv");
  scribe_.note_values(population_.stats->get_num_past_colonizations_dist(std::make_pair(25, 49), 20),
                      "num_past_colonizations_dist_25_to_50.csv");
  scribe_.note_values(population_.stats->get_num_serotypes_dist(std::make_pair(0, 4), 5),
                      "num_serotypes_dist_under_5.csv");
}

bool Simulation::is_relevant_for_trial(const Event &event) const {
  // relevant means host is alive, in a trial, and the event is not a Birthday
  auto host_id = event.get_host_id();
  if (!population_.has_host(host_id)) {
    return false; // host has likely died
  } else {
    auto host_ptr = population_.find_host(event.get_host_id());
    return host_ptr->is_in_trial() && (event.get_type() != Event::Type::Birthday);
  }
}

void Simulation::log_for_trial(const Event &event) {
  // this method assumes host is still alive
  auto host_id = event.get_host_id();
  if (!population_.has_host(host_id)) {
    throw std::runtime_error { "Simulation::log_for_trial: Host not found." };
  }
  auto host_ptr = population_.find_host(host_id);
  auto trial_arm_name = host_ptr->get_trial_arm_name();
  
  // log event in a file
  trial_scribe_.note_value(event.as_csv_string(true), trial_arm_name + ".csv");
}
