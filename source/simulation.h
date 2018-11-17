#pragma once
#include <memory>
#include <vector>
#include <ctime>
#include <string>
#include <unordered_map>

#include "event.h"
#include "event_queue.h"
#include "population.h"
#include "random_number_generator.h"
#include "scribe.h"
#include "trial_arm.h"
#include "vaccination_program.h"

class Simulation {
public:
  struct Snapshot {
    int num_hosts_under_5;
    std::vector<double> num_colonized_under_5; // indexed by serotype
  };

  enum class Mode {
    Demographic,
    Epidemiologic,
    Simulation
  };

  Simulation(int seed, std::string config_path, std::string output_dir);
  void set_population_beta(double beta);
  double get_population_beta() const;
  std::vector<Snapshot> run();

private:
  void configure(std::string config_path);

  void initialize_population();
  void initialize_colonizations();
  void initialize_trial_arm(TrialArm arm);
  Mode get_mode(int day) const;
  int run_one_day(Simulation::Mode mode, int current_day); // returns number of events processed
  void schedule_birthday(int host_id, int birthday);
  void schedule_colonization(int host_id, int day_of_colonization, int serotype);
  void schedule_recovery(int host_id, int day_of_recovery, int serotype, int day_of_colonization);
  void schedule_death(int host_id);
  void schedule_vaccination(int host_id, int day_of_vaccination, std::string vaccine);
  void schedule_life_events(int host_id);
  void process_event(const Event& event, int current_day);
  Snapshot take_snapshot() const;
  void save_statistics();
  
  // vaccination trial-related
  bool is_relevant_for_trial(const Event& event) const; // check if event is relevant (host is alive, in trial)
  void log_for_trial(const Event& event);

  int seed_;
  std::shared_ptr<Calendar> calendar_ptr_;
  std::shared_ptr<RandomNumberGenerator> rng_ptr_;
  Population population_;
  EventQueue event_queue_;
  
  int num_past_colonizations_;
  int num_days_demographic_;
  int num_days_epidemiologic_;
  int num_days_simulation_;
  int num_days_total_;
  int num_initial_hosts_;
  std::vector<VaccinationProgram> vaccination_programs_;
  std::unordered_map<std::string, TrialArm> trial_arms_;
  
  std::time_t start_time_;
  Scribe scribe_;
  Scribe trial_scribe_;
  Scribe trial_prevalence_scribe_;
  Scribe trial_exposure_scribe_;
};
