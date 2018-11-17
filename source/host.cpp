#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

#include "configuration.h"
#include "colonization.h"
#include "host.h"
#include "random_number_generator.h"
#include "serotype.h"
#include "vaccine.h"

int Host::next_unused_id { 0 };
double Host::mu_;
double Host::epsilon_;
//double Host::sigma_;
std::vector<double> Host::sigmas_;
double Host::sigma_sg_;
double Host::decay_st_;
double Host::decay_sg_;

class ColonizationRankCompare {
public:
  bool operator()(const Colonization& c1, const Colonization& c2) {
    return Serotype::get_rank(c1.get_serotype()) < Serotype::get_rank(c2.get_serotype());
  }
};

Host::Host(int day_of_birth, int day_of_death,
           std::shared_ptr<const Calendar> calendar_ptr,
           std::shared_ptr<RandomNumberGenerator> rng_ptr,
           std::string trial_arm_name)
  : day_of_birth_             { day_of_birth }
  , day_of_death_             { day_of_death }
  , calendar_ptr_             { calendar_ptr }
  , rng_ptr_                  { rng_ptr }
  , trial_arm_name_           { trial_arm_name }
  , in_trial_                 { trial_arm_name != "" }
  , total_past_colonizations_ { 0 }
  , num_colonizations_        { std::vector<int>   (Serotype::get_num_serotypes(), 0) }
  , num_past_colonizations_   { std::vector<int>   (Serotype::get_num_serotypes(), 0) }   
  , serotypes_cleared_        { std::vector<bool>  (Serotype::get_num_serotypes(), false) }
  , serotype_clearance_dates_ { std::vector<int>   (Serotype::get_num_serotypes(), -1) }
  , susceptibilities_         { std::vector<double>(Serotype::get_num_serotypes(), 1.0) }
  , infectiousness_           { std::vector<double>(Serotype::get_num_serotypes(), 1.0) }
  , duration_factors_         { std::vector<double>(Serotype::get_num_serotypes(), 1.0) }
{
  id_ = next_unused_id++;
}

// Demographic-related
int Host::get_id() const {
  return id_;
}
int Host::get_day_of_birth() const {
  return day_of_birth_;
}
int Host::get_day_of_death() const {
  return day_of_death_;
}
int Host::get_lifespan() const {
  return (day_of_death_ - day_of_birth_) / 365;
}
int Host::calculate_age() const {
  return (calendar_ptr_->get_current_day() - day_of_birth_) / 365; // age in years
} 
int Host::calculate_age_in_days() const {
  return calendar_ptr_->get_current_day() - day_of_birth_;
}

// Colonization-related
int Host::get_num_colonizations(int serotype) const {
  return num_colonizations_[serotype];
}
int Host::get_num_colonizations() const {
  return colonizations_.size();
}
int Host::get_num_past_colonizations(int serotype) const {
  return num_past_colonizations_[serotype];
}
int Host::get_num_past_colonizations() const{
  return total_past_colonizations_;
}
int Host::get_num_serotypes() const {
  int num_serotypes = 0;
  for (auto& n: num_colonizations_) {
    num_serotypes += (int) (n > 0);
  }
  return num_serotypes;
}
void Host::verify_num_colonizations() const {
  auto counts = std::vector<int>(Serotype::get_num_serotypes(), 0);
  for (const Colonization& c : colonizations_) {
    counts[c.get_serotype()]++;
  }
  for (int s = 0; s < Serotype::get_num_serotypes(); s++) {
    if (num_colonizations_[s] != counts[s]) {
      throw std::runtime_error {"Host::verify_num_colonizations: Number of colonizations may be incorrect."};
    }
  }
}

double Host::get_susceptibility(int serotype) const {
  return susceptibilities_[serotype];
}
double Host::get_infectiousness(int serotype) const {
  return infectiousness_[serotype];
}
double Host::get_duration_factor(int serotype) const {
  return duration_factors_[serotype];
}
void Host::acquire_past_colonizations(int serotype, int n) {
  num_past_colonizations_[serotype] += n;
  serotypes_cleared_[serotype] = true;
  serogroups_cleared_.insert(Serotype::get_serogroup(serotype));
  total_past_colonizations_ += n;
  update_susceptibilities();
} 
Colonization Host::sample_colonization() const { // caller: first check num_colonizations > 0
  if (colonizations_.empty())
    throw std::runtime_error {"Host::sample_colonization: No colonizations to sample."};
  int index = rng_ptr_->draw_uniform_int(0, colonizations_.size() - 1);
  return colonizations_[index];
}
Colonization Host::acquire_colonization(int serotype) {
  int day_of_colonization = calendar_ptr_->get_current_day();
  int day_of_recovery = day_of_colonization + calculate_colonization_duration(serotype);
  Colonization c { serotype, day_of_colonization, day_of_recovery, rng_ptr_ };
  colonizations_.push_back(c);
  num_colonizations_[serotype]++;
  update_susceptibilities();
  return c;
}
std::vector<int> Host::remove_expired_colonizations() {
  std::vector<int> removed_serotypes;
  int current_day = calendar_ptr_->get_current_day();
  auto it = colonizations_.begin();
  while (it != colonizations_.end()) {
    if (it->get_day_of_recovery() <= current_day) {
      int serotype = it->get_serotype();
      it = colonizations_.erase(it); // erase advances the iterator

      num_colonizations_[serotype]--;
      num_past_colonizations_[serotype]++;
      total_past_colonizations_++;
      serotypes_cleared_[serotype] = true;
      serotype_clearance_dates_[serotype] = current_day;
      serogroups_cleared_.insert(Serotype::get_serogroup(serotype));
      serogroup_clearance_dates_[Serotype::get_serogroup(serotype)] = current_day;
      removed_serotypes.push_back(serotype);
    } else {
      it++; // advance the iterator ourselves
    }
  }
  update_susceptibilities();
  return removed_serotypes;
}

// Vaccination-related
void Host::acquire_vaccination(std::string vaccine_name) {
  if (vaccination_history_.count(vaccine_name) == 0) {
    vaccination_history_.insert(std::make_pair(vaccine_name, 1));
  } else {
    ++vaccination_history_[vaccine_name];
  }
  update_susceptibilities();
  update_infectiousness();
  update_duration_factors();
}
int Host::get_num_vaccinations(std::string vaccine_name) const {
  if (vaccination_history_.count(vaccine_name) == 0) {
    return 0;
  } else {
    return vaccination_history_.at(vaccine_name);
  }
}

// Vaccine-trial related
std::string Host::get_trial_arm_name() const {
  return trial_arm_name_;
}
bool Host::is_in_trial() const {
  return in_trial_;
}

// private
int Host::calculate_colonization_duration(int serotype) {
  double gamma = Serotype::get_duration(serotype); // intrinsic mean duration absent immunity
  double kappa = Serotype::get_base_duration();    // minimum mean duration with immunity
  double fraction = exp(-epsilon_ * total_past_colonizations_) * duration_factors_[serotype];
  double mean_duration = kappa + (gamma - kappa) * fraction;
  int num_days = (int) rng_ptr_->draw_exponential(mean_duration) + (int) rng_ptr_->draw_bernoulli(0.5); // round up or down with equal probability
  return num_days;
}

// assumes carriage is up-to-date. 
void Host::update_susceptibilities() {
  // determine exclusion by fittest strain
  double omega;
  if (colonizations_.empty()) {
    omega = 0.0;
  } else {
    ColonizationRankCompare cmp;
    auto it = std::min_element(colonizations_.begin(), colonizations_.end(), cmp);
    double min_f = Serotype::get_rank(it->get_serotype());
    int Z = Serotype::get_num_serotypes(); 
    omega = mu_ * (1.0 - (min_f - 1.0) / (Z - 1.0)); // Cobey 2012, S5. 
  }

  // determine specific immunity
  for (int s = 0; s < Serotype::get_num_serotypes(); s++) {
    // serotype-specific immunity
    double st_immunity = 0.0;
    if (serotypes_cleared_[s]) {
      int num_days_ago = calendar_ptr_->get_current_day() - serotype_clearance_dates_[s];
      st_immunity = sigmas_[s] * std::exp(-decay_st_ * num_days_ago);
    }
    
    // serogroup-specific immunity
    double sg_immunity = 0.0;
    if (serogroups_cleared_.count(Serotype::get_serogroup(s)) > 0) {
      int num_days_ago = calendar_ptr_->get_current_day() - serogroup_clearance_dates_[Serotype::get_serogroup(s)];
      sg_immunity = sigma_sg_ * std::exp(-decay_sg_ * num_days_ago);
    }

    // vaccine-induced immunity
    double vc_immunity = find_max_reduction(Vaccine::EffectType::Susceptibility, s);
    
    // take maximum from all three types of specific immunity
    double omicron = std::max(vc_immunity, std::max(st_immunity, sg_immunity));
    
    // susceptibility is function of both exclusion and immunity
    susceptibilities_[s] = (1 - omega) * (1 - omicron); // modified from Cobey 2012, S4.
  }
}

void Host::update_infectiousness() {
  for (int s = 0; s < Serotype::get_num_serotypes(); s++) {
    infectiousness_[s] = 1.0 - find_max_reduction(Vaccine::EffectType::Infectiousness, s);
  }
}

void Host::update_duration_factors() {
  for (int s = 0; s < Serotype::get_num_serotypes(); s++) {
    duration_factors_[s] = 1.0 - find_max_reduction(Vaccine::EffectType::Duration, s);
  }
}

double Host::find_max_reduction(Vaccine::EffectType effect_type, int serotype) const {
  double max_reduction = 0.0;
  for (auto &kv : vaccination_history_) {
    auto& vaccine_name = kv.first;
    auto& num_doses = kv.second;
    if (Vaccine::targets_serotype(vaccine_name, effect_type, Serotype::get_name(serotype))) {   
      max_reduction = std::max(
        max_reduction,
        Vaccine::get_reduction(vaccine_name, effect_type, num_doses)
      );
    }
  }
  return max_reduction;
}

// Static functions
void Host::configure(std::string config_file) {
  auto config_folder = boost::filesystem::path(config_file).parent_path();
  boost::property_tree::ptree ptree = load_config(config_file, "host");
  // strength of competitive exclusion
  mu_ = ptree.get<double>("mu");
  if ( mu_ < 0.0 || mu_ > 1.0 )
    throw std::runtime_error { "mu must be between 0 and 1 (inclusive)." };
  
  // acquisition rate of duration-reducing immunity
  epsilon_ = ptree.get<double>("epsilon");
  if ( epsilon_ < 0.0 )
    throw std::runtime_error { "epsilon cannot be negative " };
  
  // strength of serotype-specific immunity (may vary by serotype)
  auto sigmas_path = (config_folder / ptree.get<std::string>("sigmas_file")).string();
  sigmas_ = load_vector<double>(sigmas_path, "sigmas");
  
//  sigma_ = ptree.get<double>("sigma");
//  if ( sigma_ < 0.0 || sigma_ > 1.0 )
//    throw std::runtime_error { "sigma must be between 0 and 1 (inclusive)." };
  
  // strength of serogroup-specific immunity
  sigma_sg_ = ptree.get<double>("sigma_sg");
  if ( sigma_sg_ < 0.0 || sigma_sg_ > 1.0 )
    throw std::runtime_error { "sigma_sg must be between 0 and 1 (inclusive)." };

  // waning of serotype-specific immunity
  auto halflife_st = ptree.get<double>("halflife_st");
  if ( halflife_st < 0.0 && halflife_st != -1.0 )
    throw std::runtime_error { "halflife_st must be positive, or -1 to indicate no waning." };
  decay_st_ = (halflife_st > 0.0) ? std::log(2) / halflife_st : 0.0;
  
  // waning of serogroup-specific immunity
  auto halflife_sg = ptree.get<double>("halflife_sg");
  if ( halflife_sg < 0.0 && halflife_sg != -1.0 )
    throw std::runtime_error { "halflife_sg must be positive, or -1 to indicate no waning." };
  decay_sg_ = (halflife_sg > 0.0) ? std::log(2) / halflife_sg : 0.0;
}

double Host::get_mu(){
  return Host::mu_;
}
double Host::get_epsilon(){
  return Host::epsilon_;
}
double Host::get_sigma(int serotype){
  return Host::sigmas_[serotype];
}