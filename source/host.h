#pragma once
#include <memory>
#include <vector>
#include <string>
#include <set>
#include <unordered_map>
#include "boost/container/flat_set.hpp"

#include "calendar.h"
#include "colonization.h"
#include "serotype.h"
#include "random_number_generator.h"
#include "vaccine.h"

class Host {
public:
  Host(int day_of_birth, int day_of_death,
       std::shared_ptr<const Calendar> calendar_ptr,
       std::shared_ptr<RandomNumberGenerator> rng_ptr,
       std::string trial_arm_name="");

  // Getters and setters
  int get_id() const;
  int get_day_of_birth() const;
  int get_day_of_death() const;
  int get_lifespan() const; // in years

  // Demographic member functions
  int calculate_age() const; // in years
  int calculate_age_in_days() const;

  // Colonization-related member functions
  int get_num_colonizations(int serotype) const;
  int get_num_colonizations() const;
  int get_num_past_colonizations(int serotype) const;
  int get_num_past_colonizations() const;
  int get_num_serotypes() const;
  void verify_num_colonizations() const;

  double get_susceptibility(int serotype) const;
  double get_infectiousness(int serotype) const;
  double get_duration_factor(int serotype) const;
  
  void acquire_past_colonizations(int serotype, int n);
  Colonization sample_colonization() const; // a random sample of the hosts's colonizations
  Colonization acquire_colonization(int serotype);
  std::vector<int> remove_expired_colonizations();  // returns serotypes removed

  // Vaccination-related
  void acquire_vaccination(std::string vaccine_name);
  int get_num_vaccinations(std::string vaccine_name) const;

  // Vaccine trial-related
  std::string get_trial_arm_name() const;
  bool is_in_trial() const;
  
  // Configuration
  static void configure(std::string config_file);
  static double get_mu();
  static double get_epsilon();
  static double get_sigma(int serotype);
  static const int null_id = -1;

private:
  int calculate_colonization_duration(int serotype); 
  void update_susceptibilities(); // run every time there's a change in carriage or vaccination status
  void update_infectiousness();
  void update_duration_factors();
  double find_max_reduction(Vaccine::EffectType effect_type, int serotype) const;

  // demographic details
  int id_;
  int day_of_birth_;
  int day_of_death_;
  std::shared_ptr<const Calendar> calendar_ptr_;
  std::shared_ptr<RandomNumberGenerator> rng_ptr_;

  // colonization details
  std::vector<Colonization> colonizations_;
  std::vector<int> num_colonizations_;
  std::vector<int> num_past_colonizations_;
  std::vector<bool> serotypes_cleared_;
  std::vector<int> serotype_clearance_dates_;
  boost::container::flat_set<int> serogroups_cleared_;
  std::unordered_map<int, int> serogroup_clearance_dates_; // based on most recently cleared colonization (for a given serogroup)
  int total_past_colonizations_;
  std::vector<double> susceptibilities_;
  std::vector<double> infectiousness_;
  std::vector<double> duration_factors_;

  // vaccination history
  std::unordered_map<std::string, int> vaccination_history_;
  
  // vaccine trial-related
  std::string trial_arm_name_;
  bool in_trial_;
  
  // static members
  static int next_unused_id; 
  static double mu_;       // scales the resistance to acquisition conferred by most fit std::string
  static double epsilon_;  // non-serotype-specific, colonization duration-reducing immunity
//  static double sigma_;    // anti-capsular immmunity (serotype-specific immunity)
  static std::vector<double> sigmas_;
  static double sigma_sg_; // serogroup-specific immunity
  static double decay_st_; // scales exponential decay of serotype-specific immunity
  static double decay_sg_; // scales exponential decay of serogroup-specific immunity

};