#pragma once
#include <functional>
#include <limits>
#include <memory>
#include <set>
#include <vector>
#include <boost/container/flat_set.hpp>

#include "age_group.h"
#include "calendar.h"
#include "host.h"
#include "host_container.h"
#include "random_number_generator.h"

class Population {
  class Statistics;

public:
  Population(std::shared_ptr<const Calendar> calendar_ptr, 
             std::shared_ptr<RandomNumberGenerator> rng_ptr,
             std::string config_path);

  bool has_host(int host_id) const;
  std::shared_ptr<const Host> find_host(int host_id) const;

  // needed to run a simulation
  std::vector<int> initialize(int num_hosts); // create initial population, returns host ids
  std::vector<std::pair<int, int>> seed_colonizations(); // return host_id, serotype pairs
  std::vector<std::pair<int, int>> determine_colonizations();
  void seed_past_colonizations(int num_past_colonizations); // number of past colonizations per host
  void verify() const;

  // events handlers
  int birth_host(std::string trial_arm_name=""); // returns newborn host id (needed for event scheduling)
  void remove_host(int host_id);
  int colonize_host(int host_id, int serotype); // returns day of recovery
  void recover_host(int host_id);
  void vaccinate_host(int host_id, std::string vaccine_name);
  void update_host_age_group(int host_id);
  
  // configuration
  int get_max_age() const;
  std::vector<double> get_lifespan_pmf() const;
  double get_beta() const;
  double get_omega() const;
  double get_initial_colonization_prob() const;

  // for fitting
  void set_beta(double beta);
  
  // statistics
  const std::shared_ptr<const Statistics> stats;

private:
  void configure(std::string config_path);
  int choose_random_day(std::vector<double> age_dist); // choose random age in days

  bool initialized_;
  
  HostMultiIndexContainer hosts_;
  std::shared_ptr<const Calendar> calendar_ptr_;
  std::shared_ptr<RandomNumberGenerator> rng_ptr_;
  
  // caching statistics for faster FOI calculations
  void add_num_colonized_u5(std::shared_ptr<Host> host_ptr);
  void remove_num_colonized_u5(std::shared_ptr<Host> host_ptr, bool check_age=true);
  std::vector<double> num_colonized_u5_;
  boost::container::flat_set<int> host_ids_u5_;
  int num_effective_hosts_; // i.e. all ages, not in a trial
  
  // configurable parameters
  int max_age_;
  std::vector<double> lifespan_pmf_;
  std::vector<std::vector<double>> mixing_matrix_;
  std::vector<AgeGroup> age_groups_;
  
  double beta_;  // contact rate shared by all serotypes
  double omega_; // low immigration rate of each serotype (shared by all serotypes?)
  double initial_colonization_prob_; // shared by all serotypes

  static std::vector<double> load_pmf(std::string config_file, int expected_size=-1, double epsilon=0.001);
};
