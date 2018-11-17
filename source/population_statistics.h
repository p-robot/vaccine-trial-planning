#pragma once

#include "population.h"

class Population::Statistics {
  
public:
  Statistics(Population& population);
  std::vector<int> sample_hosts(std::pair<int, int> age_range) const; // sample one colonization from each host
  int get_num_hosts() const;
  int get_num_hosts_u5() const;
  int get_num_colonizations(int serotype) const;
  int get_num_colonized() const;
  int get_num_colonized(std::pair<int, int> age_range) const;
  std::vector<double> get_num_colonized_u5() const; // indexed by serotype, hosts contributes at most 1
  
  int get_num_vaccinated(std::string vaccine_name, int num_doses, std::pair<int, int> age_range_in_days) const;
  int get_num_vaccinated(std::string vaccine_name, int num_doses) const;
  int get_num_hosts(std::pair<int, int> age_range_in_days) const;
  
  // new methods for trial arms. they're slow! iterates over all hosts
  int get_num_participants(std::string arm_name) const;
  int get_num_colonized_participants(std::string arm_name) const;
  //
  
  // age-stratified statistics
  std::vector<int> get_num_hosts_by_age() const; // indexed by age (1-year intervals)
  std::vector<int> get_num_colonizations_by_age(int serotype) const; // indexed by age, host may be counted may than once for each serotype
  std::vector<double> get_num_colonized_by_age(int serotype, bool divide_host=false) const; // indexed by age, host contributes at most once to each serotype
  std::vector<int> get_num_previously_colonized_by_age(int serotype) const;
  
  // distributions for an age stratum of the general population
  std::vector<int> get_num_colonizations_dist(std::pair<int, int> age_range, int max_count) const;
  std::vector<int> get_num_past_colonizations_dist(std::pair<int, int> age_range, int max_count) const;
  std::vector<int> get_num_serotypes_dist(std::pair<int, int> age_range, int max_count) const;  // indexed by number of serotypes
  
  // distributions for a trial arm
  std::vector<int> get_num_colonizations_dist(std::string arm_name, int max_count) const;
  std::vector<int> get_num_past_colonizations_dist(std::string arm_name, int max_count) const;
  std::vector<int> get_num_serotypes_dist(std::string arm_name, int max_count) const;  // indexed by number of serotypes
  
private:
  std::vector<int> get_dist_over_hosts(std::string arm_name, std::pair<int, int> age_range, int max_count,
                                       std::function<int(std::shared_ptr<const Host>)> lambda) const;
  const Population& pop_;
};
