#include <cmath>
#include <functional>
#include <memory>
#include <unordered_map>
#include <stdexcept>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

#include "age_group.h"
#include "calendar.h"
#include "colonization.h"
#include "configuration.h"
#include "host.h"
#include "host_container.h"
#include "population.h"
#include "population_statistics.h"
#include "random_number_generator.h"
#include "serotype.h"

namespace {
  double weighted_sum(std::vector<double> weights, std::vector<double> values) {
    double sum = 0.0;
    for (int i = 0; i < weights.size(); i++) {
      sum += (weights[i] * values[i]);
    }
    return sum;
  }
}

Population::Population(std::shared_ptr<const Calendar> calendar_ptr, 
                       std::shared_ptr<RandomNumberGenerator> rng_ptr,
                       std::string config_path)
  : calendar_ptr_     { calendar_ptr }
  , rng_ptr_          { rng_ptr }
  , num_colonized_u5_ { std::vector<double>(Serotype::get_num_serotypes(), 0.0) }
  , host_ids_u5_      { }
  , initialized_      { false }
  , stats             { std::make_shared<Statistics>(*this) }
{
  configure(config_path);
}

bool Population::has_host(int host_id) const {
  auto& index = hosts_.get<id_tag>();
  return index.find(host_id) != index.end();
}
std::shared_ptr<const Host> Population::find_host(int host_id) const {
  auto itr = hosts_.get<id_tag>().find(host_id);
  return *itr;
}

std::vector<int> Population::initialize(int num_hosts) {
  if (initialized_) {
    throw std::runtime_error { "Population::initialize: Population already initialized." };
  }
  std::vector<int> host_ids;
  for (int i = 0; i < num_hosts; i++) {
    int lifespan = choose_random_day(lifespan_pmf_);
    int age_in_days = rng_ptr_->draw_uniform_int(0, lifespan);
    auto host_ptr = std::make_shared<Host>(-age_in_days, -age_in_days + lifespan, calendar_ptr_, rng_ptr_);
    hosts_.insert(host_ptr);
    host_ids.push_back(host_ptr->get_id());
    
    auto group_itr = AgeGroup::find_by_age(age_groups_.begin(), age_groups_.end(), host_ptr->calculate_age());
    group_itr->add_host(host_ptr);
  }
  num_effective_hosts_ = num_hosts;
  initialized_ = true;
  return host_ids;
}

std::vector<std::pair<int, int>> Population::seed_colonizations() {
  std::vector<std::pair<int, int>> host_serotype_pairs;
  auto& index = hosts_.get<id_tag>();
  for (auto& host_ptr : index) {
    for (int s = 0; s < Serotype::get_num_serotypes(); s++) {
      if (rng_ptr_->draw_bernoulli(initial_colonization_prob_)) {
        host_serotype_pairs.push_back( std::make_pair(host_ptr->get_id(), s) );
      }
    }
  }
  return host_serotype_pairs;
}

std::vector<std::pair<int, int>> Population::determine_colonizations() {
  std::vector<std::pair<int, int>> host_serotype_pairs;
  auto num_hosts = num_effective_hosts_;
  for (int s = 0; s < Serotype::get_num_serotypes(); s++) {
    // get a vector of number of effective colonizations, indexed by age group
    std::vector<double> colonizations_by_group;
    for (auto& g: age_groups_) {
      colonizations_by_group.push_back(g.get_num_effective_colonizations(s));
    }

    // determine colonizations in each age group
    auto& index = hosts_.get<id_tag>();
    for (int i = 0; i < age_groups_.size(); i++) {
      // if no one is in the group, skip it
      auto group = age_groups_[i];
      if (group.size() == 0) {
        continue;
      }

      // calculate contact rate (for entire age group) based on weighted contributions from other age groups
      auto summation = weighted_sum(mixing_matrix_[i], colonizations_by_group);
      double contact_rate = (beta_ * summation / num_hosts + omega_) * group.size();

      // draw random hosts from this age group for possible colonization
      int num_contacts = rng_ptr_->draw_poisson(contact_rate);
      for (int j = 0; j < num_contacts; j++) {
        auto host_id = group.draw_random_host_id();
        auto host_itr = index.find(host_id);
        if (host_itr != index.end()) {
          auto host_ptr = *(host_itr);
          if (rng_ptr_->draw_bernoulli(host_ptr->get_susceptibility(s))) {
            host_serotype_pairs.push_back( std::make_pair(host_id, s) );
          }
        }
      }
    }
  }
  return host_serotype_pairs;
}

void Population::seed_past_colonizations(int num_past_colonizations) {
  auto& index = hosts_.get<id_tag>();
  for (auto itr = index.begin(); itr != index.end(); itr++) {
    for (int s = 0; s < Serotype::get_num_serotypes(); s++ ) {
      index.modify(itr, [&](std::shared_ptr<Host>& h) {
        h->acquire_past_colonizations(s, num_past_colonizations);
      });
    }
  }
}
void Population::verify() const {
  // check number of hosts (in the general population) is correct
  int genpop_count = 0;
  for (auto& host_ptr : hosts_.get<id_tag>()) {
    genpop_count += (int) (!host_ptr->is_in_trial());
  }
  if (genpop_count != num_effective_hosts_) {
    throw std::runtime_error { "Population::verify: Number of effective hosts may be incorrect." };
  }
  
  // assign each host to an age group
  std::unordered_map<int, std::size_t> host2group;
  for (std::size_t i = 0; i < age_groups_.size(); ++i) {
    for (auto& host_id : age_groups_[i].get_host_ids()) {
      host2group[host_id] = i;
    }
  }
  
  // check that we have accounted for all host_ids
  for (auto& h : hosts_.get<id_tag>()) {
    if (host2group.count(h->get_id()) == 0) {
      throw std::runtime_error { "Population::verify: Host not found in any age group." };
    }
  }
  
  // count all the colonizations we have (for hosts not in a trial arm)
  auto counts = std::vector<std::vector<double>>(age_groups_.size(), std::vector<double>(Serotype::get_num_serotypes(), 0.0));
  for (auto& h : hosts_.get<id_tag>()) {
    h->verify_num_colonizations();
    if (!h->is_in_trial()) {
      for (int s = 0; s < Serotype::get_num_serotypes(); s++) {
        auto c = h->get_num_colonizations(s) * h->get_infectiousness(s);
        counts[host2group[h->get_id()]][s] += c;
      }
    }
  }

  // compare to what is recorded in the age groups
  for (std::size_t i = 0; i < age_groups_.size(); i++) {
    for (int s = 0; s < Serotype::get_num_serotypes(); s++) {
      auto x1 = age_groups_[i].get_num_effective_colonizations(s);
      auto x2 = counts[i][s];
      if (std::fabs(x1 - x2) > 1e-5) {
        throw std::runtime_error { "Population::verify_num_colonizations: Number of effective colonizations may be incorrect." };
      }
    }
  }
  
  // count the number of under-5 colonized by serotype (only those not in a trial arm)
  auto counts_u5 = std::vector<double>(Serotype::get_num_serotypes(), 0.0);
  for (auto& h : hosts_.get<id_tag>()) {
    auto total = (double) h->get_num_colonizations();
    if (total > 0 && (h->calculate_age() < 5) && !h->is_in_trial()) {
      for (int s = 0; s < Serotype::get_num_serotypes(); s++) {
        double frac = h->get_num_colonizations(s) / total;
        counts_u5[s] += frac;
      }
    }
  }
  
  // compare to our vector we've been using to track number of under-5 colonized
  for (int s = 0; s < Serotype::get_num_serotypes(); s++) {
    auto x1 = num_colonized_u5_[s];
    auto x2 = counts_u5[s];
    if (std::fabs(x1 - x2) > 1e-5) {
      throw std::runtime_error { "Population::verify_num_colonizations: Number of under-5 colonized may be incorrect." };
    }
  }
  
  // check that we've got the right number of under-5's
  int count_u5 = 0;
  for (auto& h : hosts_.get<id_tag>()) {
    if (!h->is_in_trial()) {
      count_u5 += (int) (h->calculate_age() < 5);
    }
  }
  if (count_u5 != (int) host_ids_u5_.size()) {
    throw std::runtime_error { "Population::verify_num_colonizations: Number of under-5 hosts may be incorrect." };
  }
    
}

// Event handlers
int Population::birth_host(std::string trial_arm_name) {
  int day_of_birth = calendar_ptr_->get_current_day();
  int day_of_death = day_of_birth + choose_random_day(lifespan_pmf_);
  auto host_ptr =  std::make_shared<Host>(day_of_birth, day_of_death, calendar_ptr_, rng_ptr_, trial_arm_name);
  hosts_.insert(host_ptr);
  age_groups_[0].add_host(host_ptr);
  add_num_colonized_u5(host_ptr);
  if (!host_ptr->is_in_trial()) {
    host_ids_u5_.insert(host_ptr->get_id());
    num_effective_hosts_++;
  }
  return host_ptr->get_id();
}

void Population::remove_host(int host_id) {
  auto& index = hosts_.get<id_tag>();
  auto host_itr = index.find(host_id);
  if (host_itr != index.end()) {
    auto host_ptr = *(host_itr);
    auto group_itr = AgeGroup::find_by_age(age_groups_.begin(), age_groups_.end(), host_ptr->calculate_age());
    remove_num_colonized_u5(host_ptr);
    if (!host_ptr->is_in_trial()) {
      if (host_ptr->calculate_age() < 5) {
        host_ids_u5_.erase(host_ptr->get_id());
      }
      num_effective_hosts_--;
    }
    group_itr->remove_host(host_ptr);
    index.erase(host_id);
  }
}
int Population::colonize_host(int host_id, int serotype) {
  int day_of_recovery;
  auto host_itr = hosts_.get<id_tag>().find(host_id);
  if (host_itr != hosts_.get<id_tag>().end()) {
    hosts_.modify(host_itr, [&](std::shared_ptr<Host>& h) { 
      remove_num_colonized_u5(h);
      
      Colonization c = h->acquire_colonization(serotype);
      day_of_recovery = c.get_day_of_recovery();
      
      auto s = c.get_serotype();
      auto group_itr = AgeGroup::find_by_age(age_groups_.begin(), age_groups_.end(), h->calculate_age());
      group_itr->update(h, s, "colonized");
      add_num_colonized_u5(h);
    });
  }
  return day_of_recovery;
}

void Population::recover_host(int host_id) {
  auto host_itr = hosts_.get<id_tag>().find(host_id);
  if (host_itr != hosts_.get<id_tag>().end()) {
    hosts_.modify(host_itr, [=](std::shared_ptr<Host>& h) { 
      remove_num_colonized_u5(h);
      
      auto serotypes_removed = h->remove_expired_colonizations();
      for (auto& s : serotypes_removed) {
        auto group_itr = AgeGroup::find_by_age(age_groups_.begin(), age_groups_.end(), h->calculate_age());
        group_itr->update(h, s, "recovered");
      }
      add_num_colonized_u5(h);
    });
  }
}

void Population::vaccinate_host(int host_id, std::string vaccine_name) {
  auto host_itr = hosts_.get<id_tag>().find(host_id);
  if (host_itr != hosts_.get<id_tag>().end()) {
    hosts_.modify(host_itr, [=](std::shared_ptr<Host>& h) { 
      auto group_itr = AgeGroup::find_by_age(age_groups_.begin(), age_groups_.end(), h->calculate_age());
      // remove host and add back host after vaccination,
      // effectively updating the number of colonizations in the age group
      group_itr->remove_host(h);
      h->acquire_vaccination(vaccine_name);
      group_itr->add_host(h);
    });
  }
}

void Population::update_host_age_group(int host_id) {
  auto host_itr = hosts_.get<id_tag>().find(host_id);
  if (host_itr != hosts_.get<id_tag>().end()) {
    auto host_ptr = *host_itr;
    
    // remove from current group and add to new group
    auto current_group_itr = AgeGroup::find_by_host_id(age_groups_.begin(), age_groups_.end(), host_ptr->get_id());
    auto new_group_itr = AgeGroup::find_by_age(age_groups_.begin(), age_groups_.end(), host_ptr->calculate_age());
    
    current_group_itr->remove_host(host_ptr);
    new_group_itr->add_host(host_ptr);
    
    if (host_ptr->calculate_age() == 5) { // host just turned 5...
      remove_num_colonized_u5(host_ptr, false); // remove its contribution, false = no age checking
      if (!host_ptr->is_in_trial()) {
        host_ids_u5_.erase(host_ptr->get_id());
      }
    }
  }
}

// Statistics


//  Helpers
int Population::choose_random_day(std::vector<double> age_dist) {
  return 365 * rng_ptr_->draw_categorical(age_dist) + rng_ptr_->draw_uniform_int(0, 364);
}

void Population::add_num_colonized_u5(std::shared_ptr<Host> host_ptr) {
  auto total = (double) host_ptr->get_num_colonizations();
  if (total > 0 && (host_ptr->calculate_age() < 5) && !host_ptr->is_in_trial()) {
    for (int s = 0; s < Serotype::get_num_serotypes(); s++) {
      double frac = host_ptr->get_num_colonizations(s) / total;
      num_colonized_u5_[s] += frac;
    }
  }
}

void Population::remove_num_colonized_u5(std::shared_ptr<Host> host_ptr, bool check_age) {
  auto total = (double) host_ptr->get_num_colonizations();
  if (total > 0 && !(check_age && (host_ptr->calculate_age() >= 5)) && !host_ptr->is_in_trial()) {
    for (int s = 0; s < Serotype::get_num_serotypes(); s++) {
      double frac = host_ptr->get_num_colonizations(s) / total;
      num_colonized_u5_[s] -= frac;
    }
  }
}



void Population::configure(std::string config_path) {
  namespace pt = boost::property_tree;
  namespace fs = boost::filesystem;

  pt::ptree ptree = load_config(config_path, "population");
  // real-valued parameters
  max_age_                   = ptree.get<int>("max_age");
  beta_                      = ptree.get<double>("beta");
  omega_                     = ptree.get<double>("omega");
  initial_colonization_prob_ = ptree.get<double>("initial_colonization_prob");

  // vector-valued parameters
  auto config_folder = fs::path(config_path).parent_path();
  auto lifespan_pmf_file = (config_folder / ptree.get<std::string>("lifespan_pmf_file")).string();
  lifespan_pmf_ = load_pmf(lifespan_pmf_file, max_age_);

  // age groups
  auto age_mixing_file = (config_folder / ptree.get<std::string>("age_mixing_file")).string();
  auto age_cutoffs = load_vector<int>(age_mixing_file, "age_cutoffs");
  age_cutoffs.insert(begin(age_cutoffs), 0);
  age_cutoffs.push_back(max_age_);
  for (int i = 0; i < age_cutoffs.size() - 1; i++) {
    if (age_cutoffs[i] >= age_cutoffs[i + 1]) {
      throw std::runtime_error {"Population::configure: age_cutoffs must be strictly increasing and between 0 and the maximum age"};
    }
    int min_age = age_cutoffs[i];
    int max_age = age_cutoffs[i + 1];
    age_groups_.emplace_back(std::make_pair(min_age, max_age), Serotype::get_num_serotypes(), rng_ptr_);
  }

  // mixing matrix
  auto mixing_mat_pt = load_config(age_mixing_file, "mixing_matrix");
  for (auto& row : mixing_mat_pt) {
    std::vector<double> row_rates;
    for (auto& col : row.second) {
      row_rates.push_back(col.second.get_value<double>());
    }
    if (row_rates.size() != age_groups_.size()) {
      throw std::runtime_error {"Population::configure: Mixing matrix does not have incorrect shape."};
    }
    mixing_matrix_.push_back(row_rates);
  }
  if (mixing_matrix_.size() != age_groups_.size()) {
    throw std::runtime_error {"Population::configure: Mixing matrix has incorrect shape."};
  }
}

std::vector<double> Population::load_pmf(std::string config_file, int expected_size, double epsilon) {
  auto pmf = load_vector<double>(config_file, "pmf");
  double sum = std::accumulate(pmf.begin(), pmf.end(), 0.0);
  if (std::abs(sum - 1.0) > epsilon) {
    throw std::runtime_error { "Population::load_pmf: Probability mass function does not integrate to 1, but to "
      + std::to_string(sum) + " instead." };
  }
  return pmf;
}


// get configuration parameters
int Population::get_max_age() const {
  return max_age_;
}
std::vector<double> Population::get_lifespan_pmf() const {
  return lifespan_pmf_;
}

double Population::get_beta() const {
  return beta_;
}
double Population::get_omega() const {
  return omega_;
}
double Population::get_initial_colonization_prob() const {
  return initial_colonization_prob_;
}

// set beta for fitting
void Population::set_beta(double beta) {
  beta_ = beta;
}
