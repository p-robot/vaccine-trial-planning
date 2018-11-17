#include "host_container.h"
#include "population.h"
#include "population_statistics.h"



Population::Statistics::Statistics(Population& population)
  : pop_ { population }
{
  
}

std::vector<int> Population::Statistics::sample_hosts(std::pair<int, int> age_range) const {
  auto counts = std::vector<int>(Serotype::get_num_serotypes(), 0);
  for (auto& host_ptr : pop_.hosts_.get<id_tag>()) {
    if (!host_ptr->is_in_trial()) {
      auto age = host_ptr->calculate_age();
      if (age >= age_range.first && age <= age_range.second) {
        if (host_ptr->get_num_colonizations() > 0) {
          auto c = host_ptr->sample_colonization();
          counts[c.get_serotype()]++;
        }
      }
    }
  }
  return counts;
}

int Population::Statistics::get_num_hosts() const {
  return pop_.num_effective_hosts_;
}

int Population::Statistics::get_num_colonized() const {
  int count = 0;
  for (auto& host_ptr : pop_.hosts_.get<id_tag>()) {
    if (!host_ptr->is_in_trial()) {
      count += (int) (host_ptr->get_num_colonizations() > 0);
    }
  }
  return count;
}

int Population::Statistics::get_num_colonized(std::pair<int, int> age_range) const {
  int count = 0;
  for (auto& host_ptr : pop_.hosts_.get<id_tag>()) {
    if (!host_ptr->is_in_trial()) {
      auto age = host_ptr->calculate_age();
      if (age >= age_range.first && age <= age_range.second) {
        count += (int) (host_ptr->get_num_colonizations() > 0);
      }
    }
  }
  return count;
}

std::vector<double> Population::Statistics::get_num_colonized_u5() const {
  return pop_.num_colonized_u5_;
}

int Population::Statistics::get_num_colonizations(int serotype) const {
  int count = 0;
  for (auto& host_ptr : pop_.hosts_.get<id_tag>()) {
    if (!host_ptr->is_in_trial()) {
      count += host_ptr->get_num_colonizations(serotype);
    }
  }
  return count;
}

int Population::Statistics::get_num_vaccinated(std::string vaccine_name, int num_doses, std::pair<int, int> age_range_in_days) const {
  int count = 0;
  for (auto& host_ptr : pop_.hosts_.get<id_tag>()) {
    if (!host_ptr->is_in_trial()) {
      auto age = host_ptr->calculate_age_in_days();
      if (age >= age_range_in_days.first && age <= age_range_in_days.second) {
        count += (int) (host_ptr->get_num_vaccinations(vaccine_name) >= num_doses);
      }
    }
  }
  return count;
}

int Population::Statistics::get_num_vaccinated(std::string vaccine_name, int num_doses) const {
  return get_num_vaccinated(vaccine_name, num_doses, std::make_pair(0, pop_.get_max_age()));
}

int Population::Statistics::get_num_hosts(std::pair<int, int> age_range_in_days) const {
  int count = 0;
  for (auto& host_ptr : pop_.hosts_.get<id_tag>()) {
    if (!host_ptr->is_in_trial()) {
      auto age = host_ptr->calculate_age_in_days();
      if (age >= age_range_in_days.first && age <= age_range_in_days.second) {
        count += 1;
      }
    }
  }
  return count;
}

int Population::Statistics::get_num_hosts_u5() const {
  return pop_.host_ids_u5_.size();
}

// new methods for trial arms. they're slow! iterates over all hosts
int Population::Statistics::get_num_participants(std::string arm_name) const {
  int count = 0;
  for (auto& h : pop_.hosts_) {
    if (h->get_trial_arm_name() == arm_name) {
      count += 1;
    }
  }
  return count;
}

int Population::Statistics::get_num_colonized_participants(std::string arm_name) const {
  int count = 0;
  for (auto& h : pop_.hosts_) {
    if (h->get_trial_arm_name() == arm_name && h->get_num_colonizations() > 0) {
      count += 1;
    }
  }
  return count;
}

// age-stratified
std::vector<int> Population::Statistics::get_num_hosts_by_age() const {
  auto counts = std::vector<int>(pop_.max_age_ + 1, 0);
  for (auto& host_ptr : pop_.hosts_.get<id_tag>()) {
    if (!host_ptr->is_in_trial()) {
      counts[host_ptr->calculate_age()]++;
    }
  }
  return counts;
}

std::vector<int> Population::Statistics::get_num_colonizations_by_age(int serotype) const {
  auto counts = std::vector<int>(pop_.max_age_ + 1, 0);
  for (auto& host_ptr : pop_.hosts_.get<id_tag>()) {
    if (!host_ptr->is_in_trial()) {
      counts[host_ptr->calculate_age()] += host_ptr->get_num_colonizations(serotype);
    }
  }
  return counts;
}

std::vector<double> Population::Statistics::get_num_colonized_by_age(int serotype, bool divide_host) const {
  auto counts = std::vector<double>(pop_.max_age_ + 1, 0);
  for (auto& host_ptr : pop_.hosts_.get<id_tag>()) {
    if (!host_ptr->is_in_trial() && host_ptr->get_num_colonizations(serotype) > 0) {
      if (divide_host) {
        double frac = host_ptr->get_num_colonizations(serotype) / (double) host_ptr->get_num_colonizations();
        counts[host_ptr->calculate_age()] += frac;
      } else {
        counts[host_ptr->calculate_age()] += 1;
      }
    }
  }
  return counts;
}

std::vector<int> Population::Statistics::get_num_previously_colonized_by_age(int serotype) const {
  auto counts = std::vector<int>(pop_.max_age_ + 1, 0);
  for (auto& host_ptr : pop_.hosts_.get<id_tag>()) {
    if (!host_ptr->is_in_trial() && host_ptr->get_num_past_colonizations(serotype) > 0) {
      counts[host_ptr->calculate_age()] += 1;
    }
  }
  return counts;
}

// helper function for all the distribution functions
std::vector<int> Population::Statistics::get_dist_over_hosts(std::string arm_name, std::pair<int, int> age_range, int max_count,
                                                 std::function<int(std::shared_ptr<const Host>)> host_func) const {
  auto dist = std::vector<int>(max_count + 2, 0);
  for (auto& host_ptr : pop_.hosts_.get<id_tag>()) {
    if (host_ptr->get_trial_arm_name() == arm_name) {
      int age = host_ptr->calculate_age();
      if (age >= age_range.first && age <= age_range.second) {
        int n = host_func(host_ptr);
        int bin = (n <= max_count) ? n : max_count + 1; // determine bin
        dist[bin]++;
      }
    }
  }
  return dist;
}

// distributions for an age stratum of the general population
std::vector<int> Population::Statistics::get_num_colonizations_dist(std::pair<int, int> age_range, int max_count) const {
  return get_dist_over_hosts("", age_range, max_count,
                             [](std::shared_ptr<const Host> h) { return h->get_num_colonizations(); });
}

std::vector<int> Population::Statistics::get_num_past_colonizations_dist(std::pair<int, int> age_range, int max_count) const {
  return get_dist_over_hosts("", age_range, max_count,
                             [](std::shared_ptr<const Host> h) { return h->get_num_past_colonizations(); });
}
std::vector<int> Population::Statistics::get_num_serotypes_dist(std::pair<int, int> age_range, int max_count) const {
  return get_dist_over_hosts("", age_range, max_count,
                             [](std::shared_ptr<const Host> h) { return h->get_num_serotypes(); });
}


// distributions for a trial arm
std::vector<int> Population::Statistics::get_num_colonizations_dist(std::string arm_name, int max_count) const {
  return get_dist_over_hosts(arm_name, std::make_pair(0, pop_.max_age_), max_count,
                             [](std::shared_ptr<const Host> h) { return h->get_num_colonizations(); });
}
std::vector<int> Population::Statistics::get_num_past_colonizations_dist(std::string arm_name, int max_count) const {
  return get_dist_over_hosts(arm_name, std::make_pair(0, pop_.max_age_), max_count,
                             [](std::shared_ptr<const Host> h) { return h->get_num_past_colonizations(); });
}
std::vector<int> Population::Statistics::get_num_serotypes_dist(std::string arm_name, int max_count) const {
  return get_dist_over_hosts(arm_name, std::make_pair(0, pop_.max_age_), max_count,
                             [](std::shared_ptr<const Host> h) { return h->get_num_serotypes(); });
}
















