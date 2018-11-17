#include <unordered_set>
#include <vector>

#include "age_group.h"
#include "host.h"
#include "random_number_generator.h"

AgeGroup::AgeGroup(std::pair<int, int> age_range, int num_serotypes,
                   std::shared_ptr<RandomNumberGenerator> rng_ptr)
  : age_range_                   { age_range }
  , effective_size_              { 0 }
  , num_effective_colonizations_ { std::vector<double>(num_serotypes, 0.0) }
  , rng_ptr_                     { rng_ptr }
{}

std::pair<int, int> AgeGroup::get_age_range() const {
  return age_range_;
}
bool AgeGroup::contains_age(int age) const {
  return (age >= age_range_.first) && (age < age_range_.second);
}
bool AgeGroup::contains_host_id(int host_id) const {
  auto it = std::find(host_ids_.begin(), host_ids_.end(), host_id);
  return it != host_ids_.end();
}

int AgeGroup::size() const {
  return host_ids_.size();
}
int AgeGroup::effective_size() const {
  return effective_size_;
}
double AgeGroup::get_num_effective_colonizations(int serotype) const {
  return num_effective_colonizations_.at(serotype);
}

const std::vector<int> AgeGroup::get_host_ids() const {
  return host_ids_;
}
int AgeGroup::draw_random_host_id() const {
  int r = rng_ptr_->draw_uniform_int(0, host_ids_.size() - 1);
  return host_ids_[r];
}


void AgeGroup::add_host(std::shared_ptr<const Host> host_ptr) {
  host_ids_.push_back(host_ptr->get_id());
  if (!host_ptr->is_in_trial()) {
    effective_size_ += 1;
    for (int s = 0; s < num_effective_colonizations_.size(); s++) {
      double c = host_ptr->get_infectiousness(s) * host_ptr->get_num_colonizations(s);
      num_effective_colonizations_[s] += c;
    }
  }
}
void AgeGroup::remove_host(std::shared_ptr<const Host> host_ptr) {
  auto it = std::find(host_ids_.begin(), host_ids_.end(), host_ptr->get_id());
  if (it != host_ids_.end()) {
    host_ids_.erase(it);
  }
  if (!host_ptr->is_in_trial()) {
    effective_size_ -= 1;
    for (int s = 0; s < num_effective_colonizations_.size(); s++) {
      double c = host_ptr->get_infectiousness(s) * host_ptr->get_num_colonizations(s);
      num_effective_colonizations_[s] -= c;
    }
  }
}
void AgeGroup::update(std::shared_ptr<const Host> host_ptr, int serotype, std::string event) {
  if (!host_ptr->is_in_trial()) {
    double delta = host_ptr->get_infectiousness(serotype);
    if (event == "colonized") {
      num_effective_colonizations_.at(serotype) += delta;
    } else if (event == "recovered") {
      num_effective_colonizations_.at(serotype) -= delta;
    } else {
      throw std::runtime_error("AgeGroup::update_num_effective_colonizations: Event '"
                               + event + "' not recognized.");
    }
  }
}

