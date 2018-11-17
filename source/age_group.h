#pragma once
#include <string>
#include <unordered_set>
#include <vector>

#include "host.h"
#include "random_number_generator.h"

class AgeGroup {
public:
  AgeGroup(std::pair<int, int> age_range, int num_serotypes,
           std::shared_ptr<RandomNumberGenerator> rng_ptr);
  
  std::pair<int, int> get_age_range() const;
  bool contains_age(int age) const;
  bool contains_host_id(int host_id) const;
  
  int size() const;           // all hosts
  int effective_size() const; // only hosts *not* in a trial arm
  double get_num_effective_colonizations(int serotype) const; // summed over hosts *not* in a trial arm
  
  const std::vector<int> get_host_ids() const;
  int draw_random_host_id() const;

  void add_host(std::shared_ptr<const Host> host_ptr);
  void remove_host(std::shared_ptr<const Host> host_ptr);
  void update(std::shared_ptr<const Host> host_ptr, int serotype, std::string event);
  
  template<class InputIterator>
  static InputIterator find_by_age(InputIterator first, InputIterator last, int age) {
    while (first != last) {
      if (first->contains_age(age)) 
        return first;
      ++first;
    }
    return last;
  }

  template<class InputIterator>
  static InputIterator find_by_host_id(InputIterator first, InputIterator last, int host_id) {
    while (first != last) {
      if (first->contains_host_id(host_id)) 
        return first;
      ++first;
    }
    return last;
  }

private:
  std::vector<int> host_ids_;
  std::pair<int, int> age_range_; // in years, [min, max)
  std::shared_ptr<RandomNumberGenerator> rng_ptr_;
  int effective_size_;
  std::vector<double> num_effective_colonizations_;
};