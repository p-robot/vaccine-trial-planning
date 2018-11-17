#pragma once
#include <memory>
#include <string>
#include "random_number_generator.h"

class Colonization {
public:
  Colonization(int serotype, int day_of_colonization, int day_of_recovery,
               std::shared_ptr<RandomNumberGenerator> rng_ptr);
   
  // Self-explanatory getters
  int get_serotype() const;
  int get_day_of_colonization() const;
  int get_day_of_recovery() const;
  int get_duration() const;

private:
  int serotype_; 
  int day_of_colonization_;
  int day_of_recovery_;  
  std::shared_ptr<RandomNumberGenerator> rng_ptr_; 
};
