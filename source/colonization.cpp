#include <memory>
#include <stdexcept>
#include <string>
#include "colonization.h"
#include "random_number_generator.h"

Colonization::Colonization(int serotype, int day_of_colonization, int day_of_recovery,
                           std::shared_ptr<RandomNumberGenerator> rng_ptr) 
    : serotype_            { serotype }
    , day_of_colonization_ { day_of_colonization }
    , day_of_recovery_     { day_of_recovery }
    , rng_ptr_             { rng_ptr }
{
  if (day_of_recovery_ < day_of_colonization_) {
    throw std::runtime_error { "Colonization::Colonization: Day of recovery cannot be less than day of colonization." };
  }
}     

int Colonization::get_serotype() const {
  return serotype_;
}
int Colonization::get_day_of_colonization() const {
  return day_of_colonization_;
}
int Colonization::get_day_of_recovery() const {
  return day_of_recovery_;
}
int Colonization::get_duration() const {
  return day_of_recovery_ - day_of_colonization_;
}






