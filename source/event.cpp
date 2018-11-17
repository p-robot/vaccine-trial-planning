#include <stdexcept>
#include <string>
#include "event.h"
#include "random_number_generator.h"
#include "serotype.h"

RandomNumberGenerator Event::rng_ {0};

Event::Event(int host_id, int day, Type type, bool no_offset)
  : host_id_(host_id)
  , day_(day)
  , type_(type) 
{
  if (type_ != Type::Birthday && type_ != Type::Death) {
    throw std::runtime_error { "Event::Event: Use alternate constructor." };
  }
  day_offset_ = no_offset ? 0 : rng_.draw_uniform_real(0, 1);
}

Event::Event(int host_id, int day, Type type, int serotype)
  : host_id_(host_id)
  , day_(day)
  , type_(type)
  , serotype_(serotype)
{
  if (type_ != Type::Colonization) {
    throw std::runtime_error { "Event::Event: Use alternate constructor." };
  }
  day_offset_ = rng_.draw_uniform_real(0, 1);
}

Event::Event(int host_id, int day, Type type, int serotype, int day_of_colonization)
: host_id_(host_id)
, day_(day)
, type_(type)
, serotype_(serotype)
, day_of_colonization_(day_of_colonization)
{
  if (type_ != Type::Recovery) {
    throw std::runtime_error { "Event::Event: Use alternate constructor." };
  }
  day_offset_ = rng_.draw_uniform_real(0, 1);
}

Event::Event(int host_id, int day, Type type, std::string vaccine_name)
  : host_id_(host_id)
  , day_(day)
  , type_(type) 
  , vaccine_name_(vaccine_name)
{
  if (type_ != Type::Vaccination) {
    throw std::runtime_error { "Event::Event: Use alternate constructor." };
  }
  day_offset_ = rng_.draw_uniform_real(0, 1);
}


int Event::get_host_id() const {
  return host_id_;
}
int Event::get_day() const {
  return day_;
}
Event::Type Event::get_type() const {
  return type_;
}
int Event::get_serotype() const {
  if (type_ != Type::Colonization && type_ != Type::Recovery) {
    throw std::runtime_error { "Event::get_serotype: This event has no associated serotype." };
  }
  return serotype_;
}
std::string Event::get_vaccine_name() const {
  if (type_ != Type::Vaccination)
    throw std::runtime_error { "Event::get_vaccine: This event has no associated vaccine." };
  return vaccine_name_;
}
double Event::get_day_offset() const {
  return day_offset_;
}

std::string Event::as_csv_string(bool abbreviated) const {
  std::ostringstream csv_string;
  csv_string << host_id_ << "," << day_ << ",";
  if (abbreviated) {
    csv_string << Event::to_letter_code(type_);
  } else {
    csv_string << Event::to_string(type_);
  }
  if (type_ == Event::Type::Colonization) {
    csv_string << "," << Serotype::get_name(serotype_) << ",,";
  }
  if (type_ == Event::Type::Recovery) {
    csv_string << "," << Serotype::get_name(serotype_) << "," << day_of_colonization_ << ",";
  }
  if (type_ == Event::Type::Vaccination) {
    csv_string << ",,," << vaccine_name_;
  }
  return csv_string.str();
}

std::string Event::to_string(Event::Type type) {
  switch (type) {
    case Event::Type::Birthday     : return "Birthday";
    case Event::Type::Colonization : return "Colonization";
    case Event::Type::Vaccination  : return "Vaccination";
    case Event::Type::Recovery     : return "Recovery";
    case Event::Type::Death        : return "Death";
  }
}

std::string Event::to_letter_code(Event::Type type) {
  switch (type) {
    case Event::Type::Birthday     : return "B";
    case Event::Type::Colonization : return "C";
    case Event::Type::Vaccination  : return "V";
    case Event::Type::Recovery     : return "R";
    case Event::Type::Death        : return "D";
  }
}

