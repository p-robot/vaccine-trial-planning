#pragma once
#include <string>
#include "random_number_generator.h"

class Event {
public:
  enum class Type {
    Birthday,
    Colonization,
    Vaccination,
    Recovery,
    Death
  };

  Event(int host_id, int day, Type type, bool no_offset=false);     // for birthday and death events
  Event(int host_id, int day, Type type, int serotype);             // for colonization events
  Event(int host_id, int day, Type type, int serotype, int day_of_colonzation); // for recovery events
  Event(int host_id, int day, Type type, std::string vaccine_name); // for vaccination events

  int get_host_id() const;
  int get_day() const;
  Type get_type() const;
  int get_serotype() const;
  int get_day_of_colonization() const;
  std::string get_vaccine_name() const;
  double get_day_offset() const;
  std::string as_csv_string(bool abbreviated=false) const;
  
private:
  int host_id_;
  int day_;
  Type type_;
  int serotype_;
  int day_of_colonization_;
  std::string vaccine_name_;
  double day_offset_;

  static RandomNumberGenerator rng_;
  static std::string to_string(Event::Type type);
  static std::string to_letter_code(Event::Type type);
};
