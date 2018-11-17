#pragma once

#include <string>
#include <vector>

#include "event.h"
#include "population.h"
#include "scribe.h"

class TrialArm {
public:
  TrialArm(std::string name, int num_subjects, int start_day,
           std::string vaccine_name, std::vector<int> schedule);
  std::string get_name() const;
  int get_start_day() const;
  int get_num_subjects() const;
  std::string get_vaccine_name() const;
  std::vector<int> get_schedule() const;
  int get_last_vaccination_day() const;
  
private:
  std::string name_;
  int num_subjects_;
  int start_day_;
  int last_vaccination_day_;
  std::string vaccine_name_;
  std::vector<int> schedule_;
};
