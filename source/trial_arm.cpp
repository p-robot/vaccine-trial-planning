#include <algorithm>
#include <string>
#include <vector>

#include "trial_arm.h"

TrialArm::TrialArm(std::string name, int num_subjects, int start_day,
         std::string vaccine_name, std::vector<int> schedule)
  : name_                 { name }
  , num_subjects_         { num_subjects }
  , start_day_            { start_day }
  , vaccine_name_         { vaccine_name }
  , schedule_             { schedule }
  , last_vaccination_day_ { *std::max_element (schedule.begin(), schedule.end())}
{

}

std::string TrialArm::get_name() const {
  return name_;
}
int TrialArm::get_start_day() const {
  return start_day_;
}
int TrialArm::get_num_subjects() const {
  return num_subjects_;
}
std::string TrialArm::get_vaccine_name() const {
  return vaccine_name_;
}
std::vector<int> TrialArm::get_schedule() const {
  return schedule_;
}

int TrialArm::get_last_vaccination_day() const {
  return last_vaccination_day_;
}
