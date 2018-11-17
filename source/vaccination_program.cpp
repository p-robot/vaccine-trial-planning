#include <string>
#include <vector>
#include "vaccination_program.h"

VaccinationProgram::VaccinationProgram(std::string vaccine_name, int start_day, int end_day,
                                       std::vector<int> dose_schedule,
                                       std::vector<double> annual_coverage_)
  : vaccine_name_ { vaccine_name }
  , start_day_ { start_day }
  , end_day_ { end_day }
  , dose_schedule_ { dose_schedule }
  , annual_coverage_ { annual_coverage_ }
{}

std::string VaccinationProgram::get_vaccine_name() const {
  return vaccine_name_;
}
int VaccinationProgram::get_start_day() const {
  return start_day_;
}
int VaccinationProgram::get_end_day() const {
  return end_day_;
}
double VaccinationProgram::get_coverage(int day) const {
  if (day < start_day_ || day > end_day_) {
    return 0.0; // program is not in effect
  } else {
    int year = (day - start_day_) / 365;
    if (year < annual_coverage_.size()) {
      return annual_coverage_[year];
    } else {
      return annual_coverage_.back();
    }
  }
}
std::vector<int> VaccinationProgram::get_schedule() const {
  return dose_schedule_;
}


