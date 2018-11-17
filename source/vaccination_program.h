#pragma once
#include <string>
#include <vector>

class VaccinationProgram {
public:
  VaccinationProgram(std::string vaccine_name, int start_day, int end_day,
                     std::vector<int> dose_schedule,
                     std::vector<double> annual_coverage);

  std::string get_vaccine_name() const;
  int get_start_day() const;
  int get_end_day() const;
  double get_coverage(int day) const;
  std::vector<int> get_schedule() const;

private:
  std::string vaccine_name_;
  int start_day_;
  int end_day_;
  std::vector<double> annual_coverage_;
  std::vector<int> dose_schedule_;
};