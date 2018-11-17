#include <stdexcept>
#include <vector>
#include "incidence_logger.h"

IncidenceLogger::IncidenceLogger(std::set<std::string> groups, int num_serotypes, std::string output_folder)
  : groups_        { groups }
  , num_serotypes_ { num_serotypes }
  , scribe_        { output_folder }
{
  for (auto& g : groups) {
    incidents_[g] = std::vector<int>(num_serotypes_, 0);
  }
}
  
void IncidenceLogger::log(std::string group, int serotype) {
  if (groups_.count(group) == 0)
    throw std::runtime_error { "IncidenceLogger::log: '" + group + "' was not specified as a group." };
  if (serotype >= num_serotypes_ || serotype < 0)
    throw std::runtime_error { "IncidenceLogger::log: Invalid serotype index: " + std::to_string(serotype) + "." };

  incidents_[group].at(serotype)++;
}

void IncidenceLogger::flush(int current_day) {
  scribe_.note_value(current_day, "incidence_samplings_days.csv");
  for (auto& g : groups_) {
    scribe_.note_values(incidents_[g], "incidence_" + g + ".csv");
    incidents_[g] = std::vector<int>(num_serotypes_, 0); // reset counts
  }
}