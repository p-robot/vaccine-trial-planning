#pragma once
#include <unordered_map>
#include <set>
#include <vector>
#include "population.h"
#include "scribe.h"

class IncidenceLogger {
  public:
    IncidenceLogger(std::set<std::string> groups, int num_serotypes, std::string output_folder);
    void log(std::string group, int serotype);
    void flush(int current_day);

  private:
    std::set<std::string> groups_;
    int num_serotypes_;
    std::unordered_map<std::string, std::vector<int>> incidents_;
    Scribe scribe_;
};