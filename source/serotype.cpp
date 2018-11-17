#include <algorithm>
#include <regex>
#include <set>
#include <string>
#include <stdexcept>
#include <vector>
#include <boost/filesystem.hpp>
#include "configuration.h"
#include "serotype.h"

std::vector<Serotype> Serotype::serotypes_;
int Serotype::num_serotypes_;
double Serotype::base_duration_;
double Serotype::min_duration_;
double Serotype::max_duration_;
double Serotype::min_rank_;
double Serotype::max_rank_;

Serotype::Serotype(const std::string &name, const int serogroup, const double &rank) 
  : name_ {name}
  , serogroup_ {serogroup}
  , rank_ {rank}
{
  update_duration();
}
std::string Serotype::get_name() const {
  return name_;
}
int Serotype::get_serogroup() const {
  return serogroup_;
}
double Serotype::get_rank() const {
  return rank_;
}
double Serotype::get_duration() const {
  return duration_;
}
void Serotype::set_rank(double rank) {
  if (rank < get_min_rank() || rank > get_max_rank()) {
    throw std::runtime_error { "Serotype::set_rank: Rank " + std::to_string(rank) + "outside of allowed range"};
  } 
  rank_ = rank;
  update_duration();
}
void Serotype::update_duration() {
  duration_ = min_duration_ 
    + (max_duration_ - min_duration_) * (max_rank_ - rank_) / (max_rank_ - min_rank_);
}

// Static functions
int Serotype::get_num_serotypes() {
  return num_serotypes_;
}
std::string Serotype::get_name(int index) {
  return serotypes_[index].get_name();
}
int Serotype::get_serogroup(int index) {
  return serotypes_[index].get_serogroup();
}
double Serotype::get_rank(int index) {
  return serotypes_[index].get_rank();
}
double Serotype::get_duration(int index) {
  return serotypes_[index].get_duration();
}
void Serotype::set_rank(int index, double rank) {
  serotypes_[index].set_rank(rank);
}
void Serotype::set_ranks(std::vector<double> ranks) {
  if (ranks.size() != num_serotypes_) {
    throw std::runtime_error {"Serotype::set_ranks: Wrong number of ranks given."};
  }
  for (int i = 0; i < num_serotypes_; i++) {
    serotypes_[i].set_rank(ranks[i]);
  }
}

double Serotype::get_base_duration() {
  return base_duration_;
}
double Serotype::get_min_rank() {
  return min_rank_;
}
double Serotype::get_max_rank() {
  return max_rank_;
}

int Serotype::extract_serogroup(std::string serotype) {
  std::regex regex {"^([0-9]+)[A-Z]*$"};
  std::smatch match;
  bool matched = std::regex_match(serotype, match, regex);
  if (!matched)
    throw std::runtime_error {"Serotype " + serotype + " does not follow pattern of number followed by an optional capital letters, e.g. '12', '12A', '12AB'." };
  return std::stoi(match[1].str()); // return first group matched (the number portion)
}

void Serotype::configure(const std::string config_file) {
  auto config_folder = boost::filesystem::path(config_file).parent_path();
  auto ptree = load_config(config_file, "serotype");
  
  // get and check duration parameters
  base_duration_ = ptree.get<double>("base_duration");
  min_duration_  = ptree.get<double>("min_duration");
  max_duration_  = ptree.get<double>("max_duration");
  if (min_duration_ < base_duration_) {
    throw std::runtime_error {"Minimum intrinsic duration cannot be less than the base duration."};
  }
  if (max_duration_ < min_duration_) {
    throw std::runtime_error {"Maximum intrinsic duration cannot be less than minimum intrinsic duration."};
  }

  // read in the list of serotypes
  auto types_path = (config_folder / ptree.get<std::string>("serotypes_file")).string();
  auto types = load_vector<std::string>(types_path, "serotypes");

  auto ranks_path = (config_folder / ptree.get<std::string>("ranks_file")).string();
  auto ranks = load_vector<double>(ranks_path, "ranks");

  if (ranks.size() != types.size()) {
    throw std::runtime_error {"Number of fitness ranks does not match number of serotypes."};
  }

  // initialize static variables first
  num_serotypes_ = (int) types.size();
  min_rank_ = 1.0;
  max_rank_ = num_serotypes_;

  // then instantiate serotypes
  for (int i = 0; i < types.size(); i++) {
    Serotype s { types[i], extract_serogroup(types[i]), ranks[i] };
    serotypes_.push_back(s);
  }

  // check that their ranks are valid
  for (auto &s : serotypes_) {
    if (s.get_rank() < min_rank_ || s.get_rank() > max_rank_) {
      throw std::runtime_error {"Serotype::configure: Rank for serotype " + s.get_name() + " outside of allowed range."};
    }
  }
}
