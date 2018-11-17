#include <iostream>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#include "configuration.h"
#include "vaccine.h"

std::unordered_map<std::string, Vaccine> Vaccine::vaccines_;

/*************************
 * Vaccine::Effect class *
 *************************/
Vaccine::Effect::Effect()
  : targeted_serotypes_  {  }
  , max_reduction_       { 0.0 }
  , relative_reductions_ { 0.0 }
{
  calculate_reductions();
}

Vaccine::Effect::Effect(const std::set<std::string> serotypes, double max_reduction, 
                        const std::vector<double>& relative_reductions)
  : targeted_serotypes_  { serotypes }
  , max_reduction_       { max_reduction }
  , relative_reductions_ { relative_reductions }
{
  calculate_reductions();
}

bool Vaccine::Effect::targets_serotype(std::string serotype_name) const {
  return (targeted_serotypes_.count(serotype_name) > 0) || (targeted_serotypes_.count("all") > 0);
}

double Vaccine::Effect::get_reduction(int num_doses) const {
  if (num_doses < 0)
    throw std::runtime_error { "Vaccine::Effect::get_reduction: num_doses cannot be negative." };
  if (reductions_.size() == 0)
    throw std::runtime_error { "Vaccine::Effect::get_reduction: no reductions specified." };
  return (num_doses < reductions_.size()) ? reductions_[num_doses]: reductions_.back();
}

void Vaccine::Effect::set_max_reduction(double max_reduction) {
  max_reduction_ = max_reduction;
  calculate_reductions();
}

void Vaccine::Effect::calculate_reductions() {
  reductions_.clear();
  reductions_.push_back(0); // no effect at 0 doses
  for (auto r : relative_reductions_) {
    reductions_.push_back(r * max_reduction_);
  }
  reductions_.push_back(max_reduction_);
}

/*****************
 * Vaccine class *
 *****************/
Vaccine::Vaccine(const std::string &name)
  : name_ { name }
{
  // initialize with effects that do nothing
  for (auto& type : { EffectType::Susceptibility, EffectType::Infectiousness, EffectType::Duration }) {
    effects_.insert(std::make_pair(type, Effect()));
  }
}

// instance methods
void Vaccine::set_effect(Vaccine::EffectType type, const std::set<std::string> serotypes, 
                        double max_reduction, const std::vector<double>& relative_reductions) {
  effects_.find(type)->second = Effect(serotypes, max_reduction, relative_reductions);
}
std::string Vaccine::get_name() const {
  return name_;
}
bool Vaccine::targets_serotype(EffectType type, std::string serotype_name) const {
  return effects_.find(type)->second.targets_serotype(serotype_name);
}
double Vaccine::get_reduction(EffectType type, int num_doses) const {
  return effects_.find(type)->second.get_reduction(num_doses);
}
void Vaccine::set_max_reduction(EffectType type, double max_reduction) {
  effects_.find(type)->second.set_max_reduction(max_reduction);
}


// static methods
std::vector<std::string> Vaccine::get_names() {
  std::vector<std::string> keys;
  for (auto &kv: vaccines_) {
    keys.push_back(kv.first);
  }
  return keys;
}

bool Vaccine::targets_serotype(std::string name, EffectType type, std::string serotype_name) {
  return vaccines_.find(name)->second.targets_serotype(type, serotype_name);
}
double Vaccine::get_reduction(std::string name, EffectType type, int num_doses) {
  return vaccines_.find(name)->second.get_reduction(type, num_doses);
}
void Vaccine::set_max_reduction(std::string name, EffectType type, double max_reduction) {
  vaccines_.find(name)->second.set_max_reduction(type, max_reduction);
}

void Vaccine::configure(std::string config_file) {
  auto vaccines_array = load_config(config_file, "vaccines");
  for (auto &v : vaccines_array) {
    // get name of vaccine and create a Vaccine
    auto name = v.second.get<std::string>("name");
    Vaccine vaccine { name };

    // add effects
    std::unordered_map<std::string, EffectType> key_to_type;
    key_to_type["susceptibility_reduction"] = EffectType::Susceptibility;
    key_to_type["infectiousness_reduction"] = EffectType::Infectiousness;
    key_to_type["duration_reduction"]       = EffectType::Duration;

    bool has_no_effect = true;
    for (auto& kv : key_to_type) {
      auto key = kv.first;
      auto type = kv.second;
      if (v.second.count(key) > 0) {
        auto& node = v.second.get_child(key);
        vaccine.set_effect(
          type, 
          as_set<std::string>(node.get_child("serotypes")),
          node.get<double>("max"), 
          as_vector<double>(node.get_child("relative"))
        );
        has_no_effect = false;
      }
    }
    if (has_no_effect) {
      std::cerr << "Warning: At least one vaccine has no effect.";
      std::cerr << "Specify 'susceptibility_reduction', 'infectiousness_reduction', or 'duration_reduction'.";
      std::cerr << std::endl;
    }

    // add vaccine to our list of vaccines
     vaccines_.insert(std::make_pair(name, vaccine));
  }
}