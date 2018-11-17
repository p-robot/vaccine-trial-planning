#pragma once
#include <set>
#include <string>
#include <map>
#include <unordered_map>
#include <vector>


class Vaccine {
public:
  enum class EffectType {
    Susceptibility, 
    Infectiousness,
    Duration
  };
  
  static std::vector<std::string> get_names();
  static bool targets_serotype(std::string name, EffectType type, std::string serotype_name); 
  static double get_reduction(std::string name, EffectType type, int num_doses); // constant across serotypes, dose-dependent
  static void set_max_reduction(std::string name, EffectType type, double max_reduction);

  static void configure(std::string config_file);

private:
  class Effect {
  public:
    Effect();
    Effect(const std::set<std::string> serotypes, double max_reduction, const std::vector<double>& relative_reductions);
    bool targets_serotype(std::string serotype_name) const;
    double get_reduction(int num_doses) const;
    void set_max_reduction(double max_reduction);
  private:
    void calculate_reductions();
    double max_reduction_;
    std::vector<double> relative_reductions_; // indexed by number of doses received
    std::vector<double> reductions_;          // indexed by number of doses received
    std::set<std::string> targeted_serotypes_;
  };

  Vaccine(const std::string &name);
  void set_effect(EffectType type, const std::set<std::string> serotypes, double max_reduction, 
                  const std::vector<double>& relative_reductions);

  std::string get_name() const;
  bool targets_serotype(EffectType type, std::string serotype_name) const;
  double get_reduction(EffectType type, int num_doses) const;
  void set_max_reduction(EffectType type, double max_reduction);

  std::string name_;
  std::map<EffectType, Effect> effects_;

  static std::unordered_map<std::string, Vaccine> vaccines_;
};