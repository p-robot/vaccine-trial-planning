#pragma once
#include <ctime>
#include <string>
#include <vector>
#include "simulation.h"

class PrevalenceMatcher {
public:  
  enum class Mode {
    FitTransmission,
    FitVaccine
  };

  PrevalenceMatcher(Mode mode, int seed, std::string config_file, std::string output_dir);
  void match();

private:
  void configure(std::string config_file);

  std::vector<double> calculate_prevalences(std::vector<unsigned int> counts);
  std::vector<double> calculate_prevalences(std::vector<Simulation::Snapshot> snapshots, unsigned int max_samples);
  std::vector<double> calculate_errors(std::vector<double> estimated, std::vector<double> target);

  void update_transmission_parameters();
  void update_vaccine_parameters();

  void adjust_weight(double &weight, const double target, const double error, const double previous_error);
  void adjust_weights(std::vector<double>& weights, const std::vector<double>& targets, 
                      const std::vector<double>& errors, const std::vector<double>& previous_errors);
  void adjust_beta(double &beta, const double weight, const double error);
  void adjust_ranks(std::vector<double> &ranks, const std::vector<double>& weights, const std::vector<double>& errors);
  void adjust_max_susc_reduction(double &max_susceptibility_reduction, const double weight, const double error);
  void save_snapshot();

  Mode mode_;
  int seed_;
  std::string config_file_;
  std::string output_folder_;

  unsigned int max_iterations_;
  unsigned int max_samples_;
  
  double warm_up_;
  double cool_down_;
  double temperature_threshold_;
  
  std::vector<unsigned int> observed_counts_;
  std::vector<double> target_prevalences_;
  std::time_t start_time_;

  // updated every iteration
  unsigned int iteration_;
  double beta_;
  double beta_weight_;

  std::vector<double> ranks_;
  std::vector<double> rank_weights_;

  std::string vaccine_to_fit_;
  double max_susceptibility_reduction_;
  double max_susceptibility_reduction_weight_;

  std::vector<double> estimated_prevalences_;
  std::vector<double> prevalence_errors_;
  std::vector<double> previous_prevalence_errors_;
  double loglikelihood_;

  Scribe scribe_;
};