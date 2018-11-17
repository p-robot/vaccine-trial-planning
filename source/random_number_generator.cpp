#include <stdexcept>
#include <string>
#include <vector>
#include "boost/random.hpp"
#include "random_number_generator.h"

RandomNumberGenerator::RandomNumberGenerator(int seed) 
  : engine_(seed) {}

bool RandomNumberGenerator::draw_bernoulli(double p) {
  if (p < 0 || p > 1)
    throw std::runtime_error {"RandomNumberGenerator::draw_bernoulli: Bernoulli probability must be between 0 and 1. (Got " + std::to_string(p) + ")" };
  boost::random::bernoulli_distribution<> bernoulli (p);
  return bernoulli(engine_);
}

int RandomNumberGenerator::draw_uniform_int(int a, int b) {
  boost::random::uniform_int_distribution<int> uniform_int (a, b);
  return uniform_int(engine_);
}

double RandomNumberGenerator::draw_uniform_real(double a, double b) {
  boost::random::uniform_real_distribution<double> uniform_real (a, b);
  return uniform_real(engine_);
}

int RandomNumberGenerator::draw_categorical(std::vector<double> weights) {
  boost::random::discrete_distribution<int> discrete (weights.begin(), weights.end());
  return discrete(engine_);
}

std::vector<int> RandomNumberGenerator::draw_multinomial(int num_trials, 
                                                         std::vector<double> weights) {
  int num_categories = weights.size();
  std::vector<int> counts (num_categories, 0);
  for (int i = 0; i < num_trials; i++) {
    int c = draw_categorical(weights);
    counts[c]++;
  }
  return counts;
}

int RandomNumberGenerator::draw_poisson(double lambda) {
  boost::random::poisson_distribution<int> poisson (lambda);
  return poisson(engine_);
}

double RandomNumberGenerator::draw_exponential(double mean) {
  boost::random::exponential_distribution<double> exponential (1.0 / mean);
  return exponential(engine_);
}