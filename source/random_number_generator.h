#pragma once
#include <boost/random.hpp>

// Provides random draws from commonly-used distributions.
class RandomNumberGenerator {
public:
  // Seed the underlying pseudo random number generator for reproducible results.
  // Or, seed with random noise.
  RandomNumberGenerator(int seed);
  
  // Draw from Bernoulli(p), i.e. returns true with probability p
  bool draw_bernoulli(double p);
  
  // Draw from Uniform distribution over integers in [a, b]
  int draw_uniform_int(int a, int b);
  
  // Draw from Uniform distribution over reals in [a, b)
  double draw_uniform_real(double a, double b);

  // Draw from Categorical (aka Multinoulli), distribution i.e. returns k with probability weights[k]. 
  // Weights will be normalized to sum to 1.
  int draw_categorical(std::vector<double> weights);
  
  // Draw from Multinomial distribution, with probability vector proportional to weights.
  // Weights will be normalized to sum to 1.
  std::vector<int> draw_multinomial(int num_trials, std::vector<double> weights);

  // Draw from Poisson distribution with rate lambda
  int draw_poisson(double lambda);
  
  // Draw from Exponential distribution with specified mean.
  double draw_exponential (double mean);

private:
  boost::random::mt19937 engine_; // A Mersenne twister pseudorandom number generator.
};
