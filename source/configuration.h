#pragma once
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/property_tree/ptree.hpp>

// Expects JSON file
boost::property_tree::ptree load_config(const std::string config_file, const std::string key="");

// For a given boost property tree node, converts its list of child nodes to a vector
template<typename T> 
std::vector<T> as_vector(const boost::property_tree::ptree &node) {
  std::vector<T> v;
  for (auto& child : node)
    v.push_back(child.second.get_value<T>()); // child.first is name of node, 
                                              // child.second is a ptree
  return v;
}

template<typename T> 
std::set<T> as_set(const boost::property_tree::ptree &node) {
  std::set<T> s;
  for (auto& child : node)
    s.insert(child.second.get_value<T>()); // child.first is name of node, 
                                              // child.second is a ptree
  return s;
}

template<typename T>
std::vector<T> load_vector(const std::string config_file, const std::string key, int expected_size=-1) {
  std::vector<T> vector = as_vector<T>(load_config(config_file, key));
  if (expected_size > 0 && vector.size() != expected_size) {
    throw std::runtime_error { "Configuration::load_vector: " + std::to_string(vector.size()) + " elements specified in " + config_file 
      + ". Expected " + std::to_string(expected_size) + " counts." };
  } 
  return vector;
}