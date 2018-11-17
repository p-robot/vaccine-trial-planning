#include <cstdlib>
#include <numeric>
#include <string>
#include <boost/property_tree/json_parser.hpp>
#include "configuration.h"

boost::property_tree::ptree load_config(const std::string config_file, const std::string key) {
  boost::property_tree::ptree ptree;
  boost::property_tree::read_json(config_file, ptree);
  return key == "" ? ptree : ptree.get_child(key);
}



