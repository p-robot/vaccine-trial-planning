#include <fstream>
#include <random>
#include <set>
#include <string>
#include <vector>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include "configuration.h"
#include "serotype.h"
#include "vaccine.h"
#include "host.h"
#include "simulation.h"
#include "prevalence_matcher.h"

bool get_command_line_options(int argc, char** argv, 
  std::string &config_file, std::string &output_folder, std::string &task, std::vector<int> &seeds, int& num_trials) {

  // set up command line options  
  namespace po = boost::program_options;
  po::options_description desc {"Allowed options"};
  desc.add_options()
    ("help,h",                                                "show help message")
    ("config-file,c",   po::value<std::string>()->required(), "specify path to configuration file")
    ("output-folder,o", po::value<std::string>()->required(), "specify path to output folder")
    ("task,t",          po::value<std::string>()->required(), "specify task to be either 'simulate' or 'fit'")
    ("seed,s",          po::value<int>(),                     "set the seed for the random number generator")
    ("num-trials,n",    po::value<int>()->default_value(1),   "set the number of trials to run")
  ;

  po::variables_map vm;
  try {
    // read into a variables_map
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
      std::cout << desc << std::endl;
      return false;
    } else {
      po::notify(vm); // do this *after* we check for the help option
    }

    // store parameters in our variables
    config_file = vm["config-file"].as<std::string>();
    output_folder = vm["output-folder"].as<std::string>();
    task = vm["task"].as<std::string>();

    num_trials = vm["num-trials"].as<int>();
    if (num_trials <= 0) {
      std::cerr << "Specify a positive number for --num-trials option." << std::endl;
      return false;
    }

    std::random_device rd;
    for (int i = 0; i < num_trials; i++) {
      int s = vm.count("seed") ? vm["seed"].as<int>() + i : rd();
      seeds.push_back(s);
    }
    return true;
  } catch (std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return false;
  }
}

std::set<std::string> get_referenced_files(boost::property_tree::ptree pt) {
  std::set<std::string> referenced;
  for (auto it = pt.begin(); it != pt.end(); it++) {
    auto key = it->first;
    if (boost::algorithm::ends_with(key, "file")) {
      referenced.insert(it->second.get_value<std::string>());
    } else {
      auto more = get_referenced_files(it->second);
      referenced.insert(more.begin(), more.end());
    }
  }
  return referenced;
}

int main(int argc, char** argv) {
  namespace fs = boost::filesystem;

  // get command line options
  std::string config_file, output_folder, task;
  int num_trials = 0;
  std::vector<int> seeds;
  bool got_options = get_command_line_options(argc, argv, config_file, output_folder, task, seeds, num_trials);
  if (!got_options) 
    return 1;

  // create the output directory for this run
  if (fs::exists(output_folder)) 
    throw std::runtime_error { "Output directory already exists." };
  fs::create_directory(output_folder);

  // save options for future reference
  std::ofstream command_file { (fs::path(output_folder) / "command.txt").string() };
  for (int i = 0; i < argc; i++) {
    command_file << argv[i] << std::endl;
  }
  command_file.close();

  // make copies of relevant configuration files
  auto working_config_folder = fs::path(output_folder) / "configuration";
  fs::create_directory(working_config_folder);

  fs::copy_file(config_file, working_config_folder / "configuration.json");  

  auto referenced_files = get_referenced_files(load_config(config_file));
  auto starting_config_folder = fs::path(config_file).parent_path();
  for (auto &filename : referenced_files) {
    fs::copy_file(starting_config_folder / filename, working_config_folder / filename);
  }

  // configure
  config_file = (working_config_folder / "configuration.json").string(); // reset config_file be our working copy
  Serotype::configure(config_file);
  Vaccine::configure(config_file);
  Host::configure(config_file);

  // perform task
  for(int i = 0; i < num_trials; i++) {
    std::string output_folder_i;
    output_folder_i = (fs::path(output_folder) / ("trial-" + std::to_string(i))).string(); 
     
    if (task == "simulate") {
      Simulation simulation {seeds[i], config_file, output_folder_i};
      simulation.run();
    } else if (task == "fit-transmission") {
      PrevalenceMatcher prevalence_matcher { PrevalenceMatcher::Mode::FitTransmission, seeds[i], config_file, output_folder_i};
      prevalence_matcher.match();
    } else if (task == "fit-vaccine") {
      PrevalenceMatcher prevalence_matcher { PrevalenceMatcher::Mode::FitVaccine, seeds[i], config_file, output_folder_i};
      prevalence_matcher.match();
    } else {
      std::cerr << "Unrecognized task " << "'" << task << "'" << std::endl;
      return 0;
    }
  }
  return 0;
}
