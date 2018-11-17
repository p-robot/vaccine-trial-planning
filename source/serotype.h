#pragma once
#include <set>
#include <string>
#include <vector>

class Serotype {
public:  
  static int get_num_serotypes();
  
  static std::string get_name(int index);
  static int get_serogroup(int index);
  static double get_rank(int index);
  static double get_duration(int index);
  static void set_rank(int index, double rank);
  static void set_ranks(std::vector<double> ranks);
  
  static double get_base_duration();
  static double get_min_rank();
  static double get_max_rank();
  
  static void configure(std::string config_file);

private: 
  Serotype(const std::string &name, const int serogroup, const double &rank); 
//  Serotype(const Serotype&) = default;            
//  Serotype& operator=(const Serotype&) = default; 

  std::string get_name() const;
  int get_serogroup() const;
  double get_rank() const;
  double get_duration() const;
  void set_rank(double rank);
  void update_duration();

  std::string name_;
  int serogroup_;
  double rank_;
  double duration_;

  static int extract_serogroup(std::string serotype);

  static std::vector<Serotype> serotypes_;
  static int num_serotypes_;
  static double base_duration_;
  static double min_duration_;
  static double max_duration_;
  static double min_rank_;
  static double max_rank_;
  static int max_serogroup_;
};