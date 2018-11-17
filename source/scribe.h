#pragma once
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <stdexcept>
#include <string>
#include <vector>

class Scribe {
public:
  Scribe(std::string dir);
  std::ofstream& get_filestream(std::string filename);
  void copy_file(std::string filepath);
  void copy_folder(std::string folderpath);

  template<typename T>
  void note_value(T value, std::string filename) {
    auto& f = get_filestream(filename);
    f << value << '\n';
    f.flush();
  }

  template<typename T>
  void note_values(std::vector<T> values, std::string filename) {
    auto& f = get_filestream(filename);
    write_vector_as_csv(values, f);
    f.flush();
  }

  template<typename T>
  static void write_vector_as_csv(std::vector<T> const &vector, std::ofstream &out) {
    for (auto it = vector.begin(); it != vector.end(); it++) {
      if (it != vector.begin())
        out << ", ";
      out << *it;
    }
    out << '\n';
  }

private:
  std::string dir_;
  std::unordered_map<std::string, std::ofstream> filestreams_;
};