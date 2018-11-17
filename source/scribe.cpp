#include <fstream>
#include <string>
#include <vector>
#include <boost/filesystem.hpp>

#include "scribe.h"

namespace fs = boost::filesystem;

Scribe::Scribe(std::string dir) 
  : dir_         {dir}
  , filestreams_ {}
{
  if (fs::exists(dir_)) 
    throw std::runtime_error {"Scribe::Scribe: " + dir + " already exists." };
  fs::create_directories(dir_);
}

std::ofstream& Scribe::get_filestream(std::string filename) {
  if (filestreams_.count(filename) > 0) {
    return filestreams_.at(filename);
  } else {
    filestreams_.emplace(
      std::piecewise_construct, 
      std::make_tuple(filename),
      std::make_tuple() // so we can construct filestream in place with no arguments
    );
    filestreams_[filename].open(dir_ + '/' + filename, std::ios::app);
    return filestreams_[filename];
  }
}

void Scribe::copy_file(std::string filepath) {
  fs::path source {filepath};
  fs::path dest_dir {dir_};
  fs::copy_file(source, dest_dir / source.filename());
}

void Scribe::copy_folder(std::string folderpath) {
  fs::path source = fs::path {folderpath};
  fs::path target = fs::path {dir_} / source.filename();
  fs::create_directory(target);
  for (auto&& file : fs::directory_iterator { source }) {
    fs::copy_file(file.path(), target / file.path().filename());
  }
}
