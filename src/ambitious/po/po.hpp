#ifndef PO_HPP
#define PO_HPP

#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(po_hpp, "$Revision: 1283 $")

#include "ambitious/debuglog/debuglog.hpp"
#include <iostream>
#include "ambitious/ostreamer/ostreamer.hpp"
#include <string.h>
#include "ambitious/po_tool/po_tool.hpp"
#include "ambitious/to_string/to_string.hpp"

namespace po = boost::program_options;

class Commandline : public CommandlineBase {
protected:
  std::string raw_datafile_path_;
  std::string python_mesh_prefix_;
  std::string break_info_path_;
  int64_t epoch_calyear_obs_;
  int64_t start_calyear_;
  int64_t end_calyear_;
  int64_t start_dayoffset_;
  int64_t end_dayoffset_;
  int64_t macro_level_;
  int64_t micro_level_;

  void add_my_options();

public:
  Commandline(int ac, char* av[], const std::string& config_file)
    : CommandlineBase(ac, av, config_file) {
    add_my_options();
  }
  Commandline(int ac, char* av[],
              const std::string& config_file,
              const std::string& program_name,
              const std::string& program_version)
      : CommandlineBase(ac, av, config_file,
                        program_name, program_version) {
    add_my_options();
  }

  const std::string & raw_datafile_path() { return raw_datafile_path_; }
  const std::string & python_mesh_prefix() { return python_mesh_prefix_; }
  const std::string & break_info_path() { return break_info_path_; }
  const int64_t & epoch_calyear_obs() { return epoch_calyear_obs_; }
  const int64_t & start_calyear() { return start_calyear_; }
  const int64_t & end_calyear() { return end_calyear_; }
  const int64_t & start_dayoffset() { return start_dayoffset_; }
  const int64_t & end_dayoffset() { return end_dayoffset_; }
  const int64_t & macro_level() { return macro_level_; }
  const int64_t & micro_level() { return micro_level_; }

};

#endif
