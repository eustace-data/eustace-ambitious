#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(po_cpp, "$Revision: 1283 $")

#include "ambitious/debuglog/debuglog.hpp"
#include <iostream>
#include "ambitious/ostreamer/ostreamer.hpp"
#include <string.h>
#include "ambitious/po_tool/po_tool.hpp"
#include "ambitious/to_string/to_string.hpp"

#include "ambitious/po/po.hpp"

namespace po = boost::program_options;

void Commandline::add_my_options() {  
  // Declare a group of options that will be allowed only on command
  // line
  add_generic();
  
  // Declare a group of options that will be allowed both on command
  // line and in config file
  add_config()
      ("raw-datafile-path,D", 
       po::value< std::string >(&raw_datafile_path_)->default_value("data/rawbinary8"), 
       "raw datafile path")
      ("python-mesh-prefix", 
       po::value< std::string >(&python_mesh_prefix_)->default_value("data/eustace_mesh_level"), 
       "python mesh filename prefix")
      ("break-info-path", 
       po::value< std::string >(&break_info_path_)->default_value(
           "data/eustace_stations_global_R001127_daily_status.nc"), 
       "land station status filename path")
      ("epoch-calyear-obs", 
       po::value< int64_t >(&epoch_calyear_obs_)->default_value(1850), 
       "epoch calendar year for the observations")
      ("start-calyear", 
       po::value< int64_t >(&start_calyear_)->default_value(1850), 
       "start calendar year")
      ("end-calyear", 
       po::value< int64_t >(&end_calyear_)->default_value(2021), 
       "end calendar year")
      ("start-dayoffset", 
       po::value< int64_t >(&start_dayoffset_)->default_value(0), 
       "start offset, in days")
      ("end-dayoffset", 
       po::value< int64_t >(&end_dayoffset_)->default_value(-1), 
       "end offset, in days, may be negative")
      ("macro-level", 
       po::value< int64_t >(&macro_level_)->default_value(0), 
       "macro-level")
      ("micro-level", 
       po::value< int64_t >(&micro_level_)->default_value(4), 
       "micro-level")
      ;
  
  // Hidden options, will be allowed both on command line and
  // in config file, but will not be shown to the user.
  add_hidden();
}
