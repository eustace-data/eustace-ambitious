#include "ambitious/debuglog/debuglog.hpp"
#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(read_one_day_test, "$Revision: 1270 $")

#include <iostream>
#include "ambitious/timer/timer.hpp"
#include "ambitious/to_string/to_string.hpp"
#include "ambitious/ostreamer/ostreamer.hpp"
#include <string.h>
#include "ambitious/po/po.hpp"
#include "ambitious/observe/observe.hpp"
#include <eustace/analysis/inputmanager.h>
#include <Eigen/Dense>

using namespace EUSTACE;





int main(int argc, char* argv[])
{
  TimerHierarchy timer;
  timer.tic("Total");

  int64_t start_day;
  int64_t end_day;
  
  std::string config_filename_default("config.cfg");
  LOG_(std::endl);
  Commandline cmdline(argc, argv, config_filename_default.c_str());
  LOG_(std::endl);
  cmdline.set_name_version("obs", "$Rev: 1270 $");
  LOG_(std::endl);
  try{
    if (cmdline.handle_base_options()) {
      VERSIONS_PURGE();
      return 0;
    }
    LOG_(std::endl);
  }
  catch(std::exception& e)
  {
    std::cout << e.what() << std::endl;
    VERSIONS_PURGE();
    return 1;
  }    
  LOG_(std::endl);
  
  start_day = DateTimeHelper(cmdline.epoch_calyear_obs()).get_day(cmdline.start_calyear(), 1, 1);
  LOG_(std::endl);
  end_day = DateTimeHelper(cmdline.epoch_calyear_obs()).get_day(cmdline.end_calyear(), 12, 31);
  LOG_(std::endl);

  //    if (cmdline.vm().count("raw_datafile_path")) {
  std::cout << "Raw datafile path is: " << cmdline.raw_datafile_path() << std::endl;
  std::cout << "Python mesh filename prefix is: " << cmdline.python_mesh_prefix() << std::endl;
  //    }

  // Start and end day numbers since 1850-01-01
  start_day = 57374;
  end_day = 57374;
  //  end_day = 57401;

  timer.tic("Setup input manager");
  // Make an input manager with specified base path for input data
  ModelObservationSources obs_sources(cmdline.raw_datafile_path().c_str(),
                                      NULL,
                                      cmdline.epoch_calyear_obs(),
                                      &timer);

  obs_sources
      .add(AnalysisInput(
          AnalysisInput::source_satellite_ocean, 
          AnalysisInput::observable_tmean, 
          start_day, end_day))
      .add(AnalysisInput(
          AnalysisInput::source_satellite_land, 
          AnalysisInput::observable_tmin,
          start_day, end_day))
      .add(AnalysisInput(
          AnalysisInput::source_satellite_land, 
          AnalysisInput::observable_tmax,
          start_day, end_day))
      .add(AnalysisInput(
          AnalysisInput::source_satellite_ice, 
          AnalysisInput::observable_tmin,
          start_day, end_day))
      .add(AnalysisInput(
          AnalysisInput::source_satellite_ice, 
          AnalysisInput::observable_tmax,
          start_day, end_day))
      .add(AnalysisInput(
          AnalysisInput::source_insitu_land, 
          AnalysisInput::observable_tmin, 
          start_day, end_day))
      .add(AnalysisInput(
          AnalysisInput::source_insitu_land, 
          AnalysisInput::observable_tmax, 
          start_day, end_day))
      .add(AnalysisInput(
          AnalysisInput::source_insitu_ocean, 
          AnalysisInput::observable_tmean, 
          start_day, end_day));
  
  timer.toc();

  //-----------------------------------------------------------------
  // Analysis phase - use the data

  timer.tic("Retrieve data");
  // Retrieve a specific day of data

  for (int day=start_day; day <= end_day; ++day) {
    std::cout << "*********" << std::endl;
    timer.tic("Day " + to_string(day));
    obs_sources.retrieve_day(day);
    obs_sources.summary_log();
    obs_sources.display_day();
    timer.toc();
  }
    
  timer.toc();

  timer.toc(-1);
  std::cout << "\nTimings\n";
  timer.print();

  exit(0);
}


