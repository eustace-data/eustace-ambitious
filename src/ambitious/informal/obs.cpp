#include "ambitious/debuglog/debuglog.hpp"
#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(obs_cpp, "$Revision: 1273 $")

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
//  VersionRegistry::log();
  //  int64_t epoch_calyear_obs;
  int64_t start_day;
  int64_t end_day;

  std::string config_filename_default("config.cfg");
  LOG_(std::endl);
  Commandline cmdline(argc, argv, config_filename_default,
                      std::string("obs"), std::string("$Rev: 1273 $"));
  LOG_(std::endl);
  if (cmdline.handle_options()) {
    VERSIONS_PURGE();
    return cmdline.ret();
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


  TimerHierarchy timer;
  timer.tic("Total");

  timer.tic("Setup input manager");
  // Make an input manager with specified base path for input data
  ModelObservationSources obs_sources(cmdline.raw_datafile_path().c_str(),
                                      NULL,
                                      cmdline.epoch_calyear_obs(),
                                      &timer);

  // Start and end day numbers since 01/01/1850
  //  int64_t startday = 57374;
  //  int64_t endday = 57374;
  //  int64_t endday = 57401;
  int64_t startday = start_day;
  int64_t endday = end_day;
  // 53113 missing in some source?
  //  startday = 53110;
  //  endday = 53119;

  obs_sources
      .add(AnalysisInput(
          AnalysisInput::source_satellite_ocean, 
          AnalysisInput::observable_tmean, 
          startday, endday))
      .add(AnalysisInput(
          AnalysisInput::source_satellite_land, 
          AnalysisInput::observable_tmin,
          startday, endday))
      .add(AnalysisInput(
          AnalysisInput::source_satellite_land, 
          AnalysisInput::observable_tmax,
          startday, endday))
      .add(AnalysisInput(
          AnalysisInput::source_satellite_ice, 
          AnalysisInput::observable_tmin,
          startday, endday))
      .add(AnalysisInput(
          AnalysisInput::source_satellite_ice, 
          AnalysisInput::observable_tmax,
          startday, endday))
      .add(AnalysisInput(
          AnalysisInput::source_insitu_land, 
          AnalysisInput::observable_tmin, 
          startday, endday))
      .add(AnalysisInput(
          AnalysisInput::source_insitu_land, 
          AnalysisInput::observable_tmax, 
          startday, endday))
      .add(AnalysisInput(
          AnalysisInput::source_insitu_ocean, 
          AnalysisInput::observable_tmean, 
          startday, endday));
  timer.toc();

  timer.tic("Read multiple days");
  //  for (int day=startday; day <= endday; ++day) {
  for (int day=startday; day <= -1; ++day) {
    std::cout << "********* " << "Day " << to_string(day) << std::endl;
    //    timer.tic("Day " + to_string(day));
    obs_sources.retrieve_day(day);
    obs_sources.summary_log();

  }
  timer.toc();

  timer.toc(-1);
  std::cout << "\nTimings\n";
  timer.print();

  VERSIONS_PURGE();
  return 0;
}


