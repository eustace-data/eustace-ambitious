#include "ambitious/debuglog/debuglog.hpp"

#include <eustace/analysis/inputmanager.h>
#include <iostream>
#include "ambitious/timer/timer.hpp"

#include "ambitious/ostreamer/ostreamer.hpp"

#include <string.h>
#include "ambitious/po_tool/po_tool.hpp"
#include "ambitious/po/po.hpp"

using namespace EUSTACE;





//namespace po = boost::program_options;

/*
class Commandline : public CommandlineBase {
protected:
  std::string raw_datafile_path_;

public:
  Commandline(int ac, char* av[], const char * config_file)
    : CommandlineBase(ac, av, config_file) {

    // Declare a group of options that will be allowed only on command
    // line
    add_generic();
    
    // Declare a group of options that will be allowed both on command
    // line and in config file
    add_config()
      ("raw-datafile-path,D", 
       po::value< std::string >(&raw_datafile_path_)->default_value("/work/scratch/eustace/rawbinary3"), 
       "raw datafile path")
      ;
    
    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    add_hidden();
  }

  const std::string & raw_datafile_path() { return raw_datafile_path_; }

};
*/












// Display metadata for a specific day of data
void display_day(AnalysisInputRetrieval& retrieval, std::string name, int64_t day)
{
  PLOG_("Metadata for day " << day << ", " << name
        << ":\t#Obs =\t"
        << retrieval.NumObservations()
        << ",\t#CorrErrors = \t"
        << retrieval.NumLocalCorrelationRanges()
        << std::endl);
}





int main(int argc, char* argv[])
{
  Commandline cmdline(argc, argv, "config.cfg");
  cmdline.set_name_version("meta_data_in_input", "$Revision: 1285 $");
  try{
    if (cmdline.handle_base_options()) {
      return 0;
    }
    
    //    if (cmdline.vm().count("raw_datafile_path")) {
    std::cout << "Raw datafile path is: " << cmdline.raw_datafile_path() << std::endl;
    //    }
  }
  catch(std::exception& e)
  {
    std::cout << e.what() << std::endl;
    return 1;
  }    
    


  TimerHierarchy timer;
  timer.tic("Total");

  timer.tic("Setup input manager");
  // Make an input manager with specified base path for input data
  AnalysisInputManager inputmanager(cmdline.raw_datafile_path().c_str());

  //-----------------------------------------------------------------
  // Specification phase - say which data we're interested in

  // Start and end day numbers since 01/01/1850
  int64_t startday = 57374;
  int64_t endday = startday + 4;
  //  int64_t endday = 57401;

  // In this example look at satellite-derived data for ocean
  inputmanager.Specify(AnalysisInput(
    AnalysisInput::source_satellite_ocean, 
    AnalysisInput::observable_tmean, 
    startday, endday));

  // In this example look at satellite-derived data for land
  inputmanager.Specify(AnalysisInput(
    AnalysisInput::source_satellite_land, 
    AnalysisInput::observable_tmax,
    startday, endday));

  // In this example look at satellite-derived data for ice
  inputmanager.Specify(AnalysisInput(
    AnalysisInput::source_satellite_ice, 
    AnalysisInput::observable_tmax,
    startday, endday));

  // And in-situ land data over the same time period (daily minimum)
  inputmanager.Specify(AnalysisInput(
    AnalysisInput::source_insitu_land, 
    AnalysisInput::observable_tmin, 
    startday, endday));

  // And in-situ land data over the same time period (daily maximum)
  inputmanager.Specify(AnalysisInput(
    AnalysisInput::source_insitu_land, 
    AnalysisInput::observable_tmax, 
    startday, endday));
  timer.toc();

  // And in-situ ocean data over the same time period (daily mean)
  inputmanager.Specify(AnalysisInput(
    AnalysisInput::source_insitu_ocean, 
    AnalysisInput::observable_tmean, 
    startday, endday));

  timer.toc();

  
  std::cout << std::endl;
  timer.tic("Retrieve information");
  // Retrieve a specific day of data
  AnalysisInputRetrieval retrieval;

  for (int day=startday; day <= endday; ++day) {
    std::cout << "*********" << std::endl;
    timer.tic("Day " + std::to_string((long long)day));
    timer.tic("Read ocean data");
    inputmanager.RetrieveDay(retrieval, 
			     AnalysisInput::source_satellite_ocean,
			     AnalysisInput::observable_tmean,
			     day);
    timer.toc().tic("Display metadata");
    
    display_day(retrieval, "ocean", day);
    timer.toc();

    timer.tic("Read land data");
    inputmanager.RetrieveDay(retrieval, 
			     AnalysisInput::source_satellite_land,
			     AnalysisInput::observable_tmax,
			     day);
    timer.toc().tic("Display metadata");
    
    display_day(retrieval, "land", day);
    timer.toc();



    timer.tic("Read ice data");
    inputmanager.RetrieveDay(retrieval, 
			     AnalysisInput::source_satellite_ice,
			     AnalysisInput::observable_tmax,
			     day);
    timer.toc().tic("Display metadata");
    
    display_day(retrieval, "ice", day);
    timer.toc();



    timer.tic("Read ship data");
    inputmanager.RetrieveDay(retrieval, 
			     AnalysisInput::source_insitu_ocean,
			     AnalysisInput::observable_tmean,
			     day);
    timer.toc().tic("Display metadata");
    
    display_day(retrieval, "ship", day);
    timer.toc();

    timer.tic("Read ship data");
    inputmanager.RetrieveDay(retrieval, 
			     AnalysisInput::source_insitu_land,
			     AnalysisInput::observable_tmin,
			     day);
    timer.toc().tic("Display metadata");
    display_day(retrieval, "land", day);
    timer.toc();

    timer.toc();
  }
    
  timer.toc();

  timer.toc(-1);
  std::cout << "\nTimings\n";
  timer.print();

  exit(0);
}


