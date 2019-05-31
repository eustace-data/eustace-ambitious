
#include <eustace/analysis/inputmanager.h>
#include <iostream>
#include <stdlib.h>

using namespace EUSTACE;

int main(int argc, char* argv[])
{
  // Get input directory from commandline
  if (argc != 2)
  {
    std::cerr << "Usage: example_analysis_input datadirectory" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Make an input manager with specified base path for input data
  AnalysisInputManager inputmanager(argv[1]);

  //-----------------------------------------------------------------
  // Specification phase - say which data we're interested in

  // Start and end day numbers since 01/01/1850
  int64_t startday = 57374;
  int64_t endday = 57401;

  // In this example look at satellite-derived data for ocean
  inputmanager.Specify(AnalysisInput(
    AnalysisInput::source_satellite_ocean, 
    AnalysisInput::observable_tmean, 
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

  // And in-situ ocean data over the same time period (daily mean)
  inputmanager.Specify(AnalysisInput(
    AnalysisInput::source_insitu_ocean, 
    AnalysisInput::observable_tmean, 
    startday, endday));

  //-----------------------------------------------------------------
  // Analysis phase - use the data

  // Print summary meta data
  std::cout << std::endl;

  std::cout << "Summary" << std::endl;
  std::cout << "-------" << std::endl;
  std::cout << "Total locations satellite-derived ocean: " << inputmanager.TotalLocations(AnalysisInput::source_satellite_ocean) << std::endl;
  std::cout << "Total locations in-situ land stations: " << inputmanager.TotalLocations(AnalysisInput::source_insitu_land) << std::endl;
  std::cout << "-------" << std::endl;
  std::cout << "Total obs satellite-derived ocean (Tmean): " << inputmanager.TotalObservations(AnalysisInput::source_satellite_ocean, AnalysisInput::observable_tmean) << std::endl;
  std::cout << "Total obs in-situ land (Tmin): " << inputmanager.TotalObservations(AnalysisInput::source_insitu_land, AnalysisInput::observable_tmin) << std::endl;
  std::cout << "Total obs in-situ land (Tmax): " << inputmanager.TotalObservations(AnalysisInput::source_insitu_land, AnalysisInput::observable_tmax) << std::endl;
  std::cout << "Total obs in-situ ocean (Tmean): " << inputmanager.TotalObservations(AnalysisInput::source_insitu_ocean, AnalysisInput::observable_tmean) << std::endl;
  std::cout << "-------" << std::endl;
  std::cout << "Total locations in-situ ocean (Tmean) at day 57377: " << inputmanager.TotalDailyMobileLocations(AnalysisInput::source_insitu_ocean, AnalysisInput::observable_tmean, 57377) << std::endl;
  std::cout << "-------" << std::endl;

  // Retrieve a specific day of land data
  AnalysisInputRetrieval insitu_land_retrieval;
  inputmanager.RetrieveDay(insitu_land_retrieval, 
			   AnalysisInput::source_insitu_land,
			   AnalysisInput::observable_tmin,
			   57377);
  
  // Announce what we're doing to console and say total observations read
  std::cout << std::endl;
  std::cout << "In-situ land data" << std::endl;
  std::cout << "Loaded observations: " << insitu_land_retrieval.NumObservations() << std::endl;

  // Have a look at some numbers
  for(AnalysisInputIterator insitu_land_iterator(insitu_land_retrieval); insitu_land_iterator.HasElement(); insitu_land_iterator.Next())
  {
    // Location ID
    uint64_t location_id = insitu_land_iterator.Identifier();
    
    // Demonstrate automatic location lookup to get latitude and longitude
    // This is the same as insitu_land_retrieval.ComputeLocationLookup(location_id)
    const LocationEntry& location = insitu_land_iterator.Location();

    // Measurement value
    double measurement = insitu_land_iterator.Measurement();

    // Uncorrelated uncertainty
    double uncertainty_uncorrelated = insitu_land_iterator.UncertaintyUncorrelated();

    // Print some info to console
    std::cout 
      << "lat: " << location.latitude
      << " lon: " << location.longitude
      << " temperature: " << measurement
      << " uncertainty (uncorrelated): " << uncertainty_uncorrelated
      << std::endl;
  }


  // Retrieve a specific day of satellite-derived ocean data
  AnalysisInputRetrieval satellite_ocean_retrieval;
  inputmanager.RetrieveDay(satellite_ocean_retrieval, 
			   AnalysisInput::source_satellite_ocean,
			   AnalysisInput::observable_tmean,
			   57377);
  
  // Announce what we're doing to console
  std::cout << std::endl;
  std::cout << "Satellite-derived ocean data" << std::endl;
  std::cout << "Loaded observations: " << satellite_ocean_retrieval.NumObservations() << std::endl;

  // Have a look at some numbers
  for(AnalysisInputIterator satellite_ocean_iterator(satellite_ocean_retrieval); satellite_ocean_iterator.HasElement(); satellite_ocean_iterator.Next())
  {
    // Same as above for land but also get a locally correlated uncertainty component
    const LocationEntry& location = satellite_ocean_iterator.Location();
    double measurement = satellite_ocean_iterator.Measurement();
    double uncertainty_uncorrelated = satellite_ocean_iterator.UncertaintyUncorrelated();
    double uncertainty_locally_correlated = satellite_ocean_iterator.UncertaintyLocallyCorrelated(0);

    // Print some info to console
    std::cout 
      << "lat: " << location.latitude
      << " lon: " << location.longitude
      << " temperature: " << measurement
      << " uncertainty (uncorrelated): " << uncertainty_uncorrelated
      << " uncertainty (locally correlated): " << uncertainty_locally_correlated
      << std::endl;
  }

  // Retrieve a specific day of in-situ ocean data
  AnalysisInputRetrieval insitu_ocean_retrieval;
  inputmanager.RetrieveDay(insitu_ocean_retrieval, 
			   AnalysisInput::source_insitu_ocean,
			   AnalysisInput::observable_tmean,
			   57377);
  
  // Announce what we're doing to console
  std::cout << std::endl;
  std::cout << "In-situ ocean data" << std::endl;
  std::cout << "Loaded observations: " << insitu_ocean_retrieval.NumObservations() << std::endl;

  // Have a look at some numbers
  for(AnalysisInputIterator insitu_ocean_iterator(insitu_ocean_retrieval); insitu_ocean_iterator.HasElement(); insitu_ocean_iterator.Next())
  {
    // Same as above for land but also get a locally correlated uncertainty component
    const LocationEntry& location = insitu_ocean_iterator.Location();
    double measurement = insitu_ocean_iterator.Measurement();
    double uncertainty_uncorrelated = insitu_ocean_iterator.UncertaintyUncorrelated();

    // Print some info to console
    std::cout 
      << "lat: " << location.latitude
      << " lon: " << location.longitude
      << " temperature: " << measurement
      << " uncertainty (uncorrelated): " << uncertainty_uncorrelated
      << std::endl;
  }

  // Choose three land station IDs for which there is data (could be any three)
  std::vector<uint64_t> location_id_subset;
  AnalysisInputIterator insitu_land_iterator(insitu_land_retrieval);
  insitu_land_iterator.Next();
  insitu_land_iterator.Next();
  insitu_land_iterator.Next();
  location_id_subset.push_back( insitu_land_iterator.Identifier() );
  insitu_land_iterator.Next();
  location_id_subset.push_back( insitu_land_iterator.Identifier() );
  insitu_land_iterator.Next();
  insitu_land_iterator.Next();
  location_id_subset.push_back( insitu_land_iterator.Identifier() );

  // Retrieve all days of data for these
  AnalysisInputRetrievalCollection collection;
  inputmanager.RetrieveLocations(collection, AnalysisInput::source_insitu_land,  AnalysisInput::observable_tmin, location_id_subset);
  
  // Show this info
  for (uint64_t index = 0; index < 3; index++)
  {
    std::cout << "time sequence for in-situ land station: " << location_id_subset[index] << std::endl;
    const AnalysisInputRetrieval& station_time_sequence = *( collection.Retrievals()[index] );
    for(AnalysisInputIterator time_iterator(station_time_sequence); time_iterator.HasElement(); time_iterator.Next())
    {
      std::cout << "  day: " << time_iterator.Identifier() << "  temperature: " << time_iterator.Measurement() << std::endl;
    }
  }
}
