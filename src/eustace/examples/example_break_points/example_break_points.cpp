#include <eustace/analysis/break_points/insitu_land_break_points.h>
#include <iostream>
#include <stdlib.h>

using namespace EUSTACE;

int main(int argc, char* argv[])
{
  // Get input directory from commandline
  if (argc != 2)
  {
    std::cerr << "Usage: example_break_points FILENAME" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Make the reader 
  BreakPointsReader reader(argv[1]);

  //-----------------------------------------------------------------
  // Return the observation lookup: lat-lon coordinates identifying the insitu source position
  // 0 - lat_0,lon_0 
  // 1 - lat_1, lon_1
  //    .....
  // ntot - lat_ntot, lon_ntot

  std::vector<std::vector<float> > location_lookup = reader.ObservationLocationLookup();
  size_t size = location_lookup[0].size();
  std::cout << "Location look up" << std::endl;
  std::cout << "-------" << std::endl;
  std::cout << "0 - lat_0,lon_0 : (" << location_lookup[0][0] << "," <<  location_lookup[1][0] << ")" <<  std::endl;
  std::cout << "1 - lat_1,lon_1 : (" <<  location_lookup[0][1] << "," <<  location_lookup[1][1] << ")" << std::endl;
  std::cout << "-------" << std::endl;
  std::cout << "-------" << std::endl;
  std::cout << "ntot - lat_ntot, lon_ntot : (" << location_lookup[0][size-1] << "," <<  location_lookup[1][size-1] << ")" << std::endl;



  // Returnig pointer to storage class for breaking points
  ObservationBreakingPoints* break_points = reader.observation();
  std::cout << "Breaking Points" << std::endl;
  std::cout << break_points->break_time[0] << std::endl;
  std::cout << "-------" << std::endl;
  std::cout << break_points->break_time[1] << std::endl;
  std::cout << break_points->break_time[break_points->number_of_observations()-1] << std::endl;

  std::cout << "-------" << std::endl;
  std::cout << "Number of observations = " << break_points->number_of_observations() << std::endl;
  std::cout << "-------" << std::endl;

  //Stations analysis  
  /*
  std::vector<int32_t> stations = break_points->stations_collection();
  std::vector<int32_t> counter = break_points->stations_count();
  std::cout << "Some of the stations that detected break points" << std::endl;

  for (unsigned i = 0; i <5 ; i++){

  std::cout << "Station index = " << stations[i] << std::endl;
  std::cout << "Number of detections = " << counter[i] << std::endl;

  }
  std::cout << "-------" << std::endl;
  */

  //Filtering break points
  /*
  int filter_station = 4;
  std::vector<int32_t> filtered = break_points->filtered_breakpoints(filter_station);
  std::cout << "Breaking Points filtered by station " << filter_station << std::endl;
  std::cout << "-------" << std::endl;
  for (std::vector<int32_t>::const_iterator it = filtered.begin(); it != filtered.end(); ++it)
    std::cout << *it  << std::endl;
  */

  //Applying rejection policies
  // Caveat: these are purely empirical/arbitrary rejection policies, based on the merged_likelihood attribute: no specific policy has been given by WP1 science owners.
  /*
  break_points->apply_policy(0, 4);

  std::cout << "Breaking Points after rejection" << std::endl;
  std::cout << break_points->break_time[0] << std::endl;
  std::cout << "-------" << std::endl;
  std::cout << break_points->break_time[1] << std::endl;
  std::cout << break_points->break_time[break_points->number_of_observations()-1] << std::endl;
  std::cout << "Number of observations after rejection = " << break_points->number_of_observations() << std::endl;
  std::cout << "-------" << std::endl;
  */

  break_points->~ObservationBreakingPoints();

  reader.close();

}
