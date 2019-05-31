#include "insitu_land_break_points.h"

namespace EUSTACE
{
  
  void BreakNetCDFHelper::GetDimensions(std::vector<size_t>& dimensions, const char* name) const throw (BreakPointsReaderException)
  {
    // ID of variable
    int varid(1);
    nc_inq_varid(ncid, name, &varid);

    // Number of dimensions the variable has
    int ndims(0);
    nc_inq_varndims(ncid, varid, &ndims);
    
    // ID for each dimension
    std::vector<int> dimension_ids(ndims);
    nc_inq_var(ncid, varid, NULL, NULL, NULL, &dimension_ids[0], NULL);
    
    // Count of each dimension
    for (std::vector<int>::const_iterator dimension_id_iterator(dimension_ids.begin()); dimension_id_iterator != dimension_ids.end(); dimension_id_iterator++)
      {
	size_t length(0);
	nc_inq_dimlen(ncid, *dimension_id_iterator, &length);
	dimensions.push_back(length);
      }
  }

  template<class DataT> void BreakNetCDFHelper::ReadData(std::vector<DataT>& contents, const char* name) const throw (BreakPointsReaderException)
  {
    
    int varid(0);
    int status = nc_inq_varid(ncid, name, &varid);
    if (status != NC_NOERR) 
      std::cout << "ID STATUS " << status << std::endl ;

    status = nc_get_var(ncid, varid, &contents[0]);
    if (status != NC_NOERR) 
      std::cout << "VAR STATUS " << status << std::endl ;
  
  }
  
  
  void BreakPointsReader::close() throw(BreakPointsReaderException)
  {
    
    netcdf.~BreakNetCDFHelper();
    
  }
  
  const std::vector<std::vector<float> > BreakPointsReader::ObservationLocationLookup() throw(BreakPointsReaderException)
  {
    
    std::vector<size_t> latitude_dimension;
    std::vector<size_t> longitude_dimension;
    netcdf.GetDimensions(latitude_dimension, "latitude");
    netcdf.GetDimensions(longitude_dimension, "longitude");

    if ((latitude_dimension.size() != 1 ) || (latitude_dimension.size() != 1 ))
      throw BreakPointsReaderException();
    if (latitude_dimension[0] != longitude_dimension[0])
      throw BreakPointsReaderException();

    // Number of points: it has to correspond to the collection of all the stations
    size_t npoints = latitude_dimension[0];

    // Get latitude and longitude
    latitude.reserve(npoints);
    longitude.reserve(npoints);
    netcdf.ReadData<float>(latitude, "latitude");
    netcdf.ReadData<float>(longitude, "longitude");
    
    std::vector<std::vector<float> > location_lookup(2,std::vector<float>(npoints));
    for (unsigned i = 0;i < npoints ; i++){

      // allow for longitudes expressed as [0,360) instead of [-180,180)
      if (longitude[i] >= 180.0)
	longitude[i] -= 360.0;

      location_lookup[0][i]=latitude[i];
      location_lookup[1][i]=longitude[i];
      
    }

    return location_lookup;

  }

  ObservationBreakingPoints* BreakPointsReader::observation() throw(BreakPointsReaderException)
  {

    std::vector<size_t> merged_break_dimension;
    netcdf.GetDimensions(merged_break_dimension, "merged_break_time");
    
    if (merged_break_dimension.size() != 1 )
      throw BreakPointsReaderException();
    
    // Number of points: it has to correspond to the collection of all the brek points
    size_t npoints = merged_break_dimension[0];

    std::vector<int32_t> break_time(npoints);
    std::vector<int32_t> break_station(npoints);
    std::vector<int8_t> temp_break_likelihood(npoints);

    netcdf.ReadData<int32_t>(break_time, time);
    netcdf.ReadData<int32_t>(break_station, station);
    netcdf.ReadData<int8_t>(temp_break_likelihood, likelihood);
    
    std::vector<int> break_likelihood(npoints);
    for(unsigned i = 0; i <npoints ; i++)
      break_likelihood[i]=int(temp_break_likelihood[i]);

    ObservationBreakingPoints* observation = new ObservationBreakingPoints(break_time, break_station, break_likelihood);

    return observation;
  }
}

