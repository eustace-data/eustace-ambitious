#ifndef _EUSTACE_ANALYSIS_BREAK_POINTS_INSITU_LAND_BREAK_POINTS_
#define _EUSTACE_ANALYSIS_BREAK_POINTS_INSITU_LAND_BREAK_POINTS_

#include<eustace/analysis/break_points/break_points.h>
#include<iostream>
#include <math.h>
#include <netcdf.h>
#include <stdint.h>
#include <stdint.h>
#include <vector>

namespace EUSTACE
{
  class BreakPointsReaderException
  {
  };
  
  class BreakNetCDFHelper
  {

  public:
    BreakNetCDFHelper(const char* filename) throw (BreakPointsReaderException):
    ncid(0)
      {
	// Open it and set member file id                                                                                                                                                              
	int errorcode = nc_open(filename, NC_NOWRITE, &ncid);
	if (errorcode)
	  {
	    throw BreakPointsReaderException();
	  }
      }
    
  public:
    virtual ~BreakNetCDFHelper() 
      {
    if (ncid)
    {
      nc_close(ncid);
      ncid = 0;
    }
  }
  
  public:
    void GetDimensions(std::vector<size_t>& dimensions, const char* name) const throw (BreakPointsReaderException);
    template<class DataT> void ReadData(std::vector<DataT>& contents, const char* name) const throw (BreakPointsReaderException);
    
  private:
    int ncid;
  };
  
  class BreakPointsReader
  {
    
  public:  
  BreakPointsReader(const char* filename) :
    netcdf(filename)
    {
      
      variable= "merged_break";
      time = "merged_break_time";
      station = "merged_break_station";
      likelihood = "merged_break_likelihood";
      
    }
    
  public:
    virtual ~BreakPointsReader()
      {
      }
    
  public:
    void close() throw(BreakPointsReaderException);
    const std::vector<std::vector<float> > ObservationLocationLookup() throw(BreakPointsReaderException);
    ObservationBreakingPoints* observation() throw(BreakPointsReaderException);
    
  private:    
    const char* variable;
    const char* time;
    const char* station;
    const char* likelihood;
   
    // Helper for reading netcdf file
    BreakNetCDFHelper netcdf;
    
    // For computing location lookup
    std::vector<float> latitude;
    std::vector<float> longitude;

  };
}

#endif // #ifdef _EUSTACE_ANALYSIS_BREAK_POINTS_INSITU_LAND_BREAK_POINTS_
