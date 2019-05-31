#include <cppunit/extensions/HelperMacros.h>
#include <eustace/analysis/break_points/insitu_land_break_points.h>
#include <iostream>
#include <fstream>
#include <memory.h>
#include <stdlib.h>
#include <typeinfo>
#include <vector>

using namespace EUSTACE;

// rawbinary3 uses R 001124 data and R001127 status
// const char* FILENAME = "/gws/nopw/j04/eustace/data/internal/D1.7/daily/eustace_stations_global_R001127_daily_status.nc";
const char* FILENAME = "/home/flindgre/EUSTACE/svn/platform/research/uedin-sandbox/eustace_stations_global_R001127_daily_status.nc";

class TestBreakPointsReader  : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( TestBreakPointsReader );
  CPPUNIT_TEST( testReadBreakPoints );
  CPPUNIT_TEST( testObservationLocationLookup );
  CPPUNIT_TEST( testReturnObservation );
  CPPUNIT_TEST_SUITE_END();

public:
 
  void testReadBreakPoints()
  {

    const char* wrong_name = "gababgfobgo";
    CPPUNIT_ASSERT_THROW(BreakPointsReader reader(wrong_name),BreakPointsReaderException);

  }

  void testObservationLocationLookup()
  {


    BreakPointsReader reader(FILENAME);
    std::vector<std::vector<float> > location_lookup = reader.ObservationLocationLookup();

    unsigned indices[] = { 0, 10, 35207};
    float expected_latitudes[] = { 48.4000, 48.8167, -22.2170};
    float expected_longitudes[] = { -123.4833,  -124.1333, 30.0000};

    for (unsigned i = 0; i < 3; i++)
      {
	
	CPPUNIT_ASSERT_EQUAL(expected_latitudes[i], location_lookup[0][indices[i]]);
	CPPUNIT_ASSERT_EQUAL(expected_longitudes[i], location_lookup[1][indices[i]]);
   
      }
    
    reader.close();
    
  }

  void testReturnObservation()
  {

    BreakPointsReader reader(FILENAME);
    ObservationBreakingPoints* break_points = reader.observation();
    //    size_t shape = 155943; // old data
    size_t shape = 130950;
    CPPUNIT_ASSERT_EQUAL(shape, break_points->break_time.size());
    CPPUNIT_ASSERT_EQUAL(shape, break_points->break_station.size());
    CPPUNIT_ASSERT_EQUAL(shape, break_points->break_likelihood.size());

    unsigned array_indices[] = {24, 119};
    int32_t expected_break_time[] = {37620, 48212};
    int32_t expected_break_station[] = {9, 39};
    int32_t expected_break_likelihood[] = {31, 14};

    for (unsigned i = 0; i <2 ; i++){

      CPPUNIT_ASSERT_EQUAL(expected_break_time[i], break_points->break_time[array_indices[i]]);
      CPPUNIT_ASSERT_EQUAL(expected_break_station[i], break_points->break_station[array_indices[i]]);
      CPPUNIT_ASSERT_EQUAL(expected_break_likelihood[i], break_points->break_likelihood[array_indices[i]]);

    }
 
    break_points->~ObservationBreakingPoints();
    reader.close();
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestBreakPointsReader );

 
