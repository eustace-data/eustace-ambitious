#include <eustace/analysis/break_points/break_points.h>
#include <cppunit/extensions/HelperMacros.h>
#include <iostream>
#include <fstream>
#include <stdint.h>
#include <vector>
#include <memory>

using namespace EUSTACE;
using namespace std;

class TestObservationBreakingPoints : public CppUnit::TestCase
{
  
public:
  CPPUNIT_TEST_SUITE( TestObservationBreakingPoints );
  CPPUNIT_TEST( testConstructorAndAttributes );
  CPPUNIT_TEST( testStationsRelatedMethods );
  CPPUNIT_TEST( testHardCutOff );
  CPPUNIT_TEST( testLaplaceKernel );
  CPPUNIT_TEST_SUITE_END();
  
  void testConstructorAndAttributes()
  {
    
    int size = 5;
    int32_t tmp0[] = {54786, 56978, 57343, 57708, 58074};
    std::vector<int32_t> break_time( tmp0, tmp0+size ); 
    
    int32_t tmp1[] = {4, 4, 4, 4, 4};
    std::vector<int32_t> break_station( tmp1, tmp1+size ); 
    
    int tmp2[] = {1, 1, 1, 1, 1};
    std::vector<int> break_likelihood( tmp2, tmp2+size ); 
    
    ObservationBreakingPoints A(break_time, break_station, break_likelihood);
    CPPUNIT_ASSERT_EQUAL(A.number_of_observations(),5);
    
    for (unsigned i = 0; i <  5; i++)
      {
	
	CPPUNIT_ASSERT_EQUAL(A.break_time[i], tmp0[i]);
	CPPUNIT_ASSERT_EQUAL(A.break_station[i], tmp1[i]);
	CPPUNIT_ASSERT_EQUAL(A.break_likelihood[i], tmp2[i]);
	
      }
  }
  
  void  testStationsRelatedMethods()
  {

    int size = 8;
    int32_t tmp0[] = {54786, 56978, 57343, 57708, 58074, 57343, 57708, 58074};
    std::vector<int32_t> break_time( tmp0, tmp0+size ); 
    
    int32_t tmp1[] = {4, 3, 3, 3, 2, 3, 0, 0};
    std::vector<int32_t> break_station( tmp1, tmp1+size ); 
    
    int tmp2[] = {1, 1, 1, 1, 1, 1, 1, 1};
    std::vector<int> break_likelihood( tmp2, tmp2+size ); 
    
    ObservationBreakingPoints A(break_time, break_station, break_likelihood);
    CPPUNIT_ASSERT_EQUAL(A.number_of_stations(),4);

    int32_t expected_collection[] = {0, 2, 3, 4}; 
    int32_t expected_count[] = {2, 1, 4, 1}; 
    int32_t expected_filtered[] = {56978, 57343, 57708, 57343}; 
    int32_t filter_station = 3;

    size_t vec_size = 4;
    CPPUNIT_ASSERT_EQUAL(A.stations_collection().size(),vec_size);
    CPPUNIT_ASSERT_EQUAL(A.stations_count().size(),vec_size);
    CPPUNIT_ASSERT_EQUAL(A.filtered_breakpoints(filter_station).size(),vec_size);

    for (unsigned i = 0; i <  4; i++)
      {
	
	CPPUNIT_ASSERT_EQUAL(A.stations_collection()[i], expected_collection[i]);
       	CPPUNIT_ASSERT_EQUAL(A.stations_count()[i], expected_count[i]);
	CPPUNIT_ASSERT_EQUAL(A.filtered_breakpoints(filter_station)[i], expected_filtered[i]);
   
      }
  }

  void testHardCutOff()
  {

    int size = 8;
    int32_t tmp0[] = {54786, 56978, 57343, 57708, 58074, 57343, 57708, 58074};
    std::vector<int32_t> break_time( tmp0, tmp0+size ); 
    
    int32_t tmp1[] = {4, 3, 3, 3, 2, 0, 0, 0};
    std::vector<int32_t> break_station( tmp1, tmp1+size ); 
    
    int tmp2[] = {1, 100, 112, 13, 10, 21, 2, 3};
    std::vector<int> break_likelihood( tmp2, tmp2+size ); 
    
    ObservationBreakingPoints A(break_time, break_station, break_likelihood);
    CPPUNIT_ASSERT_THROW(A.apply_policy(4),std::runtime_error);
    
    A.apply_policy(0,13);
    
    int32_t expected_break_time[] = {56978, 57343, 57708, 57343};
    int32_t expected_break_station[] = {3, 3, 3, 0};
    int32_t expected_break_likelihood[] = {100, 112, 13, 21};
    
    for (unsigned i = 0; i <  4; i++)
      {
	
	CPPUNIT_ASSERT_EQUAL(A.break_time[i], expected_break_time[i]);
	CPPUNIT_ASSERT_EQUAL(A.break_station[i], expected_break_station[i]);
	CPPUNIT_ASSERT_EQUAL(A.break_likelihood[i], expected_break_likelihood[i]);
	
      }
  }

  void testLaplaceKernel()
  {

    int size = 8;
    int32_t tmp0[] = {54786, 56978, 57343, 57708, 58074, 57343, 57708, 58074};
    std::vector<int32_t> break_time( tmp0, tmp0+size ); 
    
    int32_t tmp1[] = {4, 3, 3, 3, 2, 0, 0, 0};
    std::vector<int32_t> break_station( tmp1, tmp1+size ); 
    
    int tmp2[] = {1, 100, 112, 13, 10, 21, 2, 3};
    std::vector<int> break_likelihood( tmp2, tmp2+size ); 
    
    ObservationBreakingPoints A(break_time, break_station, break_likelihood);
    A.apply_policy(1,.5,.01);
    
    int32_t expected_break_time[] = {56978, 57343};
    int32_t expected_break_station[] = {3, 3};
    int32_t expected_break_likelihood[] = {100, 112};
    
    for (unsigned i = 0; i <  2; i++)
      {
	
	CPPUNIT_ASSERT_EQUAL(A.break_time[i], expected_break_time[i]);
	CPPUNIT_ASSERT_EQUAL(A.break_station[i], expected_break_station[i]);
	CPPUNIT_ASSERT_EQUAL(A.break_likelihood[i], expected_break_likelihood[i]);
	
      }
  }

  
};
CPPUNIT_TEST_SUITE_REGISTRATION( TestObservationBreakingPoints );

