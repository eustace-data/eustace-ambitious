#include <cppunit/extensions/HelperMacros.h>
#include <iostream>
#include <fstream>

#include <ambitious/meshes/meshes.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <string>
#include <iostream>


class TestMeshes  : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( TestMeshes );
  CPPUNIT_TEST( testCounterclockwiseTriangleIndexSwap );
  CPPUNIT_TEST( testDateTimeEpochConversion );
  CPPUNIT_TEST( testDateTimeEpochIndex );
  CPPUNIT_TEST( testDateTimeLeapYears );
  CPPUNIT_TEST( testDateTimeInvertibilityWholeToWhole );
  CPPUNIT_TEST( testDateTimeInvertibilityFracToFrac1 );
  CPPUNIT_TEST( testDateTimeInvertibilityFracToFrac2 );
  CPPUNIT_TEST( testTimeMapperMapping1 );
  CPPUNIT_TEST( testTimeMeshNumBasis );
  CPPUNIT_TEST( testTimeMapperMapping2 );
  CPPUNIT_TEST( testConverter );
  CPPUNIT_TEST_SUITE_END();

public:
  
  void testCounterclockwiseTriangleIndexSwap()
  {
    //    for (int i = 0; i<16; i++) {
    //      std::cout << "swap(" << i << ", level) =";
    //      for (int k = 0; k < 5; k++) {
    //	      std::cout << "\t" << triangle_index_swapwise(i,k);
    //      }
    //      std::cout << std::endl;
    //    }
    CPPUNIT_ASSERT_EQUAL_MESSAGE("triangle_index_swapwise(0)",
				 int64_t(0),
				 triangle_index_swapwise(0,4));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("triangle_index_swapwise(1)",
				 int64_t(2),
				 triangle_index_swapwise(1,4));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("triangle_index_swapwise(2)",
				 int64_t(1),
				 triangle_index_swapwise(2,4));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("triangle_index_swapwise(3)",
				 int64_t(3),
				 triangle_index_swapwise(3,4));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("triangle_index_swapwise(4)",
				 int64_t(8+0),
				 triangle_index_swapwise(4,4));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("triangle_index_swapwise(5)",
				 int64_t(8+2),
				 triangle_index_swapwise(5,4));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("triangle_index_swapwise(6)",
				 int64_t(8+1),
				 triangle_index_swapwise(6,4));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("triangle_index_swapwise(7)",
				 int64_t(8+3),
				 triangle_index_swapwise(7,4));
  }


  void testDateTimeLeapYears()
  {
    DateTimeHelper helper(1850);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Days of 1900",
                                 int64_t(365),
                                 helper.days_in_calyear(1900));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Days of 1996",
                                 int64_t(366),
                                 helper.days_in_calyear(1996));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Days of 2000",
                                 int64_t(366),
                                 helper.days_in_calyear(2000));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Days of 2001",
                                 int64_t(365),
                                 helper.days_in_calyear(2001));
  }

  void testDateTimeEpochConversion()
  {
    DateTimeHelper helper0(1850);
    DateTimeHelper helper1(1950);
    int64_t number0 = helper0.get_day(1900, 1, 1);
    int64_t number1 = helper1.get_day(1900, 1, 1);
    int64_t number0f1 = helper0.day_from_other(number1, helper1);
    int64_t number1f0 = helper1.day_from_other(number0, helper0);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Day nr of 1900-01-01",
                                 number0,
                                 number0f1);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Day nr of 1900-01-01",
                                 number1,
                                 number1f0);
  }

  void testDateTimeEpochIndex()
  {
    DateTimeHelper helper(1850);
    int64_t number0 = helper.get_day(1850, 1, 1);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Day nr of 1850-01-01",
                                 int64_t(0),
                                 number0);
    number0 = helper.get_day(1849, 12, 31);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Day nr of 1849-12-31",
                                 int64_t(-1),
                                 number0);
  }
    
  void testDateTimeInvertibilityWholeToWhole()
  {
    DateTimeHelper helper(1850);
    int64_t number0 = helper.get_day(1995, 1, 1);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Year of 1995-01-01",
                                 int64_t(1995),
                                 helper.get_calyear(number0));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Year of 1994-12-31",
                                 int64_t(1994),
                                 helper.get_calyear(number0 - 1));
    double frac_year = helper.fractional_year(number0, 0.0);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Year ID of 1995-01-01",
                                 double(1995 - 1850),
                                 std::floor(frac_year));
    int64_t number1 = helper.whole_day(frac_year);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Day ID number of 1995-01-01",
                                 number0,
                                 number1);
  }

  void testDateTimeInvertibilityFracToFrac1() {
    DateTimeHelper helper(1850);
    double frac_day = 12345.6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Fractional day(year(day)",
                                 frac_day,
                                 helper.fractional_day(helper.fractional_year(frac_day)),
                                 1e-100);
  }
  void testDateTimeInvertibilityFracToFrac2() {
    DateTimeHelper helper(1850);
    double frac_year = 156.7;
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Fractional year(day(year))",
                                 frac_year,
                                 helper.fractional_year(helper.fractional_day(frac_year)),
                                 1e-100);
  }
  void testTimeMapperMapping1() {
    TimeMapper mapper(NULL, 1850);

    CPPUNIT_ASSERT_EQUAL_MESSAGE("Empty time mesh #basis",
                                 int64_t(1),
                                 mapper.num_basis());
  }

  void testTimeMeshNumBasis() {
    TimeMesh mesh1(30, 60, 30, 1850, TimeUnits::Day, TimeBasis::Bspline2, false);
    TimeMesh mesh2(30, 60, 30, 1850, TimeUnits::Day, TimeBasis::Bspline2, true);
    TimeMesh mesh3(30, 60, 30, 1850, TimeUnits::Year, TimeBasis::Bspline2, false);
    TimeMesh mesh4(30, 60, 30, 1850, TimeUnits::Year, TimeBasis::Bspline2, true);
    TimeMesh mesh5(0.0, 1.0, 30, 1850, TimeUnits::Season, TimeBasis::Bspline2, false);
    TimeMesh mesh6(0.0, 1.0, 30, 1850, TimeUnits::Season, TimeBasis::Bspline2, true);
    TimeMesh mesh7(0.0, 1.0, 30, 1850, TimeUnits::Season, TimeBasis::Harmonic, false);
    TimeMesh mesh8(0.0, 1.0, 30, 1850, TimeUnits::Season, TimeBasis::Harmonic, true);

    CPPUNIT_ASSERT_EQUAL_MESSAGE("TimeMesh(30, 60, 30, 1850, Day, Bspline2, false).num_basis",
                                 int64_t(32),
                                 mesh1.num_basis());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("TimeMesh(30, 60, 30, 1850, Day, Bspline2, true).num_basis",
                                 int64_t(30),
                                 mesh2.num_basis());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("TimeMesh(30, 60, 30, 1850, Year, Bspline2, false).num_basis",
                                 int64_t(32),
                                 mesh3.num_basis());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("TimeMesh(30, 60, 30, 1850, Year, Bspline2, true).num_basis",
                                 int64_t(30),
                                 mesh4.num_basis());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("TimeMesh(0.0, 1.0, 30, 1850, Season, Bspline2, false).num_basis",
                                 int64_t(32),
                                 mesh5.num_basis());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("TimeMesh(0.0, 1.0, 30, 1850, Season, Bspline2, true).num_basis",
                                 int64_t(30),
                                 mesh6.num_basis());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("TimeMesh(0.0, 1.0, 30, 1850, Season, Harmonic, false).num_basis",
                                 int64_t(61),
                                 mesh7.num_basis());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("TimeMesh(0.0, 1.0, 30, 1850, Season, Harmonic, true).num_basis",
                                 int64_t(61),
                                 mesh8.num_basis());
  }

  void testTimeMapperMapping2() {
    TimeMesh mesh(30, 60, 30, 1850, TimeUnits::Day, TimeBasis::Bspline2, false);
    TimeMapper t_mapper(&mesh, 1850);
    SpaceTimeMapper mapper(&t_mapper, NULL, NULL, NULL, 1850);
    std::vector<Eigen::Triplet<double>> output;

    CPPUNIT_ASSERT_EQUAL_MESSAGE("Time, Daily, 30 days",
                                 int64_t(32),
                                 t_mapper.num_basis());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("SpaceTime, Daily, 30 days",
                                 int64_t(32),
                                 mapper.num_basis());

    SpaceTimeMapper::Point point;
    int64_t num;

    point = SpaceTimeMapper::Point(0.0, 0.0, 0, false);
    SpaceTimeMapper::offset_type offset(0, 0);
    num = mapper.mapping(point, false, offset, false, NULL, output);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("SpaceTime, Map day 0, num",
                                 int64_t(0),
                                 num);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("SpaceTime, Map day 0, num=size",
                                 uint64_t(num),
                                 output.size());
    output.clear();

    point = SpaceTimeMapper::Point(0.0, 0.0, 30, false);
    num = mapper.mapping(point, false, offset, false, NULL, output);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("SpaceTime, Map day 30, num",
                                 int64_t(3),
                                 num);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("SpaceTime, Map day 30, num=size",
                                 uint64_t(num),
                                 output.size());
    output.clear();

    point = SpaceTimeMapper::Point(0.0, 0.0, 30, false);
    std::vector<double> weights{2.0, 3.0};
    num = mapper.mapping(point, false, offset, false, &weights, output);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("SpaceTime, Map day 30, 2 weights, num",
                                 int64_t(3 * 2),
                                 num);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("SpaceTime, Map day 30, 2 weights, num=size",
                                 uint64_t(num),
                                 output.size());
    output.clear();
}

  typedef fmesh::Point FP;
  typedef std::pair<double, double> LL;
  void testConverter()
  {
    FP fp;
    LL ll;

    Converter::latlong_to_euclidean(LL(0.0, 0.0), fp);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("(0,0) to Euclidean, x",
                                         1.0,
                                         double(fp[0]),
                                         1e-100);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("(0,0) to Euclidean, y",
                                         0.0,
                                         double(fp[1]),
                                         1e-100);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("(0,0) to Euclidean, z",
                                         0.0,
                                         double(fp[2]),
                                         1e-100);

    Converter::euclidean_to_latlong(FP(1.0, 0.0, 0.0), ll);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("(1,0,0) to latlong, lat",
                                         0.0,
                                         double(ll.first),
                                         1e-100);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("(1,0,0) to latlong, long",
                                         0.0,
                                         double(ll.second),
                                         1e-100);
    


    Converter::latlong_to_euclidean(LL(0.0, 90.0), fp);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("(0,90) to Euclidean, x",
                                         0.0,
                                         double(fp[0]),
                                         1e-16);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("(0,90) to Euclidean, y",
                                         1.0,
                                         double(fp[1]),
                                         1e-16);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("(0,90) to Euclidean, z",
                                         0.0,
                                         double(fp[2]),
                                         1e-16);

    Converter::euclidean_to_latlong(FP(0.0, 1.0, 0.0), ll);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("(0,1,0) to latlong, lat",
                                         0.0,
                                         double(ll.first),
                                         1e-100);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("(0,1,0) to latlong, long",
                                         90.0,
                                         double(ll.second),
                                         1e-100);
    

    


    Converter::latlong_to_euclidean(LL(90.0, 0.0), fp);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("(90,0) to Euclidean, x",
                                         0.0,
                                         double(fp[0]),
                                         1e-16);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("(90,0) to Euclidean, y",
                                         0.0,
                                         double(fp[1]),
                                         1e-16);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("(90,0) to Euclidean, z",
                                         1.0,
                                         double(fp[2]),
                                         1e-16);
    
    Converter::euclidean_to_latlong(FP(0.0, 0.0, 1.0), ll);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("(0,0,1) to latlong, lat",
                                         90.0,
                                         double(ll.first),
                                         1e-100);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("(0,0,1) to latlong, long",
                                         0.0,
                                         double(ll.second),
                                         1e-100);
}



};

CPPUNIT_TEST_SUITE_REGISTRATION( TestMeshes );
