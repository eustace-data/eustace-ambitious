#include <cppunit/extensions/HelperMacros.h>
#include <eustace/analysis/inputmanager.h>
#include <iostream>
#include <fstream>

using namespace EUSTACE;
using namespace std;

class TestAnalysisInput  : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( TestAnalysisInput );
  CPPUNIT_TEST( testConstructorAndAttributes );
  CPPUNIT_TEST( testSourceNames );
  CPPUNIT_TEST( testObservableNames );
  CPPUNIT_TEST_SUITE_END();

public:
  
  void testConstructorAndAttributes()
  {
    AnalysisInput a(AnalysisInput::source_insitu_ocean, AnalysisInput::observable_tmax, 50001, 60002);
    CPPUNIT_ASSERT_EQUAL(AnalysisInput::source_insitu_ocean, a.Source());
    CPPUNIT_ASSERT_EQUAL(AnalysisInput::observable_tmax, a.Observable());
    CPPUNIT_ASSERT_EQUAL(int64_t(50001), a.StartDay());
    CPPUNIT_ASSERT_EQUAL(int64_t(60002), a.EndDay());
  }

  void testSourceNames()
  {
    CPPUNIT_ASSERT_EQUAL(std::string("surfaceairmodel_land"), 
			 std::string(AnalysisInput::sourcename[AnalysisInput::source_satellite_land]));
    CPPUNIT_ASSERT_EQUAL(std::string("surfaceairmodel_ocean"), 
			 std::string(AnalysisInput::sourcename[AnalysisInput::source_satellite_ocean]));
    CPPUNIT_ASSERT_EQUAL(std::string("surfaceairmodel_ice"), 
			 std::string(AnalysisInput::sourcename[AnalysisInput::source_satellite_ice]));
    CPPUNIT_ASSERT_EQUAL(std::string("surfaceairmodel_lake"), 
			 std::string(AnalysisInput::sourcename[AnalysisInput::source_satellite_lake]));
    CPPUNIT_ASSERT_EQUAL(std::string("insitu_land"), 
			 std::string(AnalysisInput::sourcename[AnalysisInput::source_insitu_land]));
    CPPUNIT_ASSERT_EQUAL(std::string("insitu_ocean"), 
			 std::string(AnalysisInput::sourcename[AnalysisInput::source_insitu_ocean]));
  }

  void testObservableNames()
  {
    CPPUNIT_ASSERT_EQUAL(std::string("Tmin"), std::string(AnalysisInput::observablename[AnalysisInput::observable_tmin]));
    CPPUNIT_ASSERT_EQUAL(std::string("Tmax"), std::string(AnalysisInput::observablename[AnalysisInput::observable_tmax]));
    CPPUNIT_ASSERT_EQUAL(std::string("Tmean"), std::string(AnalysisInput::observablename[AnalysisInput::observable_tmean]));
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestAnalysisInput );
