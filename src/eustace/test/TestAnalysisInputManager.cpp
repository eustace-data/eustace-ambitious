#include <cppunit/extensions/HelperMacros.h>
#include <eustace/analysis/inputmanager.h>
#include <iostream>
#include <fstream>

using namespace EUSTACE;
using namespace std;

class TestAnalysisInputManager  : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( TestAnalysisInputManager );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testPathNames );
  CPPUNIT_TEST_SUITE_END();

public:
  
  void testConstructor()
  {
    AnalysisInputManager manager("/some/base/path");
    CPPUNIT_ASSERT_EQUAL(std::string("/some/base/path"), manager.BasePath());
  }

  void testPathNames()
  {
    AnalysisInputManager manager("/a/base/path");
    
    std::string result("somenonsense");

    manager.PathNameLocalCorrelationRange(result, AnalysisInput::source_satellite_ice, AnalysisInput::observable_tmean);
    CPPUNIT_ASSERT_EQUAL(std::string("/a/base/path/surfaceairmodel_ice_Tmean_localcorrelationranges.bin"), result);

    manager.PathNameLocationLookup(result, AnalysisInput::source_satellite_ocean);
    CPPUNIT_ASSERT_EQUAL(std::string("/a/base/path/locationlookup_satellite.bin"), result);

    manager.PathNameMobileLocationLookup(result, AnalysisInput::source_insitu_ocean, AnalysisInput::observable_tmean, 61049);
    CPPUNIT_ASSERT_EQUAL(std::string("/a/base/path/insitu_ocean/2017/insitu_ocean_Tmean_mobilelocations_20170223.bin"), result);

    manager.PathNameObservation(result, AnalysisInput::source_insitu_land, AnalysisInput::observable_tmax, 61048);
    CPPUNIT_ASSERT_EQUAL(std::string("/a/base/path/insitu_land/2017/insitu_land_Tmax_20170222.bin"), result);

    manager.PathNameLocationPointers(result, AnalysisInput::source_satellite_land, AnalysisInput::observable_tmean);
    CPPUNIT_ASSERT_EQUAL(std::string("/a/base/path/surfaceairmodel_land_Tmean_locationpointers.bin"), result);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestAnalysisInputManager );
