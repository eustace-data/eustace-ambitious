#include <cppunit/extensions/HelperMacros.h>
#include <eustace/analysis/fileio/locationpointers.h>
#include <iostream>
#include <fstream>
#include <memory.h>

using namespace EUSTACE;
using namespace std;

class TestLocationPointers  : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( TestLocationPointers );
  CPPUNIT_TEST( testReadExample );
  CPPUNIT_TEST_SUITE_END();

public:

  void testReadExample()
  {
    // buffer to hold result
    std::vector< std::vector<LocationPointerRecord> > record;

    // scope to create and destroy reader object
    {
      // make input stream
      const char* pathname_locationpointers = EUSTACE_TEST_DATA_DIRECTORY "/example_rawbinary_locationpointers.bin";
      ifstream inputstream(pathname_locationpointers, ios::binary);

      // subset of location pointers to examine
      std::vector<uint64_t> location_id_subset;
      location_id_subset.push_back(1);
      location_id_subset.push_back(2);
      location_id_subset.push_back(3);
      location_id_subset.push_back(4);

      // instantiate reader and read information
      LocationPointersRawBinaryReader reader;
      reader.Read(record, inputstream, location_id_subset);
    }

    // there should be a record for every location id (even ones with no data)
    CPPUNIT_ASSERT_EQUAL(size_t(4), record.size());

    // ID 1 should occur on first day
    CPPUNIT_ASSERT_EQUAL(size_t(1), record[0].size());
    CPPUNIT_ASSERT_EQUAL(int64_t(57499), record[0][0].daynumber);

    // No data for ID 2
    CPPUNIT_ASSERT_EQUAL(size_t(0), record[1].size());

    // And others are on second day
    CPPUNIT_ASSERT_EQUAL(size_t(1), record[2].size());
    CPPUNIT_ASSERT_EQUAL(int64_t(57500), record[2][0].daynumber);
    CPPUNIT_ASSERT_EQUAL(size_t(1), record[3].size());
    CPPUNIT_ASSERT_EQUAL(int64_t(57500), record[3][0].daynumber);

    // check measurement number in observation file
    double measurement(0.0);
    {
      const char* pathname_obs = EUSTACE_TEST_DATA_DIRECTORY "/example_rawbinary_observations_57499.bin";
      ifstream inputstream_obs(pathname_obs, ios::binary);
      inputstream_obs.seekg(record[0][0].fileoffset, ifstream::beg);
      inputstream_obs.read((char*)&measurement, sizeof(measurement));
    }
    CPPUNIT_ASSERT_DOUBLES_EQUAL(288.0, measurement, 0.000001);    
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestLocationPointers );
