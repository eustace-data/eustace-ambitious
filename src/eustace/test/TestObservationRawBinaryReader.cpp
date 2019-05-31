#include <cppunit/extensions/HelperMacros.h>
#include <eustace/analysis/fileio/rawbinary.h>
#include <iostream>
#include <fstream>
#include <memory.h>

using namespace EUSTACE;
using namespace std;

class TestRawBinaryObservationReader  : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( TestRawBinaryObservationReader );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testReadExample );
  CPPUNIT_TEST_SUITE_END();

public:
  void testConstructor()
  {
    ObservationRawBinaryReader result(7);
  }

  void testReadExample()
  {
    ObservationRawBinaryReader reader(2);
    const char* pathname = EUSTACE_TEST_DATA_DIRECTORY "/example_rawbinary_observations_57500.bin";
    ifstream inputstream(pathname, ios::binary);
    ObservationEntries result;
    reader.Read(result, inputstream);

    // check expected UUID
    uint8_t uuid_bytes[] = { 0x49, 0x96, 0x79, 0x7e, 0xf0, 0xdc, 0x44, 0x95, 0xb9, 0x1d, 0x43, 0x45, 0x18, 0x9f, 0xcb, 0x71 };
    CPPUNIT_ASSERT_EQUAL(0, memcmp(uuid_bytes, result.LocationLookupID().Bytes(), 16));

    //      location = numpy.array([ 0, 1, 3, 4 ], numpy.uint64)
    //      measurements = numpy.array([ 275.0, 288.0, 204.0, 310.0 ], numpy.float32)
    //      uncorrelatederror = numpy.array([5.0, 5.0, 7.0, 8.0], numpy.float32)
    //      locallycorrelatederror = numpy.array([ [ 2.3, 3.3,  1.2, 8.9 ], [ 0.222, 7.6, 8.888, 9.999 ] ], numpy.float32)

    // check contents
    CPPUNIT_ASSERT_EQUAL(size_t(3), result.NumEntries());
    
    ObservationIterator obs(result);
    CPPUNIT_ASSERT_EQUAL(uint64_t(0), obs.Identifier());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(275.0, obs.Measurement(), 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, obs.UncertaintyUncorrelated(), 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.3, obs.UncertaintyLocallyCorrelated(0), 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.222, obs.UncertaintyLocallyCorrelated(1), 1E-6);
    obs.Next();
    CPPUNIT_ASSERT_EQUAL(uint64_t(3), obs.Identifier());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(204.0, obs.Measurement(), 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(7.0, obs.UncertaintyUncorrelated(), 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.2, obs.UncertaintyLocallyCorrelated(0), 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(8.888, obs.UncertaintyLocallyCorrelated(1), 1E-6);
    obs.Next();
    CPPUNIT_ASSERT_EQUAL(uint64_t(4), obs.Identifier());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(310.0, obs.Measurement(), 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(8.0, obs.UncertaintyUncorrelated(), 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(8.9, obs.UncertaintyLocallyCorrelated(0), 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(9.999, obs.UncertaintyLocallyCorrelated(1), 1E-6);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestRawBinaryObservationReader );
