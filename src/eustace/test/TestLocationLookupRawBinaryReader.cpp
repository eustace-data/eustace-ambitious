#include <cppunit/extensions/HelperMacros.h>
#include <eustace/analysis/fileio/rawbinary.h>
#include <iostream>
#include <fstream>
#include <memory.h>

using namespace EUSTACE;
using namespace std;

class TestLocationLookupRawBinaryReader  : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( TestLocationLookupRawBinaryReader );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testReadExample );
  CPPUNIT_TEST_SUITE_END();

public:
  void testConstructor()
  {
    LocationLookupRawBinaryReader result;
  }

  void testReadExample()
  {
    LocationLookupRawBinaryReader reader;
    const char* pathname = EUSTACE_TEST_DATA_DIRECTORY "/example_rawbinary_locationlookup.bin";
    ifstream inputstream(pathname, ios::binary);
    LocationLookup locationlookup;
    reader.Read(locationlookup, inputstream);

    // check expected UUID
    uint8_t uuid_bytes[] = { 0x49, 0x96, 0x79, 0x7e, 0xf0, 0xdc, 0x44, 0x95, 0xb9, 0x1d, 0x43, 0x45, 0x18, 0x9f, 0xcb, 0x71 };
    CPPUNIT_ASSERT_EQUAL(0, memcmp(uuid_bytes, locationlookup.ID().Bytes(), 16));

    // check contents
    CPPUNIT_ASSERT_EQUAL(size_t(5), locationlookup.Lookup().size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 45.0, locationlookup.Lookup()[0].latitude, 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( -8.9, locationlookup.Lookup()[0].longitude, 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 38.0, locationlookup.Lookup()[1].latitude, 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 22.1, locationlookup.Lookup()[1].longitude, 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 42.0, locationlookup.Lookup()[2].latitude, 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(  0.0, locationlookup.Lookup()[2].longitude, 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 58.0, locationlookup.Lookup()[3].latitude, 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( -5.0, locationlookup.Lookup()[3].longitude, 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 41.0, locationlookup.Lookup()[4].latitude, 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 10.0, locationlookup.Lookup()[4].longitude, 1E-6);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestLocationLookupRawBinaryReader );
