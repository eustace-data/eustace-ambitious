#include <cppunit/extensions/HelperMacros.h>
#include <eustace/analysis/fileio/rawbinary.h>
#include <iostream>
#include <fstream>

using namespace EUSTACE;
using namespace std;

class TestLocalCorrelationRangeRawBinaryReader  : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( TestLocalCorrelationRangeRawBinaryReader );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testReadExample );
  CPPUNIT_TEST_SUITE_END();

public:
  void testConstructor()
  {
    LocalCorrelationRangeRawBinaryReader result;
  }

  void testReadExample()
  {
    LocalCorrelationRangeRawBinaryReader reader;
    const char* pathname = EUSTACE_TEST_DATA_DIRECTORY "/example_rawbinary_localcorrelationrange.bin";
    ifstream inputstream(pathname, ios::binary);
    vector<double> ranges;
    reader.Read(ranges, inputstream);
    CPPUNIT_ASSERT_EQUAL(size_t(2), ranges.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(6.6, ranges[0], 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.5, ranges[1], 1E-6);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestLocalCorrelationRangeRawBinaryReader );
