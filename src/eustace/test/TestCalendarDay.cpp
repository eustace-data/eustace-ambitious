#include <cppunit/extensions/HelperMacros.h>
#include <eustace/timeutils/timebase.h>
#include <iostream>
#include <fstream>

using namespace EUSTACE;
using namespace std;

class TestCalendarDay  : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( TestCalendarDay );
  CPPUNIT_TEST( testConstructFromIntegers );
  CPPUNIT_TEST( testConstructFromString );
  CPPUNIT_TEST( testConvertToString );
  CPPUNIT_TEST_SUITE_END();

public:
 
  void testConstructFromIntegers()
  {
    CalendarDay c(1927, 12, 18);
    CPPUNIT_ASSERT_EQUAL(27, c.tm_year);
    CPPUNIT_ASSERT_EQUAL(11, c.tm_mon);
    CPPUNIT_ASSERT_EQUAL(18, c.tm_mday);
    CPPUNIT_ASSERT_EQUAL(0, c.tm_hour);
    CPPUNIT_ASSERT_EQUAL(0, c.tm_sec);
  }

  void testConstructFromString()
  {
    CalendarDay c(std::string("18540703"));
    CPPUNIT_ASSERT_EQUAL(-46, c.tm_year);
    CPPUNIT_ASSERT_EQUAL(6, c.tm_mon);
    CPPUNIT_ASSERT_EQUAL(3, c.tm_mday);
    CPPUNIT_ASSERT_EQUAL(0, c.tm_hour);
    CPPUNIT_ASSERT_EQUAL(0, c.tm_sec);
  }

  void testConvertToString()
  {
    CalendarDay c(2023, 9, 3);
    std::string str;
    c.Text(str);
    CPPUNIT_ASSERT_EQUAL(std::string("20230903"), str);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestCalendarDay );
