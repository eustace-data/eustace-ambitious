#include <cppunit/extensions/HelperMacros.h>
#include <eustace/timeutils/timebase.h>
#include <iostream>
#include <fstream>

using namespace EUSTACE;
using namespace std;

class TestTimeBaseDays  : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( TestTimeBaseDays );
  CPPUNIT_TEST( testNumber );
  CPPUNIT_TEST( testTime );
  CPPUNIT_TEST_SUITE_END();

public:

  void testNumber()
  {
    TimeBaseDays timebase(CalendarDay(1937, 6, 4));
    CPPUNIT_ASSERT_EQUAL(int64_t(0), timebase.Number(CalendarDay(1937, 6, 4)));
    CPPUNIT_ASSERT_EQUAL(int64_t(1), timebase.Number(CalendarDay(1937, 6, 5)));
    CPPUNIT_ASSERT_EQUAL(int64_t(365), timebase.Number(CalendarDay(1938, 6, 4)));
    CPPUNIT_ASSERT_EQUAL(int64_t(730), timebase.Number(CalendarDay(1939, 6, 4)));
    CPPUNIT_ASSERT_EQUAL(int64_t(731), timebase.Number(CalendarDay(1939, 6, 5)));
    CPPUNIT_ASSERT_EQUAL(int64_t(-3), timebase.Number(CalendarDay(1937, 6, 1)));
    CPPUNIT_ASSERT_EQUAL(int64_t(-368), timebase.Number(CalendarDay(1936, 6, 1)));
  }

  void testTime()
  {
    TimeBaseDays timebase(CalendarDay(2012, 3, 7));

    std::tm result = timebase.Time(0);
    CPPUNIT_ASSERT_EQUAL(112, result.tm_year);
    CPPUNIT_ASSERT_EQUAL(2, result.tm_mon);
    CPPUNIT_ASSERT_EQUAL(7, result.tm_mday);
    CPPUNIT_ASSERT_EQUAL(0, result.tm_hour);
    CPPUNIT_ASSERT_EQUAL(0, result.tm_min);
    CPPUNIT_ASSERT_EQUAL(0, result.tm_sec);

    result = timebase.Time(1);
    CPPUNIT_ASSERT_EQUAL(112, result.tm_year);
    CPPUNIT_ASSERT_EQUAL(2, result.tm_mon);
    CPPUNIT_ASSERT_EQUAL(8, result.tm_mday);
    CPPUNIT_ASSERT_EQUAL(0, result.tm_hour);
    CPPUNIT_ASSERT_EQUAL(0, result.tm_min);
    CPPUNIT_ASSERT_EQUAL(0, result.tm_sec);

    result = timebase.Time(368);
    CPPUNIT_ASSERT_EQUAL(113, result.tm_year);
    CPPUNIT_ASSERT_EQUAL(2, result.tm_mon);
    CPPUNIT_ASSERT_EQUAL(10, result.tm_mday);
    CPPUNIT_ASSERT_EQUAL(0, result.tm_hour);
    CPPUNIT_ASSERT_EQUAL(0, result.tm_min);
    CPPUNIT_ASSERT_EQUAL(0, result.tm_sec);
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestTimeBaseDays );
