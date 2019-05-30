#include <cppunit/extensions/HelperMacros.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <ambitious/timer/timer.hpp>


class TestTimer  : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( TestTimer );
  CPPUNIT_TEST( testTimerConsistency );
  CPPUNIT_TEST_SUITE_END();

public:
  
  void testTimerConsistency()
  {
    Timer timer;

    std::vector<unsigned int> delays(4);
    delays[0] = 1000000;
    delays[1] = 2000000;
    delays[2] = 5000000;
    delays[3] =  500000;

    for (unsigned int i=0; i < delays.size(); ++i) {
      timer.tic();
      usleep(delays[i]);
      timer.toc();

      // 2e-3 chosen to increase likelihood of OK test under valgrind
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("timer.count - delay",
					   0.0,
					   timer.count() - double(delays[i])/1e6,
					   2e-3);
    }
  }


};

CPPUNIT_TEST_SUITE_REGISTRATION( TestTimer );
