#include <cppunit/extensions/HelperMacros.h>
#include <iostream>
#include <fstream>
#include <ambitious/poq/poq.hpp>


class TestPoq  : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( TestPoq );
  CPPUNIT_TEST( testPoqCdf );
  CPPUNIT_TEST_SUITE_END();

public:
  
  void testPoqCdf()
  {
    poq_distribution<double>::Vector theta(5);
    theta << -0.5, 0.1, 0.2, 0.3, 0.4;
    poq_distribution<double> poq(theta);

    double p0 = 0.6;
    double q = poq.quantile(p0);
    double p = poq.cdf(q);

    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("cdf(quantile(p)) - p", 0.0, p-p0, 1e-15);
  }


};

CPPUNIT_TEST_SUITE_REGISTRATION( TestPoq );
