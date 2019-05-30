#include <cppunit/extensions/HelperMacros.h>
#include <iostream>
#include <fstream>

#include <ambitious/meshes/meshes.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <string>
#include <iostream>


class TestBlocks  : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( TestBlocks );
  CPPUNIT_TEST( testZeroedWeight );
  CPPUNIT_TEST_SUITE_END();

public:
  
  void testZeroedWeight()
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Scaled to zero (double)",
				 double(0.0),
				 double(-0.5) * double(0.0));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Scaled to zero (bool)",
				 true,
				 double(-0.5) * double(0.0) == double(0.0));
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestBlocks );
