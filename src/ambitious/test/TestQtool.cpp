#include <cppunit/extensions/HelperMacros.h>
#include <iostream>
#include <fstream>
#include <ambitious/qtool/qtool.hpp>


#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <string>
#include <iostream>


typedef double MyScalar;
using Vector = Eigen::Matrix<MyScalar, Eigen::Dynamic, 1>;
using Matrix = Eigen::Matrix<MyScalar, Eigen::Dynamic, Eigen::Dynamic>;
using SparseMatrix = Eigen::SparseMatrix<MyScalar, Eigen::ColMajor>;

int pow2(unsigned int exponent)
{
  if (exponent == 0) {
    return 1;
  } else {
    return 2 << (exponent-1);
  }
}



class TestQtool  : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( TestQtool );
  CPPUNIT_TEST( testPartialInversion );
  CPPUNIT_TEST_SUITE_END();

public:
  
  void testPartialInversion()
  {
    Eigen::SparseMatrix<double> Q(3,3);
    Q.insert(0,0) = 2.0/4;
    Q.insert(0,1) = -1.0/4;
    Q.insert(1,0) = -1.0/4;
    Q.insert(1,1) = 2.0/4;
    Q.insert(1,2) = -1.0/4;
    Q.insert(2,1) = -1.0/4;
    Q.insert(2,2) = 2.0/4;

    Eigen::SparseMatrix<double> S(3,3);
    S.insert(0,0) = 3.0;
    S.insert(0,1) = 2.0;
    S.insert(1,0) = 2.0;
    S.insert(1,1) = 4.0;
    S.insert(1,2) = 2.0;
    S.insert(2,1) = 2.0;
    S.insert(2,2) = 3.0;

    QTool<double> Q_tool;
    Q_tool.Q(Q);
    Q_tool.analyzePattern();
    Q_tool.factorize();
    Q_tool.partialInversion();
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("norm(Qtool.S - S)",
					 0.0,
					 (Q_tool.S() - S).norm(),
					 2.5e-15);

    Q_tool.Q(Q);
    Q_tool.partialInversion();
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("norm(Qtool.S - S) after automatic factorisation",
					 0.0,
					 (Q_tool.S() - S).norm(),
					 2.5e-15);
  }


};

CPPUNIT_TEST_SUITE_REGISTRATION( TestQtool );
