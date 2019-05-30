#include <cppunit/extensions/HelperMacros.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <ambitious/ncdf/ncdf.hpp>
#include <ambitious/to_string/to_string.hpp>

#define CPPUNIT_ASSERT_NEQ_MESSAGE(msg, avoid, actual) \
  CPPUNIT_ASSERT_MESSAGE(std::string( \
      std::string("not equal assertion failed\n") +         \
      "- Avoid this  : " + to_string(avoid) + "\n" +    \
      "- Actual value: " + to_string(actual) + "\n" +                   \
      "- " + msg).c_str(),                                              \
                         actual != avoid)
#define CPPUNIT_ASSERT_GEQ_MESSAGE(msg, thelimit, actual) \
  CPPUNIT_ASSERT_MESSAGE(std::string( \
      std::string("greater than or equal assertion failed\n") +         \
      "- Lower limit : " + to_string(thelimit) + "\n" +    \
      "- Actual value: " + to_string(actual) + "\n" +      \
      "- " + msg).c_str(),                                              \
                         actual >= thelimit)
#define CPPUNIT_ASSERT_LEQ_MESSAGE(msg, thelimit, actual) \
  CPPUNIT_ASSERT_MESSAGE(std::string( \
      std::string("less than or equal assertion failed\n") +    \
      "- Upper limit : " + to_string(thelimit) + "\n" +         \
      "- Actual value: " + to_string(actual) + "\n" +           \
      "- " + msg).c_str(),                                      \
                         actual <= thelimit)


template <class FFF>
class TestNCDF_EigenScaling  : public CppUnit::TestCase
{
 public:
  CPPUNIT_TEST_SUITE( TestNCDF_EigenScaling );
  CPPUNIT_TEST( testFromLimits<short> );
  CPPUNIT_TEST( testFromLimits<ushort> );
  CPPUNIT_TEST( testFromLimits<int> );
  CPPUNIT_TEST( testFromLimits<unsigned int> );
  CPPUNIT_TEST( testFromLimits<long> );
  CPPUNIT_TEST( testFromLimits<unsigned long> );
  CPPUNIT_TEST( testFromLimits<long long> );
  CPPUNIT_TEST( testFromLimits<unsigned long long> );

  CPPUNIT_TEST( testInversion<short> );
  CPPUNIT_TEST( testInversion<ushort> );
  CPPUNIT_TEST( testInversion<int> );
  CPPUNIT_TEST( testInversion<unsigned int> );
  CPPUNIT_TEST( testInversion<long> );
  CPPUNIT_TEST( testInversion<unsigned long> );
  CPPUNIT_TEST( testInversion<long long> );
  CPPUNIT_TEST( testInversion<unsigned long long> );
  CPPUNIT_TEST_SUITE_END();
  
 public:

  template <class TTT>
  void check_extremes(const std::string& msg,
                      TTT lower, TTT upper, TTT middle,
                      TTT mini, TTT maxi)
  {
    CPPUNIT_ASSERT_GEQ_MESSAGE(msg + ", lower", mini, lower);
    CPPUNIT_ASSERT_LEQ_MESSAGE(msg + ", upper", maxi, upper);
    CPPUNIT_ASSERT_LEQ_MESSAGE(msg + ", lower should be small", middle, lower);
    CPPUNIT_ASSERT_GEQ_MESSAGE(msg + ", upper should be large", middle, upper);
  }

  template <class TTT>
  void check_not_missing(const std::string& msg,
                         TTT lower, TTT upper, TTT middle,
                         TTT missing)
  {
    CPPUNIT_ASSERT_NEQ_MESSAGE(msg + ", lower shouldn't be missing",
                               missing,
                               lower);
    CPPUNIT_ASSERT_NEQ_MESSAGE(msg + ", upper shouldn't be missing",
                               missing,
                               upper);
    CPPUNIT_ASSERT_NEQ_MESSAGE(msg + ", middle shouldn't be missing",
                               missing,
                               middle);
  }


  template <class TTT>
  void testFromLimits()
  {
    //    if (sizeof(FFF) <= sizeof(TTT)) {
    //  return;
    //}
    
    //    Eigen::Matrix<double, Eigen::Dynamic, 1> vec;
    //    std::pair<double, FFF> lim{-100.0, 100.0};
    std::vector<std::pair<FFF, FFF> > lims;
    lims.push_back(std::pair<FFF, FFF>(-10.0, 10.0));
    lims.push_back(std::pair<FFF, FFF>(-50.0, 50.0));
    lims.push_back(std::pair<FFF, FFF>(-50.0+273.15, 50.0+273.15));
    lims.push_back(std::pair<FFF, FFF>(-100.0, 100.0));
    lims.push_back(std::pair<FFF, FFF>(-100.0+273.15, 100.0+273.15));
    
    for (typename std::vector<std::pair<FFF, FFF> >::iterator it =  lims.begin();
         it != lims.end();
         ++it) {
      std::pair<FFF, FFF> lim = *it;
      
      CompressedVector<FFF, TTT> cvec;
      {
        cvec.set_scaling(lim, false);

        TTT lower = cvec.compress(lim.first);
        TTT upper = cvec.compress(lim.second);
        TTT middle = cvec.compress((lim.first + lim.second) / 2);

        check_extremes("no missing",
                       lower, upper, middle,
                       (TTT)std::numeric_limits<TTT>::min(),
                       (TTT)std::numeric_limits<TTT>::max());

        std::pair<FFF, FFF> safe_lim = cvec.limits_estimate();

        lower = cvec.compress(safe_lim.first);
        upper = cvec.compress(safe_lim.second);
        
        check_extremes("no missing, safe",
                       lower, upper, middle,
                       (TTT)std::numeric_limits<TTT>::min(),
                       (TTT)std::numeric_limits<TTT>::max());
      }
      
      {
        cvec.set_scaling(lim, true);
        
        TTT lower = cvec.compress(lim.first);
        TTT upper = cvec.compress(lim.second);
        TTT middle = cvec.compress((lim.first + lim.second) / 2);
        
        if (std::numeric_limits<TTT>::is_signed) {
          check_not_missing("missing allowed",
                            lower, upper, middle,
                            std::numeric_limits<TTT>::min());

          check_extremes("missing allowed",
                         lower, upper, middle,
                         (TTT)(std::numeric_limits<TTT>::min() + 1),
                         (TTT)(std::numeric_limits<TTT>::max()));
          
          CPPUNIT_ASSERT_EQUAL_MESSAGE("no missing, middle",
                                       (TTT)(0),
                                       middle);
        } else {
          check_not_missing("missing allowed",
                            lower, upper, middle,
                            std::numeric_limits<TTT>::max());
          
          check_extremes("missing allowed",
                         lower, upper, middle,
                         (TTT)(0),
                         (TTT)(std::numeric_limits<TTT>::max() - 1));
        }
      }
    }
  }
  
  template <class TTT>
  void testInversion()
  {
    //    Eigen::Matrix<FFF, Eigen::Dynamic, 1> vec;
    //    std::pair<FFF, FFF> lim{-100.0, 100.0};
    std::vector<std::pair<FFF, FFF> > lims;
    lims.push_back(std::pair<FFF, FFF>(-10.0, 10.0));
    lims.push_back(std::pair<FFF, FFF>(-50.0, 50.0));
    lims.push_back(std::pair<FFF, FFF>(-50.0+273.15, 50.0+273.15));
    lims.push_back(std::pair<FFF, FFF>(-100.0, 100.0));
    lims.push_back(std::pair<FFF, FFF>(-100.0+273.15, 100.0+273.15));
    lims.push_back(std::pair<FFF, FFF>(-10.0, 10.0));
    lims.push_back(std::pair<FFF, FFF>(-50.0, 50.0));
    lims.push_back(std::pair<FFF, FFF>(-50.0+273.15, 50.0+273.15));
    lims.push_back(std::pair<FFF, FFF>(-100.0, 100.0));
    lims.push_back(std::pair<FFF, FFF>(-100.0+273.15, 100.0+273.15));
    
    for (typename std::vector<std::pair<FFF, FFF> >::iterator it =  lims.begin();
         it != lims.end();
         ++it) {
      std::pair<FFF, FFF> lim = *it;
      CompressedVector<FFF, TTT> cvec;
      
      for (int allow_missing = 0; allow_missing < 2; ++allow_missing) {
        cvec.set_scaling(lim, allow_missing > 0);

        FFF val = lim.first / FFF(3.0) + lim.second * FFF(2.0) / FFF(3.0);
        TTT comp = cvec.compress(val);
        FFF decomp = cvec.decompress(comp);
        TTT recomp = cvec.compress(decomp);
        
        CPPUNIT_ASSERT_EQUAL_MESSAGE("missing is " + to_string(allow_missing > 0) +
                                     ", recompression comparison",
                                     comp,
                                     recomp);
      }

      for (int allow_missing = 0; allow_missing < 2; ++allow_missing) {
        cvec.set_scaling(lim, allow_missing > 0);
        typename CompressedVector<FFF,TTT>::vector_type vec;
        vec.resize(1,1);

        vec[0] = lim.first / FFF(3.0) + lim.second * FFF(2.0) / FFF(3.0);
        cvec.compress(vec);
        typename CompressedVector<FFF,TTT>::vector_type decomp;
        cvec.decompress(decomp);
        CompressedVector<FFF, TTT> recomp;
        recomp.set_scaling(lim, allow_missing > 0);
        recomp.compress(decomp);
        
        CPPUNIT_ASSERT_EQUAL_MESSAGE("missing is " + to_string(allow_missing > 0) +
                                     ", recompression comparison",
                                     cvec[0],
                                     recomp[0]);
      }
}
  }


};




CPPUNIT_TEST_SUITE_REGISTRATION( TestNCDF_EigenScaling<float> );
CPPUNIT_TEST_SUITE_REGISTRATION( TestNCDF_EigenScaling<double> );
CPPUNIT_TEST_SUITE_REGISTRATION( TestNCDF_EigenScaling<long double> );
