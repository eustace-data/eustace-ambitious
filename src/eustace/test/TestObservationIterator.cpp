#include <cppunit/extensions/HelperMacros.h>
#include <eustace/analysis/observationsource.h>
#include <iostream>
#include <fstream>

using namespace EUSTACE;
using namespace std;

struct TEST_OBS_ITER_DATA
{
  uint64_t locationID;
  double measurement;
  double uncertainty_uncorrelated;
  double uncertainty_locally_correlated[2];
};

class TestObservationIterator  : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( TestObservationIterator );
  CPPUNIT_TEST( testWithExampleData );
  CPPUNIT_TEST_SUITE_END();

public:
  
  void testWithExampleData()
  {
    // make some example data
    struct TEST_OBS_ITER_DATA testdata[5] = {
      { 23, 278.22, 0.012, { 0.177, 2.999 } },
      { 18, 250.80, 0.333, { 0.222, 0.123 } },
      {328, 290.11, 0.012, { 1.102, 3.888 } },
      { 52, 310.80, 0.001, { 5.555,18.678 } },
      {  1, 320.00, 0.111, { 0.002, 2.000 } },
    };

    // copy into ObservationEntries structure
    ObservationEntries obs;
    obs.SetDimensions(2, 5);

    // check sizes as expected
    CPPUNIT_ASSERT_EQUAL(size_t(5), obs.NumEntries());
    CPPUNIT_ASSERT_EQUAL(sizeof(testdata), obs.TotalSizeBytes());

    // copy data into obs structure
    memcpy(obs.Bytes(), testdata, sizeof(testdata));

    // Make iterator
    ObservationIterator iter(obs);

    // Should have first element
    CPPUNIT_ASSERT_EQUAL(true, iter.HasElement());

    // Iteration 0
    CPPUNIT_ASSERT_EQUAL(uint64_t(23), iter.Identifier());
    CPPUNIT_ASSERT_EQUAL(278.22, iter.Measurement());
    CPPUNIT_ASSERT_EQUAL(0.012, iter.UncertaintyUncorrelated());
    CPPUNIT_ASSERT_EQUAL(0.177, iter.UncertaintyLocallyCorrelated(0));
    CPPUNIT_ASSERT_EQUAL(2.999, iter.UncertaintyLocallyCorrelated(1));
    CPPUNIT_ASSERT_EQUAL(true, iter.HasElement());

    // Iteration 1
    iter.Next();
    CPPUNIT_ASSERT_EQUAL(uint64_t(18), iter.Identifier());
    CPPUNIT_ASSERT_EQUAL(250.80, iter.Measurement());
    CPPUNIT_ASSERT_EQUAL(0.333, iter.UncertaintyUncorrelated());
    CPPUNIT_ASSERT_EQUAL(0.222, iter.UncertaintyLocallyCorrelated(0));
    CPPUNIT_ASSERT_EQUAL(0.123, iter.UncertaintyLocallyCorrelated(1));
    CPPUNIT_ASSERT_EQUAL(true, iter.HasElement());

    // Iteration 2
    iter.Next();
    CPPUNIT_ASSERT_EQUAL(uint64_t(328), iter.Identifier());
    CPPUNIT_ASSERT_EQUAL(290.11, iter.Measurement());
    CPPUNIT_ASSERT_EQUAL(0.012, iter.UncertaintyUncorrelated());
    CPPUNIT_ASSERT_EQUAL(1.102, iter.UncertaintyLocallyCorrelated(0));
    CPPUNIT_ASSERT_EQUAL(3.888, iter.UncertaintyLocallyCorrelated(1));
    CPPUNIT_ASSERT_EQUAL(true, iter.HasElement());

    // Iteration 3
    iter.Next();
    CPPUNIT_ASSERT_EQUAL(uint64_t(52), iter.Identifier());
    CPPUNIT_ASSERT_EQUAL(310.80, iter.Measurement());
    CPPUNIT_ASSERT_EQUAL(0.001, iter.UncertaintyUncorrelated());
    CPPUNIT_ASSERT_EQUAL(5.555, iter.UncertaintyLocallyCorrelated(0));
    CPPUNIT_ASSERT_EQUAL(18.678, iter.UncertaintyLocallyCorrelated(1));
    CPPUNIT_ASSERT_EQUAL(true, iter.HasElement());

    // Iteration 4
    iter.Next();
    CPPUNIT_ASSERT_EQUAL(uint64_t(1), iter.Identifier());
    CPPUNIT_ASSERT_EQUAL(320.00, iter.Measurement());
    CPPUNIT_ASSERT_EQUAL(0.111, iter.UncertaintyUncorrelated());
    CPPUNIT_ASSERT_EQUAL(0.002, iter.UncertaintyLocallyCorrelated(0));
    CPPUNIT_ASSERT_EQUAL(2.000, iter.UncertaintyLocallyCorrelated(1));
    CPPUNIT_ASSERT_EQUAL(true, iter.HasElement());

    // Should say that there aren't any more iterations available
    iter.Next();
    CPPUNIT_ASSERT_EQUAL(false, iter.HasElement());
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestObservationIterator );
