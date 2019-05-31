#include <cppunit/extensions/HelperMacros.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <ambitious/bidirmap/bidirmap.hpp>


class TestBidirmap  : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( TestBidirmap );
  CPPUNIT_TEST( testBidirmapInsert );
  CPPUNIT_TEST_SUITE_END();

public:
  
  void testBidirmapInsert()
  {
    BidirectionalMap<int64_t,int64_t> bidirmap;

    typedef std::pair<int64_t, int64_t > Pair;

    std::vector<Pair> data;
    data.push_back(Pair(1, 2));
    data.push_back(Pair(3, 0));
    data.push_back(Pair(2, 1));

    for (unsigned int i=0; i < data.size(); ++i) {
      bidirmap.insert(data[i].first, data[i].second);
      CPPUNIT_ASSERT_EQUAL_MESSAGE("data[i].second vs bidirmap.AtoB(data[i].first)",
				   data[i].second,
				   bidirmap.AtoB(data[i].first));
      CPPUNIT_ASSERT_EQUAL_MESSAGE("data[i].first vs bidirmap.AtoB(data[i].second)",
				   data[i].first,
				   bidirmap.BtoA(data[i].second));
    }

    CPPUNIT_ASSERT_EQUAL_MESSAGE("data.size vs bidirmap.AtoB.size",
				 data.size(),
				 bidirmap.AtoB().size());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("data.size vs bidirmap.BtoA.size",
				 data.size(),
				 bidirmap.BtoA().size());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("data.size vs bidirmap.size",
				 data.size(),
				 bidirmap.size());

    for (unsigned int i=0; i < data.size(); ++i) {
      CPPUNIT_ASSERT_EQUAL_MESSAGE("data[i].second vs bidirmap.AtoB(data[i].first)",
				   data[i].second,
				   bidirmap.AtoB(data[i].first));
      CPPUNIT_ASSERT_EQUAL_MESSAGE("data[i].first vs bidirmap.AtoB(data[i].second)",
				   data[i].first,
				   bidirmap.BtoA(data[i].second));
    }
    
  }


};

CPPUNIT_TEST_SUITE_REGISTRATION( TestBidirmap );
