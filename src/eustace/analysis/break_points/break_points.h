#ifndef _EUSTACE_ANALYSIS_BREAK_POINTS_BREAK_POINTS_H_
#define _EUSTACE_ANALYSIS_BREAK_POINTS_BREAK_POINTS_H_

#include <algorithm>
#include <iostream>
#include <math.h>
#include <memory.h>
#include <stdexcept>
#include <stdint.h>
#include <string>
#include <vector>

namespace EUSTACE
{
  class ObservationBreakingPoints
  {
     
  public:
    std::vector<int32_t> break_time;
    std::vector<int32_t> break_station;
    std::vector<int> break_likelihood;
    
  public:
  ObservationBreakingPoints(std::vector<int32_t> _break_time, 
			    std::vector<int32_t> _break_station,
			    std::vector<int> _break_likelihood) :
    break_time(_break_time),
      break_station(_break_station),
      break_likelihood(_break_likelihood)
	{
	  POLICIES[0]="HARD_CUT_OFF";
	  POLICIES[1]="LAPLACE_KERNEL";
	}

  public:
    ~ObservationBreakingPoints()
      {}
  public:
    int number_of_observations()
    { return break_time.size(); }
    
    //Stations basic statistics
    std::vector<int32_t> stations_collection() const;
    int number_of_stations() const;
    std::vector<int32_t> stations_count() const;
    
    //Break points filtering
    std::vector<int32_t> filtered_breakpoints(int32_t station_index) const;

    //Break points rejection
    void apply_policy(int policy_code, double threshold = 0., double decay_constant = 0.);
    void apply_base_policy(double threshold);
    void apply_kernel_policy(double threshold, double decay_constant);
    
  private:
    static const int MAX_LIKELIHOOD = 127;
    const char *POLICIES[2];
    
    //Helper for rejecting breaking points
    template <typename T1, typename T2> 
      void reject(std::vector<T1>* member_vector, std::vector<T1> member_vector_copy, std::vector<T2> discriminant_vector, double threshold);
  };
}

#endif
