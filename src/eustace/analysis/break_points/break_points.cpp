#include "break_points.h"


namespace EUSTACE
{
  std::vector<int32_t> ObservationBreakingPoints::stations_collection() const
  {   
    
    std::vector<int32_t> vector_copy;
    vector_copy.reserve(break_station.size());
    copy(break_station.begin(),break_station.end(),back_inserter(vector_copy));

    std::sort(vector_copy.begin(), vector_copy.end());
    vector_copy.erase(std::unique( vector_copy.begin(), vector_copy.end() ), vector_copy.end() );
    
    return vector_copy;
  }

  int  ObservationBreakingPoints::number_of_stations() const
  {

    std::vector<int32_t> stations = stations_collection();
    return stations.size(); 

  }
  
  std::vector<int32_t> ObservationBreakingPoints::stations_count() const 
  {
    
    std::vector<int32_t> stations = stations_collection();
    std::vector<int32_t> count;
    
    for (std::vector<int32_t>::iterator it = stations.begin(); it != stations.end(); ++it)
      count.push_back(std::count(break_station.begin(), break_station.end(), *it));
    
    return count;
    
  } 
  
  std::vector<int32_t> ObservationBreakingPoints::filtered_breakpoints(int32_t station_index) const
  {
    
    std::vector<int32_t> filtered_break_points;

    for (std::vector<int32_t>::const_iterator it = break_station.begin(); it != break_station.end(); ++it){
      if(*it == station_index)
	filtered_break_points.push_back(break_time.at(std::distance(break_station.begin(), it)));

    }

    return filtered_break_points;
  }
  
  void ObservationBreakingPoints::apply_policy(int policy_code, double threshold, double decay_constant)
  {
    std::runtime_error Error("Wrong policy code: policy code can be either HARD_CUT_OFF (0) or LAPLACE_KERNEL (1)!");
    std::cout << POLICIES[0] << std::endl;
    if (policy_code == 0)
      apply_base_policy(threshold);
    else if (policy_code == 1)
      apply_kernel_policy(threshold, decay_constant);
    else
      throw Error;
  }
  
  void ObservationBreakingPoints::apply_base_policy(double threshold)
  {    
    
    std::vector<int> discriminant_vector =break_likelihood;
    
    std::vector<int32_t> vector_copy;
    reject<int32_t,int>(&break_time, vector_copy, discriminant_vector, threshold);
    reject<int32_t,int>(&break_station, vector_copy, discriminant_vector, threshold);
    reject<int32_t,int>(&break_likelihood, vector_copy, discriminant_vector, threshold);

  }
  
  void ObservationBreakingPoints::apply_kernel_policy(double threshold, double decay_constant)
  {

    std::vector<double> discriminant_vector;
    for (std::vector<int>::iterator it = break_likelihood.begin(); it != break_likelihood.end(); ++it)
      discriminant_vector.push_back(exp(-abs(*it-MAX_LIKELIHOOD)*decay_constant));
  
  
  
    std::vector<int32_t> vector_copy;
    reject<int32_t,double>(&break_time, vector_copy, discriminant_vector, threshold);
    reject<int32_t,double>(&break_station, vector_copy, discriminant_vector, threshold);
    reject<int32_t,double>(&break_likelihood, vector_copy, discriminant_vector, threshold);

  }
  
  template <typename T1, typename T2>
  void ObservationBreakingPoints::reject(std::vector<T1>* member_vector, std::vector<T1> member_vector_copy, std::vector<T2> discriminant_vector, double threshold)
  {
    
    for (unsigned i = 0; i < discriminant_vector.size(); i++)
      {
	
	if(discriminant_vector[i] >= threshold){
	  member_vector_copy.push_back(member_vector->at(i));
	}
      }
    
    member_vector->clear();
    member_vector->resize(member_vector_copy.size());
    std::copy(member_vector_copy.begin(), member_vector_copy.end(), member_vector->begin());
    member_vector_copy.clear();
    
  }
}
