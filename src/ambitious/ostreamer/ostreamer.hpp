#ifndef OSTREAMER_HPP
#define OSTREAMER_HPP

#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(ostreamer_hpp, "$Revision: 1297 $")

#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <iterator>

/** A helper template function for printing vectors to streams
 */
template<class T>
std::ostream& ostream_copy_vec(std::ostream& os,
			       const std::vector<T>& vec,
			       const char * sep = " ")
{
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(os, sep));
  return os;
}

/** A helper template function for printing sets to streams
 */
template<class T>
std::ostream& ostream_copy_set(std::ostream& os,
			       const std::set<T>& vec,
			       const char * sep = " ")
{
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(os, sep));
  return os;
}

/** Template operator for printing vectors to streams
 */
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
  return ostream_copy_vec(os, vec, " ");
}

/** Template operator for printing sets to streams
 */
template<class T>
std::ostream& operator<<(std::ostream& os, const std::set<T>& vec)
{
  return ostream_copy_set(os, vec, " ");
}

#endif
