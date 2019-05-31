#ifndef TO_STRING_HPP
#define TO_STRING_HPP

#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(to_string_hpp, "$Revision: 1297 $")

/** \file to_string.hpp
 *
 * Workaround for ambiguous overload of std::to_string in older compilers.
 * [Stackoverflow reference](https://stackoverflow.com/questions/10664699/stdto-string-more-than-instance-of-overloaded-function-matches-the-argument)
 */

#include <type_traits>
#include <string>
#include <vector>
#include <set>

template<typename T>
inline
typename std::enable_if<std::is_integral<T>::value && std::is_signed<T>::value, std::string>::type
to_string(T const val) {
    return std::to_string(static_cast<long long>(val));
}

template<typename T>
inline
typename std::enable_if<std::is_integral<T>::value && std::is_unsigned<T>::value, std::string>::type
to_string(T const val) {
    return std::to_string(static_cast<unsigned long long>(val));
}

template<typename T>
inline
typename std::enable_if<std::is_floating_point<T>::value, std::string>::type
to_string(T const val) {
    return std::to_string(static_cast<long double>(val));
}


#endif
