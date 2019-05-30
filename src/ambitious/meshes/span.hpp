#ifndef SPAN_HPP
#define SPAN_HPP

#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(span_hpp, "$Revision: 1290 $")

#include "ambitious/debuglog/debuglog.hpp"

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>

#include <eustace/timeutils/timebase.h>
#include "ambitious/timer/timer.hpp"
#include "ambitious/to_string/to_string.hpp"




template <typename Value>
class Span {
 public:
  typedef Value value_type;
  typedef std::pair<value_type, value_type> SpanBase;
 protected:
  SpanBase span_;
  bool infinite_;
  bool empty_;
 public:
  Span(value_type a, value_type b)
      : span_(a, b),
        infinite_(false),
        empty_(a > b)
  {
  }
  Span(bool is_infinite)
      : span_(),
        infinite_(is_infinite),
        empty_(!is_infinite)
  {
    if (is_infinite) {
      set_infinite();
    } else {
      set_empty();
    }
  }

  value_type first() const {
    return span_.first;
  }
  value_type second() const {
    return span_.second;
  }
  bool empty() const {
    return empty_;
  }
  bool infinite() const {
    return infinite_;
  }

  const SpanBase& span() const {
    return span_;
  }
  bool before(value_type day) const {
    return (!empty_ &&
            !infinite_ &&
            (day < span_.first));
  }
  bool after(value_type day) const {
    return (!empty_ &&
            !infinite_ &&
            (day > span_.second));
  }
  bool inside(value_type day) const {
    return (!empty_ &&
            (infinite_ ||
             ((day >= span_.first) &&
              (day <= span_.second))));
  }
  bool outside(value_type day) const {
    return (empty_ ||
            (!infinite_ &&
             ((day < span_.first) ||
              (day > span_.second))));
  }

  bool intersects(const Span& other) const {
    if (empty_ || other.empty_) {
      return false;
    }
    return ((infinite_ || other.infinite_) ||
            ((span_.first <= other.span_.second) &&
             (span_.second >= other.span_.first)));
  }

  bool covers(const Span& other) const {
    if (empty_ || other.empty_) {
      return false;
    }
    return (infinite_ ||
            (!other.infinite_ &&
             (span_.first <= other.span_.first) &&
             (span_.second >= other.span_.second)));
  }

  Span& copy(const Span& other) {
    span_ = other.span_;
    infinite_ = other.infinite_;
    empty_ = other.empty_;
    return *this;
  }

  /** Join two spans
   *
   * Non-overlapping spans are joined by including the gap between them.
   */
  Span& join(const Span& other) {
    if (!infinite_ && !other.empty_) {
      if (empty_) {
        copy(other);
      } else if (other.infinite_) {
        set_infinite();
      } else {
        if (span_.first > other.span_.first) {
          span_.first = other.span_.first;
        }
        if (span_.second < other.span_.second) {
          span_.second = other.span_.second;
        }
      }
    }
    return *this;
  }
  /** Compute the intersection of two spans
   */
  Span& intersect(const Span& other) {
    if (!empty_ && !other.infinite_) {
      if (other.empty_) {
        this->set_empty();
      } else {
        if (infinite_) {
          copy(other);
        } else {
          if (span_.first < other.span_.first) {
            span_.first = other.span_.first;
          }
          if (span_.second > other.span_.second) {
            span_.second = other.span_.second;
          }
          empty_ = span_.first > span_.second;
        }
      }
    }
    return *this;
  }
  
  Span& set(value_type a, value_type b) {
    span_.first = a;
    span_.second = b;
    infinite_ = false;
    empty_ = a > b;
    return *this;
  }
  Span& first(value_type a) {
    span_.first = a;
    infinite_ = false;
    empty_ = a > span_.second;
    return *this;
  }
  Span& second(value_type b) {
    span_.second = b;
    infinite_ = false;
    empty_ = span_.first > b;
    return *this;
  }
  Span& set_infinite() {
    span_.first = std::numeric_limits<value_type>::min();
    span_.second = std::numeric_limits<value_type>::min();
    infinite_ = true;
    empty_ = false;
    return *this;
  }
  Span& set_empty() {
    span_.first = std::numeric_limits<value_type>::max();
    span_.second = std::numeric_limits<value_type>::min();
    infinite_ = false;
    empty_ = true;
    return *this;
  }
};

#ifdef OLDCXX
template<class Value>
std::ostream& operator<<(std::ostream& output,
                         Span<Value>& S) {
#else
template<class Value>
std::ostream& operator<<(std::ostream& output,
                         const Span<Value>& S) {
#endif
  if (S.empty()) {
    output << "Empty";
  } else if (S.infinite()) {
    output << "(-Inf, Inf)";
  } else {
    output << "(" << S.first() << ", " << S.second() << ")";
  }
  return output;
}



#endif
