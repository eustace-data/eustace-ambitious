#ifndef TIMER_HPP
#define TIMER_HPP

#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(timer_hpp, "$Revision: 1283 $")

#include <string>
#include <iostream>
#include <chrono>
#include <list>
#include <stack>


/** Helper macro for shortcuts to a private timer pointer in a class */
#define TIMER_HELPERS \
  void timer_tic(const std::string& msg) { \
    if (timer_) timer_->tic(msg); \
  } \
  void timer_toc() { \
    if (timer_) timer_->toc(); \
  } \
  void timer_toc(int i) { \
    if (timer_) timer_->toc(i); \
  } \
  void timer_toc_tic(const std::string& msg) { \
    if (timer_) timer_->toc().tic(msg); \
  } \
  void timer_toc_tic(int i, const std::string& msg) { \
    if (timer_) timer_->toc(i).tic(msg); \
  }


class Timer {
  std::chrono::high_resolution_clock::time_point t0_;
  std::chrono::high_resolution_clock::time_point t1_;
public:
  Timer();
  Timer& tic();
  Timer& toc();
  std::chrono::duration<double> time() const;
  double count() const;
  std::chrono::high_resolution_clock::time_point start();
  std::chrono::high_resolution_clock::time_point ending();

  Timer& print(const std::string& msg, const std::string& msg2 = std::string(""));
  Timer& print(const std::string& msg, double scaling, const std::string& msg2 = std::string(""));
  Timer& printinverted(const std::string& msg, double scaling, const std::string& msg2 = std::string(""));
};

class TimerHierarchy;
class TimerHierarchy {

  class TimerRecord;

  class TimerList : public std::list<TimerRecord> {
  public:
    TimerList& print(const std::string& prefix);
  };
  
  class TimerRecord {
  public:
    TimerHierarchy& head_;
    std::string name_;
    Timer timer_;
    TimerList subtimers_;
    TimerRecord(TimerHierarchy& head, const std::string& name) :
      head_(head), name_(name), timer_(), subtimers_() {};
    TimerRecord& print(const std::string& prefix);
  };

  TimerList timers_;
  std::stack<TimerRecord*> active_;
  int verbosity_;
  std::string subprefix_;

public:
  TimerHierarchy(int verbosity=0, const std::string& subprefix="    ") :
    timers_(), active_(), verbosity_(verbosity), subprefix_(subprefix) {};
  // Return the current level; size of active_ minus 1. -1 means no active timers.
  int level() const;
  // Start a sub-timer
  TimerHierarchy& tic(const std::string& name, int verbosity);
  TimerHierarchy& tic(const std::string& name = std::string(""));
  // End the active timer
  TimerHierarchy& toc();
  // Return a ref to the active timer
  Timer& timer();
  // End the active timer and return a ref to it
  Timer& toc_timer();
  // End all sub-timers above level. -1 ends all timers.
  TimerHierarchy& toc(int level);
  TimerHierarchy& print(const std::string& prefix = std::string(""));
};

#endif
