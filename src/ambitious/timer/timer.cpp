#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(timer_cpp, "$Revision: 1283 $")

#include "ambitious/timer/timer.hpp"
#include <cassert>

Timer::Timer() {
  tic();
}
Timer& Timer::tic() {
  t0_ = std::chrono::high_resolution_clock::now();
  t1_ = t0_;
  return *this;
}
Timer& Timer::toc() {
  t1_ = std::chrono::high_resolution_clock::now();
  return *this;
}
std::chrono::duration<double> Timer::time() const {
  return std::chrono::duration_cast<std::chrono::duration<double>>(t1_ - t0_);
}
double Timer::count() const {
  return time().count();
}
std::chrono::high_resolution_clock::time_point Timer::start() {
  return t0_;
}
std::chrono::high_resolution_clock::time_point Timer::ending() {
  return t1_;
}

Timer& Timer::print(const std::string& msg, const std::string& msg2) {
  std::cout << msg << time().count() << msg2 << std::endl;
  return *this;
}

Timer& Timer::print(const std::string& msg, double scaling, const std::string& msg2) {
  std::cout << msg << time().count() * scaling << msg2 << std::endl;
  return *this;
}
Timer& Timer::printinverted(const std::string& msg, double scaling, const std::string& msg2) {
  std::cout << msg << scaling / time().count() << msg2 << std::endl;
  return *this;
}






// Return the current level; size of active_, minus 1. -1 means no active timers.
int TimerHierarchy::level() const {
  return int(active_.size()) - 1;
}
// Start a sub-timer
TimerHierarchy& TimerHierarchy::tic(const std::string& name, int verbosity) {
  if (verbosity > 0) {
    std::cout << name << std::endl;
  }
  if (active_.empty()) {
    // No active timer, add global timer
    timers_.push_back(TimerRecord(*this, name));
    active_.push(&timers_.back());
  } else {
    // Active timer, add subtimer
    (*active_.top()).subtimers_.push_back(TimerRecord(*this, name));
    active_.push(&(*active_.top()).subtimers_.back());
  }
  return *this;
}
TimerHierarchy& TimerHierarchy::tic(const std::string& name) {
  return tic(name, verbosity_);
}
// End the active timer
TimerHierarchy& TimerHierarchy::toc() {
  if (!active_.empty()) {
    (*active_.top()).timer_.toc();
    if (verbosity_ > 1) {
      // TODO: fix this to include an appropriate prefix
      (*active_.top()).timer_.print((*active_.top()).name_ + ":\t");
    }
    active_.pop();
  }
  return *this;
}
// End all sub-timers above level. -1 ends all timers.
TimerHierarchy& TimerHierarchy::toc(int level) {
  while (int(active_.size()) > level+1) {
    toc();
  }
  return *this;
}

// End the active timer and return a ref to it
Timer& TimerHierarchy::toc_timer() {
  assert(!active_.empty()); // "Attempt to access nonexistent timer!"
  Timer& timer = (*active_.top()).timer_;
  toc();
  return timer;
}

// Return a ref to the active timer
Timer& TimerHierarchy::timer() {
  assert(!active_.empty()); // "Attempt to access nonexistent timer!"
  return (*active_.top()).timer_;
}


TimerHierarchy::TimerList& TimerHierarchy::TimerList::print(const std::string& prefix) {
  for (auto iter = (*this).begin(); iter != (*this).end(); ++iter) {
    //    std::cout << prefix << *iter << std::endl;
    (*iter).print(prefix);
  }
  return *this;
}

TimerHierarchy::TimerRecord& TimerHierarchy::TimerRecord::print(const std::string& prefix) {
  timer_.print(prefix + name_ + ":\t");
  subtimers_.print(prefix + head_.subprefix_);
  return *this;
}

TimerHierarchy& TimerHierarchy::print(const std::string& prefix) {
  timers_.print(prefix);
  return *this;
}
