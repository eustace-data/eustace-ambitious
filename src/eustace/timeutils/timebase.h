#ifndef _EUSTACE_TIMEUTILS_TIMEBASE_H
#define _EUSTACE_TIMEUTILS_TIMEBASE_H

#include <ctime>
#include <string>
#include <stdint.h>

namespace EUSTACE
{
  class CalendarDay : public std::tm
  {
  public:
    CalendarDay(const std::tm& src) :
      std::tm(src)
    {
    }
    
    CalendarDay(int64_t year, int64_t month, int64_t day);

    CalendarDay(const std::string& str);

  public:

    void Text(std::string& str) const;

  private:
    void Text(char* buffer) const;
  };

  class TimeBaseDays
  {
  public:
    TimeBaseDays(const std::tm& _epoch) :
      epoch(_epoch)
    {
    }

  public:

    int64_t Number(const std::tm& time) const;

    std::tm Time(int64_t daynumber) const;

  private:
    std::tm epoch;
  };
}

#endif
