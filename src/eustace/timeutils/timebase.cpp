#include "timebase.h"
#include <stdio.h>
#include <memory.h>
#include <iostream>

namespace EUSTACE
{
  CalendarDay::CalendarDay(int64_t year, int64_t month, int64_t day)
  {
    memset(reinterpret_cast<std::tm*>(this), 0, sizeof(std::tm));
    tm_year = year - 1900;
    tm_mon = month - 1;
    tm_mday = day; 
    mktime(this);
  }

  CalendarDay::CalendarDay(const std::string& str)
  {
    int year(0), month(0), day(0);
    sscanf(str.c_str(), "%04d%02d%02d", &year, &month, &day);
    *this = CalendarDay(year, month, day);
  }

  void CalendarDay::Text(std::string& str) const
  {
    char buffer[16];
    memset(buffer, 0, sizeof(buffer));
    Text(buffer);
    str = buffer;
  }

  void CalendarDay::Text(char* buffer) const
  {
    std::strftime(buffer, 9, "%Y%m%d", this);
  }

  int64_t TimeBaseDays::Number(const std::tm& time) const
  {
    std::tm t0 = epoch;
    std::tm t1 = time;

    return std::difftime(std::mktime(&t1), std::mktime(&t0)) / 86400;
  }

  std::tm TimeBaseDays::Time(int64_t daynumber) const
  {
    std::tm epoch_offset = epoch;
    epoch_offset.tm_mday += daynumber;
    mktime(&epoch_offset);
    return epoch_offset;
  }
}
