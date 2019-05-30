#include "ambitious/debuglog/debuglog.hpp"
#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(cal_cpp, "$Revision: 1283 $")

#include <iostream>
#include "ambitious/timer/timer.hpp"
#include "ambitious/to_string/to_string.hpp"
#include "ambitious/ostreamer/ostreamer.hpp"
#include <string.h>
#include "ambitious/po/po.hpp"
#include "ambitious/observe/observe.hpp"
#include <eustace/analysis/inputmanager.h>
#include "ambitious/meshes/meshes.hpp"

using namespace EUSTACE;


int main(int argc, char* argv[])
{
  if (argc == 2) {
    std::tm tm = DateTimeHelper(1850).get_tm(atoi(argv[1]));
    char prev = std::cout.fill('0');
    std::cout <<
        std::setw(4) << tm.tm_year + 1900 << "-" <<
        std::setw(2) << tm.tm_mon + 1 << "-" <<
        std::setw(2) << tm.tm_mday <<
        std::endl;
    std::cout.fill(prev);
  } else if (argc == 4) {
    int64_t day = DateTimeHelper(1850).get_day(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
    std::cout <<
        day <<
        std::endl;
  } else {
    std::cout << "Usage: cal day, or cal year month day" << std::endl;
  }

  VERSIONS_PURGE();
  return 0;
}


