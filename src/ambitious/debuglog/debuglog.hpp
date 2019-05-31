#ifndef DEBUGLOG_HPP
#define DEBUGLOG_HPP

#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(debuglog_hpp, "$Revision: 1290 $")

#include <iostream>
#include <string>

// Define NDEBUG to disable assert
#include <cassert>

namespace debuglog {
std::string shorten_filename(const std::string& filename);
}

#ifndef WHEREAMI
#define WHEREAMI debuglog::shorten_filename(__FILE__) << "(" << __LINE__ << ")\t"
#endif

// Plain log
#ifndef PLOG_
#define PLOG_(msg) std::cout << msg;
#endif

// Verbose log
#ifndef LOG_
#define LOG_(msg) std::cout << WHEREAMI << msg;
#endif

// Error log
#ifndef ELOG_
#define ELOG_(msg) std::cerr << WHEREAMI << msg;
#endif

// Conditional debug log
#ifndef LOG
#ifdef DEBUG
#define LOG(msg) LOG_(msg)
#else
#define LOG(msg)
#endif
#endif

#endif
