#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(debuglog_cpp, "$Revision: 1292 $")

#include "ambitious/debuglog/debuglog.hpp"

namespace debuglog {

const std::string SHORTEN_PREFIX = "ambitious/";

std::string shorten_filename(const std::string& filename) {
  size_t pos = filename.find(SHORTEN_PREFIX, 0);
  if (pos != std::string::npos) {
    return filename.substr(pos, filename.length());
  } else {
    return filename;
  }
}

}
