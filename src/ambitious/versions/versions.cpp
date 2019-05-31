#include "ambitious/versions/versions.hpp"

VersionRegistry* VersionRegistry::registry_ = NULL;

VERSIONS_ADD(versions_cpp, "$Revision: 1278 $")


int64_t VersionRegistry::extract_svn_version(const std::string& revision) const {
  size_t first_blank = revision.find_first_of(" ", 0);
  if (first_blank == std::string::npos) {
    return 0;
  }
  size_t last_blank = revision.find_last_of(" ", first_blank + 1);
  if (last_blank == std::string::npos) {
    return 0;
  }
  std::string version_str = revision.substr(first_blank + 1, last_blank - first_blank - 1);
  return stoll(version_str);
}

void* VersionRegistry::internal_add(const std::string& filename,
                                    const std::string& revision) {
  map_.insert(VersionMap::value_type(filename, revision));
  return NULL;
}

std::string VersionRegistry::internal_get(const std::string filename) {
  VersionMap::iterator info;
  info = map_.find(filename);
  if (info != map_.end()) {
    return info->second;
  } else {
    return std::string("Unknown revision");
  }
}

int64_t VersionRegistry::internal_get_value(const std::string filename) {
  return extract_svn_version(internal_get(filename));
}

int64_t VersionRegistry::internal_get_latest() {
  VersionMap::iterator info;
  int64_t latest_version = 0;
  for (info = map_.begin();
       info != map_.end();
       ++info) {
    int64_t version = extract_svn_version(info->second);
    if (version > latest_version) {
      latest_version = version;
    }
  }
  return latest_version;
}

void VersionRegistry::internal_log() {
  std::cout << "Latest revision: " << internal_get_latest() << std::endl;
  for (VersionMap::iterator info = map_.begin();
       info != map_.end();
       ++info) {
    std::cout << info->first << " = " << info->second << std::endl;
  }
}

void VersionRegistry::internal_clear() {
  map_.clear();
}
