#ifndef VERSIONS_HPP
#define VERSIONS_HPP

#include <iostream>
#include <string>
#include <map>

// Plain log definition, avoiding include debuglog.hpp
#ifndef PLOG_
#define PLOG_(msg) std::cout << msg;
#endif

/** Keeping track of the revision numbers of source files */
class VersionRegistry {
 private:
  typedef std::map<std::string, std::string> VersionMap;
  static VersionRegistry* registry_;
  VersionMap map_;
  VersionRegistry() {};
  VersionRegistry(VersionRegistry const&) {};
  VersionRegistry& operator=(VersionRegistry const&) { return *this; };
  
  int64_t extract_svn_version(const std::string& revision) const;
  void* internal_add(const std::string& filename,
                     const std::string& revision);
  std::string internal_get(const std::string filename);
  int64_t internal_get_value(const std::string filename);
  int64_t internal_get_latest();
  void internal_log();
  void internal_clear();

 public:
  static void purge() {
    //    PLOG_("Purging version registry" << std::endl);
    clear();
    if (registry_) {
      delete registry_;
      registry_ = NULL;
    }
  }
  static void clear() {
    //    PLOG_("Clearing version registry" << std::endl);
    if (registry_) {
      registry_->internal_clear();
    }
  }

  static VersionRegistry* registry() {
    if (!registry_) {
      //      PLOG_("Constructing version registry" << std::endl);
      registry_ = new VersionRegistry;
    }
    return registry_;
  }
  static void* add(const std::string& filename,
                   const std::string& revision) {
    return registry()->internal_add(filename, revision);
  }
  static std::string get(const std::string filename) {
    return registry()->internal_get(filename);
  }
  static int64_t get_value(const std::string filename) {
    return registry()->internal_get_value(filename);
  }
  static int64_t get_latest() {
    return registry()->internal_get_latest();
  }
  static void log() {
    registry()->internal_log();
  }
};

#define VERSIONS_ADD(prefix,file) namespace { \
  void* VERSIONS_MAKE_NAME(prefix) = VersionRegistry::add(__FILE__, file); \
  }
#define VERSIONS_GET() VersionRegistry::get(__FILE__)
#define VERSIONS_CLEAR() VersionRegistry::clear()
#define VERSIONS_PURGE() VersionRegistry::purge()

#define VERSIONS_MAKE_NAME(prefix) VERSIONS_JOIN(versions_dummy_, prefix)
/// Trick from cppuint:
#define VERSIONS_JOIN( symbol1, symbol2 ) _VERSIONS_DO_JOIN( symbol1, symbol2 )
/// \internal
#define _VERSIONS_DO_JOIN( symbol1, symbol2 ) _VERSIONS_DO_JOIN2( symbol1, symbol2 )
/// \internal
#define _VERSIONS_DO_JOIN2( symbol1, symbol2 ) symbol1##symbol2

VERSIONS_ADD(versions_hpp, "$Revision: 1283 $")

#endif
