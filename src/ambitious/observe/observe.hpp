#ifndef OBSERVE_HPP
#define OBSERVE_HPP

#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(observe_hpp, "$Revision: 1290 $")

#include "ambitious/debuglog/debuglog.hpp"
#include <eustace/analysis/inputmanager.h>
#include <iostream>
#include <set>
#include <map>
#include <Eigen/Dense>
#include "ambitious/meshes/meshes.hpp"
#include "ambitious/timer/timer.hpp"
#include "ambitious/to_string/to_string.hpp"

using namespace EUSTACE;









class ModelObservationSource {
 public:
  typedef std::set<int64_t> MissingDays;

 protected:
  //  AnalysisInput::sourcetype sourcetype_; // source_[satellite|insitu]_[land|ocean|ice|lake], num_sources
  //  AnalysisInput::observabletype observabletype_; // observable_[tmin|tmax|tmean], num_observables

  AnalysisInputManager& inputmanager_;
  AnalysisInput input_;
  SpaceTimeMapper* mask_;
  TimerHierarchy* timer_;
  MissingDays missing_days_;
  bool missing_days_is_built_;
  int64_t day_;

  std::vector<uint64_t> locationID_;
  std::vector<LocationEntry> latlong_;
  Eigen::Matrix<double,Eigen::Dynamic,1> measurement_;
  Eigen::Matrix<double,Eigen::Dynamic,1> uncertainty_uncorrelated_;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> uncertainty_locally_correlated_;
  Eigen::Matrix<double,Eigen::Dynamic,1> prec_uncorrelated_;
 public:
  ModelObservationSource(AnalysisInputManager& inputmanager,
                         const AnalysisInput& input,
                         SpaceTimeMapper* mask,
                         TimerHierarchy* timer) :
      inputmanager_(inputmanager),
      input_(input),
      mask_(mask),
      timer_(timer),
      missing_days_(),
      missing_days_is_built_(false),
      day_(-1)
  {
    locationID_.resize(0);
    latlong_.resize(0);
    measurement_.resize(0,1);
    uncertainty_uncorrelated_.resize(0,1);
    uncertainty_locally_correlated_.resize(0,0);
    prec_uncorrelated_.resize(0,1);
  }

  std::vector<uint64_t>& locationID() {
    return locationID_;
  }
  std::vector<LocationEntry>& latlong() {
    return latlong_;
  }
  Eigen::Matrix<double,Eigen::Dynamic,1>& measurement() {
    return measurement_;
  }
  Eigen::Matrix<double,Eigen::Dynamic,1>& uncertainty_uncorrelated() {
    return uncertainty_uncorrelated_;
  }
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& uncertainty_locally_correlated() {
    return uncertainty_locally_correlated_;
  }
  Eigen::Matrix<double,Eigen::Dynamic,1>& prec_uncorrelated() {
    return prec_uncorrelated_;
  }



  
  const AnalysisInput& analysis_input() {
    return input_;
  }

  std::string source_name(const AnalysisInput& input) {
    return std::string(input.sourcename[input.Source()]);
  }
  std::string source_name() {
    return source_name(input_);
  }
  std::string observable_name(const AnalysisInput& input) {
    return std::string(input.observablename[input.Observable()]);
  }
  std::string observable_name() {
    return observable_name(input_);
  }

  bool is_mobile() const {
    return (input_.Source() == AnalysisInput::source_insitu_ocean);
  }
  bool is_cell() const {
    return ((input_.Source() != AnalysisInput::source_insitu_land) &&
            (input_.Source() != AnalysisInput::source_insitu_ocean));
  }

  void summary_log() {
    PLOG_("\tSource:\t" << source_name() <<
          "\tObservable:\t" << observable_name() <<
          "\tObservations:\t" << measurement_.rows() <<
          "\tCorrelated:\t" << uncertainty_locally_correlated_.cols() <<
          "\tRange:\t" << input_.StartDay() << "-" << input_.EndDay() <<
          "\tMissing:\t" << missing_days().size() <<
          std::endl);
  }

  bool obs_file_exists(int64_t day) {
    std::ifstream inputstream;
    inputstream.exceptions(std::ifstream::goodbit); // Don't throw exceptions
    std::string pathname_observation;
    inputmanager_.PathNameObservation(pathname_observation,
                                      input_.Source(),
                                      input_.Observable(),
                                      day);
    inputstream.open(pathname_observation.c_str(), std::ios::binary);
    if (inputstream.good()) {
      inputstream.close();
      return true;
    } else {
      inputstream.clear();
    }
    return false;
  }

  std::set<int64_t>& missing_days() {
    if (!missing_days_is_built_) {
      build_missing_days();
    }
    return missing_days_;
  }

  bool missing_day(int64_t day) {
    if (!missing_days_is_built_) {
      build_missing_days();
    }
    if ((input_.StartDay() < 0) ||
        (day < input_.StartDay()) ||
        (day > input_.EndDay())) {
      return true;
    }
    MissingDays::iterator it = missing_days_.find(day);
    return (it != missing_days_.end());
  }

  int64_t build_missing_days() {
    missing_days_.clear();
    missing_days_is_built_ = true;
    if (input_.StartDay() < 0 ) {
      return 0;
    }
    for (int64_t day = input_.StartDay();
         day <= input_.EndDay();
         ++day) {
      if (!obs_file_exists(day)) {
        missing_days_.insert(day);
        PLOG_("\tSource:\t" << source_name() <<
              "\tObservable:\t" << observable_name() <<
              "\tRange:\t" << input_.StartDay() << "-" << input_.EndDay() <<
              "\tMissing day: " << day <<
              std::endl);
      }
    }
    return missing_days_.size();
  }

    /** Find out the true temporal data extent.
     *
     * This function should normally not be called directly. Instead, use
     * ModelObservationSources::max_time_range, which caches the result.
     */
  DaySpan max_time_range() {
    int64_t max_search_day = 65000;
    
    std::string pathname_observation;
    int64_t startday;
    int64_t endday;

    bool found = false;
    for (int64_t day = 0;
         (!found) && (day <= max_search_day);
         ++day) {
      found = obs_file_exists(day);
      if (found) {
        startday = day;
      }
    }
    if (!found) {
      // No available observations.
      startday = -1;
      endday = -1;
    } else {
      found = false;
      for (int64_t day = max_search_day;
           (!found) && (day >= startday);
           --day) {
        found = obs_file_exists(day);
        if (found) {
          endday = day;
        }
      }
      if (!found) {
        LOG_("\tSource:\t" << source_name() <<
             "\tObservable:\t" << observable_name() <<
             "\tRange:\t" << startday << ":" << endday <<
             std::endl);
        LOG_("This should not happen: end date not found despite start date found.");
        startday = -1;
        endday = -1;
      }
    }

    if (startday >= 0) {
      return DaySpan(startday, endday);
    } else {
      return DaySpan(false);
    }
  }

 private:
  /** Retrieve a specific day of data
   */
  void raw_retrieve_day(int64_t day,
                        AnalysisInputRetrieval& retrieval);

 public:
  int64_t get_day() const {
    return day_;
  }
  int64_t size() const {
    return measurement_.rows();
  }
  void retrieve_day(int64_t day) {
    AnalysisInputRetrieval retrieval;
    day_ = day;
    if (!missing_day(day)) {
      timer_tic("Read day " + to_string(day));
      inputmanager_.RetrieveDay(retrieval,
                                input_.Source(),
                                input_.Observable(),
                                day);
      timer_toc();
      raw_retrieve_day(day, retrieval);
    } else {
      locationID_.resize(0);
      latlong_.resize(0);
      measurement_.resize(0,1);
      uncertainty_uncorrelated_.resize(0,1);
      uncertainty_locally_correlated_.resize(0,0);
      prec_uncorrelated_.resize(0,1);
    }
    //timer_toc();
  }

  // Display the currently held data
  void display()
  {
    if (day_ < 0) {
      return;
    }
    std::string name = source_name() + std::string(" ") + observable_name() + std::string(" day ") + to_string(day_);
    std::cout
        << "Measurements (" << name << ")" << std::endl;
    std::cout << "#obs = " << measurement_.size()
	      << "\t#ranges = " << uncertainty_locally_correlated_.cols()
	      << std::endl;
    for (int64_t i = 0; i < measurement_.size(); ++i) {
      std::cout
          << "#: " << i
          << ",\tlocID: " << locationID_[i]
          << ",\tlat: " << latlong_[i].latitude
          << ",\tlon: " << latlong_[i].longitude
          << ",\tobs: " << measurement_(i)
          << ",\tunc (uncorr): " << uncertainty_uncorrelated_(i);
      for (int j = 0; j < uncertainty_locally_correlated_.cols(); ++j) {
        std::cout
            << ",\tunc (corr, #" << j << "): " << uncertainty_locally_correlated_(i,j);
      }
      std::cout << std::endl;
    }
  }

  TIMER_HELPERS
};

class ModelObservationSources {
 public:
  typedef std::vector<ModelObservationSource*> sources_type;
  typedef std::pair<int64_t, int64_t> SrcObs;
  typedef std::map<SrcObs, DaySpan > ObsTimeRanges;
 private:
  TimerHierarchy* timer_;
  AnalysisInputManager inputmanager_;
  sources_type sources_;
  SpaceTimeMapper* mask_;
  DateTimeHelper dt_helper_;
  ObsTimeRanges max_time_ranges_;

 public:
  ModelObservationSources(const std::string& rawbinary_path,
                          SpaceTimeMapper* mask,
                          int64_t epoch_year,
                          TimerHierarchy* timer) :
      timer_(timer),
      inputmanager_(rawbinary_path.c_str()),
      sources_(0),
      mask_(mask),
      dt_helper_(epoch_year),
      max_time_ranges_()
  {
  }
  ~ModelObservationSources() {
    for (sources_type::iterator iter = sources_.begin();
         iter != sources_.end();
         ++iter) {
      delete *iter;
    }
    sources_.resize(0);
  }

  /** Find out the true temporal data extent.
   *
   * The result is cached to speed up subsequent queries.
   */
  const DaySpan& max_time_range(const AnalysisInput& input) {
    SrcObs input_info(input.Source(), input.Observable());
    ObsTimeRanges::iterator it = max_time_ranges_.find(input_info);
    if (it != max_time_ranges_.end()) {
      return (*it).second;
    }

    ModelObservationSource obs_source(inputmanager_, input, mask_, timer_);

    std::pair<ObsTimeRanges::iterator, bool> result =
        max_time_ranges_.insert(ObsTimeRanges::value_type(input_info,
                                                          obs_source.max_time_range()));
    return (*result.first).second;
  }

  ModelObservationSources& add(const AnalysisInput& input) {
    int64_t startday = input.StartDay();
    int64_t endday = input.EndDay();

    PLOG_("\tSource:\t" << std::string(input.sourcename[input.Source()]) <<
          "\tObservable:\t" << std::string(input.observablename[input.Observable()]) << std::endl);

    const DaySpan& time_range = max_time_range(input);

    if (time_range.empty()) {
      PLOG_("\t\tNo observations" << std::endl);
      return *this;
    }
    if (!time_range.intersects(DaySpan(startday, endday))) {
      PLOG_("\t\tNo observations in target range" << std::endl);
      return *this;
    }
    if (time_range.before(startday)) {
      PLOG_("\t\tStart day moved from " << startday << " to " << time_range.first() << std::endl);
      startday = time_range.first();
    }
    if (time_range.after(endday)) {
      PLOG_("\t\tEnd   day moved from " << endday << " to " << time_range.second() << std::endl);
      endday = time_range.second();
    }

    AnalysisInput the_input(input.Source(), input.Observable(), startday, endday);
    ModelObservationSource* the_source =
        new ModelObservationSource(inputmanager_, the_input, mask_, timer_);
    the_source->build_missing_days();
    sources_.push_back(the_source);
    inputmanager_.Specify(the_input);

    return *this;
  }

  sources_type& sources() {
    return sources_;
  }

  const DateTimeHelper& dt_helper() const {
    return dt_helper_;
  }

  void summary_log() {
    PLOG_(to_string(sources_.size()) + " data sources" << std::endl);
    for (sources_type::iterator iter = sources_.begin();
         iter != sources_.end();
         ++iter) {
      (*iter)->summary_log();
    }
  }

  /** Retrieve observations for all sources on a single day */
  void retrieve_day(int64_t day) {
    //    timer_tic("Read one day, " + to_string(sources_.size()) + " sources");
    for (sources_type::iterator iter = sources_.begin();
         iter != sources_.end();
         ++iter) {
      (*iter)->retrieve_day(day);
    }
    //    timer_toc();
  }

  /** Display currently loaded observations */
  void display_day() {
    for (sources_type::iterator iter = sources_.begin();
         iter != sources_.end();
         ++iter) {
      (*iter)->display();
    }
  }
  
  TIMER_HELPERS
};

class ModelComponent {
};

class ModelComponents {
};


#endif
