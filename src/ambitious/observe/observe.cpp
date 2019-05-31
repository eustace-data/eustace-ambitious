#include "ambitious/observe/observe.hpp"

#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(observe_cpp, "$Revision: 1290 $")





void ModelObservationSource::raw_retrieve_day(int64_t day,
                        AnalysisInputRetrieval& retrieval)
{
  size_t i = 0;
  AnalysisInputIterator iterator(retrieval);
  std::vector<int64_t> mask_count;
  int64_t found = 0;
  if (mask_) {
    timer_tic("Mask locations");
    mask_count.resize(retrieval.NumObservations());
    bool remember = !is_mobile();
    bool cellular = is_cell();
    for (i=0; i < retrieval.NumObservations(); ++i) {
      const LocationEntry& loc = iterator.Location();
      
      std::vector<Eigen::Triplet<double> > output(0);
      SpaceTimeMapper::offset_type offset(0, 0);
      mask_count[i] = mask_->mapping(
          SpaceTimeMapper::Point(loc.latitude,
                                 loc.longitude,
                                 day,
                                 cellular),
          true, // solar_time
          offset, // target row index
          remember, // remember
          NULL, // No extra weighting
          output);
      if (mask_count[i] > 0) {
        ++found;
      }
      iterator.Next();
    }
    timer_toc();
  } else {
    found = retrieval.NumObservations();
  }

  PLOG_("\tRetained/retrieved: " << found << "/" << retrieval.NumObservations() << std::endl);

  locationID_.resize(found);
  latlong_.resize(found);
  measurement_.resize(found,1);
  uncertainty_uncorrelated_.resize(found,1);
  uncertainty_locally_correlated_.resize(found,
                                         retrieval.NumLocalCorrelationRanges());
  prec_uncorrelated_.resize(found,1);

  if (found > 0) {
    //    timer_tic("Copy observations");
    size_t j = 0;
    iterator = AnalysisInputIterator(retrieval);
    for (i = 0; i < retrieval.NumObservations(); ++i) {
      if (!mask_ || (mask_count[i] > 0)) {
        locationID_[j] = iterator.Identifier();
        latlong_[j] = retrieval.ComputeLocationLookup(locationID_[j]);
        measurement_(j) = iterator.Measurement();
        uncertainty_uncorrelated_(j) = iterator.UncertaintyUncorrelated();
        for (size_t k = 0; k < retrieval.NumLocalCorrelationRanges(); ++k) {
          uncertainty_locally_correlated_(j,k) = iterator.UncertaintyLocallyCorrelated(k);
        }
        ++j;
      }
      iterator.Next();
    }
    prec_uncorrelated_ =
        1.0 / uncertainty_uncorrelated().array() / uncertainty_uncorrelated().array();
    //    timer_toc();
  }
  //  std::cout
  //    << "NumObservations = " << retrieval.NumObservations()
  //    << ",\tNumber of read observations = " << i
  //    << ",\tmeasurement_.size() = " << measurement_.size()
  //    << std::endl;
}
