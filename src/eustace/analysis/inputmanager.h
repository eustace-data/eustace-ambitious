#ifndef _EUSTACE_ANALYSIS_INPUTMANAGER_H_
#define _EUSTACE_ANALYSIS_INPUTMANAGER_H_

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <stdint.h>
#include "observationsource.h"
#include "fileio/rawbinary.h"

namespace EUSTACE
{
  class AnalysisInput
  {
  public:

    enum sourcetype
    {
      source_satellite_land,
      source_satellite_ocean,
      source_satellite_ice,
      source_satellite_lake,
      source_insitu_land,
      source_insitu_ocean,
      num_sources
    };

    enum observabletype
    {
      observable_tmin,
      observable_tmax,
      observable_tmean,
      num_observables
    };

  public:

    static const char* sourcename[num_sources];
    static const char* observablename[num_observables];

  public:
    
    AnalysisInput(enum sourcetype _source, 
		  enum observabletype _observable,
		  int64_t _startday,
		  int64_t _endday) :
      source(_source),
      observable(_observable),
      startday(_startday),
      endday(_endday)
    {
    }

    enum sourcetype Source() const
    { return source; }

    enum observabletype Observable() const
    { return observable; }

    int64_t StartDay() const
    { return startday; }

    int64_t EndDay() const
    { return endday; }

  private:
    enum sourcetype source;
    enum observabletype observable;
    int64_t startday;
    int64_t endday;
  };

  class AnalysisInputRetrieval
  {
    friend class AnalysisInputManager;
    friend class AnalysisInputIterator;

  public:    

    AnalysisInputRetrieval();    

    ~AnalysisInputRetrieval();

    size_t NumObservations() const
    {
      return observations.NumEntries();
    }

    size_t NumLocalCorrelationRanges() const
    {
      return local_correlation_ranges->size();
    }

    double LocalCorrelationRange(size_t local_index) const
    {
      return (*local_correlation_ranges)[local_index];
    }

    const LocationEntry& ComputeLocationLookup(uint64_t location_id) const
    { 
      return location_lookup->Lookup()[location_id];
    }

  private:

    std::vector<double>* local_correlation_ranges;
    LocationLookup* location_lookup;
    LocationLookup mobile_location_lookup;
    ObservationEntries observations;
  };

  class AnalysisInputRetrievalCollection
  {
  public:

    AnalysisInputRetrievalCollection();

    virtual ~AnalysisInputRetrievalCollection();

    AnalysisInputRetrieval* NewRetrieval();

    const std::vector<AnalysisInputRetrieval*>& Retrievals() const
    { return retrievals; }

  private:

    std::vector<AnalysisInputRetrieval*> retrievals;
  };

  class AnalysisInputIterator : public ObservationIterator
  {
  public:
    AnalysisInputIterator(const AnalysisInputRetrieval& _retrieval) :
      ObservationIterator(_retrieval.observations),
      retrieval(&_retrieval)
    {
    }

    const LocationEntry& Location() const
    {
      return retrieval->ComputeLocationLookup(Identifier());
    }

  private:
    const AnalysisInputRetrieval* retrieval;
  };


  template<class DataT> class AnalysisInputCache
  {
  public:
    AnalysisInputCache()
    {
    }

    virtual ~AnalysisInputCache()
    {
      typename std::map<std::string, DataT*>::iterator iter;
      for (iter = cache.begin(); iter != cache.end(); iter++)
	delete iter->second;
      cache.clear();
    }

  public:

    virtual DataT* Retrieve(const std::string& pathname)
    {
      typename std::map<std::string, DataT*>::iterator iter;
      iter = cache.find(pathname);
      if (iter != cache.end())
      	return iter->second;
      DataT* newitem = new DataT;
      std::ifstream inputstream(pathname.c_str(), std::ios::binary);
      Read(*newitem, inputstream);
      inputstream.close();
      cache[pathname] = newitem;
      return newitem;
    }

  protected:

    virtual void Read(DataT& result, std::istream& inputstream) = 0;

  private:
    std::map<std::string, DataT*> cache;
  };

  class CacheLocationLookup : public AnalysisInputCache<LocationLookup>
  {
  protected:
    virtual void Read(LocationLookup& result, std::istream& inputstream);
  };

  class CacheLocalCorrelationRange : public AnalysisInputCache<std::vector<double> >
  {
  protected:
    virtual void Read(std::vector<double>& result, std::istream& inputstream);
  };

  class AnalysisInputManager
  {
  public:

    AnalysisInputManager(const char* _basepath);

    void Specify(const AnalysisInput& input);

    void RetrieveDay(AnalysisInputRetrieval& retrieval,
		     enum AnalysisInput::sourcetype source,
		     enum AnalysisInput::observabletype observable,
		     int64_t daynumber);

    void RetrieveLocations(AnalysisInputRetrievalCollection& collection, 
			   enum AnalysisInput::sourcetype source,
			   enum AnalysisInput::observabletype observable,
			   std::vector<uint64_t>& location_id_subset);

    size_t TotalLocations(enum AnalysisInput::sourcetype source) const;

    size_t TotalDailyMobileLocations(enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable, int64_t daynumber) const;

    size_t TotalObservations(enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable);

    size_t TotalDailyObservations(enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable, int64_t daynumber);

  public:

    const std::string& BasePath() const
    { return basepath; }

    const std::vector<AnalysisInput>& Inputs() const
    { return inputs; }

    void PathNameObservation(std::string& pathname, enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable, int64_t daynumber) const;

    void PathNameLocationLookup(std::string& pathname, enum AnalysisInput::sourcetype source) const;

    void PathNameLocationPointers(std::string& pathname, enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable) const;

    void PathNameMobileLocationLookup(std::string& pathname, enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable, int64_t daynumber) const;
    
    void PathNameLocalCorrelationRange(std::string& pathname, enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable) const;
  
  private:
    std::map<enum AnalysisInput::sourcetype, std::string> staticlocation;
    std::string basepath;
    std::vector<AnalysisInput> inputs;
    CacheLocationLookup cache_locationlookup;
    CacheLocalCorrelationRange cache_localcorrelationrange;
  };

}

#endif // #ifndef _EUSTACE_ANALYSIS_INPUTMANAGER_H_
