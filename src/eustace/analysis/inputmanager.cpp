#include "inputmanager.h"
#include "fileio/locationpointers.h"
#include "eustace/definitions.h"

namespace EUSTACE
{
  const char* AnalysisInput::sourcename[AnalysisInput::num_sources] =
  {
    "surfaceairmodel_land",
    "surfaceairmodel_ocean",
    "surfaceairmodel_ice",
    "surfaceairmodel_lake",
    "insitu_land",
    "insitu_ocean"
  };

  const char* AnalysisInput::observablename[AnalysisInput::num_observables] =
  {
    "Tmin",
    "Tmax",
    "Tmean"
  };

  AnalysisInputRetrieval::AnalysisInputRetrieval() :
    local_correlation_ranges(0),
    location_lookup(0)
  {
  }

  AnalysisInputRetrieval::~AnalysisInputRetrieval()
  {
  }

  AnalysisInputRetrievalCollection::AnalysisInputRetrievalCollection()
  {
  }

  AnalysisInputRetrievalCollection::~AnalysisInputRetrievalCollection()
  {
    for (std::vector<AnalysisInputRetrieval*>::iterator iter(retrievals.begin()); iter != retrievals.end(); iter++)
    {
      delete *iter;
    }
    retrievals.clear();
  }

  AnalysisInputRetrieval* AnalysisInputRetrievalCollection::NewRetrieval()
  {
    AnalysisInputRetrieval* newretrieval = new AnalysisInputRetrieval;
    retrievals.push_back(newretrieval);
    return newretrieval;
  }

  AnalysisInputManager::AnalysisInputManager(const char* _basepath) :
    basepath(_basepath)
  {
    staticlocation[AnalysisInput::source_satellite_land] = "satellite";
    staticlocation[AnalysisInput::source_satellite_ocean] = "satellite";
    staticlocation[AnalysisInput::source_satellite_ice] = "satellite";
    staticlocation[AnalysisInput::source_satellite_lake] = "lake";
    staticlocation[AnalysisInput::source_insitu_land] = "insitu_land";
  }

  void AnalysisInputManager::Specify(const AnalysisInput& input)
  {
    inputs.push_back(input); 
  }

  size_t AnalysisInputManager::TotalLocations(enum AnalysisInput::sourcetype source) const
  {    
    // init to zero
    size_t num_entries(0);

    // look at fixed location types
    if (staticlocation.find(source) != staticlocation.end())
    {
      // Pathname to fixed location lookup
      std::string pathname_location_lookup;
      PathNameLocationLookup(pathname_location_lookup, source);

      // Open stream and count
      std::ifstream inputstream(pathname_location_lookup.c_str(), std::ios::binary);
      LocationLookupUUID uuid;
      LocationLookupRawBinaryReader().ReadAndCheckExtendedHeader(inputstream, uuid, num_entries);
      inputstream.close();
    }

    return num_entries;
  }

  size_t AnalysisInputManager::TotalDailyMobileLocations(enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable, int64_t daynumber) const
  {
    // init to zero
    size_t num_entries(0);

    // make pathname to daily mobile location lookup
    std::string pathname_mobile_location_lookup;
    PathNameMobileLocationLookup(pathname_mobile_location_lookup, source, observable, daynumber);

    // Read into our local mobile lookup structure
    std::ifstream inputstream(pathname_mobile_location_lookup.c_str(), std::ios::binary);
    LocationLookupUUID uuid;
    LocationLookupRawBinaryReader().ReadAndCheckExtendedHeader(inputstream, uuid, num_entries);
    inputstream.close();

    // return result
    return num_entries;
  }

  size_t AnalysisInputManager::TotalDailyObservations(enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable, int64_t daynumber)
  {
    // initialise to zero
    size_t dailycount(0);

    // correlation range numbers
    std::string pathname_local_correlation_range;
    PathNameLocalCorrelationRange(pathname_local_correlation_range, source, observable);
    size_t num_local_correlation_ranges = cache_localcorrelationrange.Retrieve(pathname_local_correlation_range)->size();

    // path to these obs
    std::string pathname_observation;
    PathNameObservation(pathname_observation, source, observable, daynumber);

    // reader object
    ObservationRawBinaryReader reader(num_local_correlation_ranges);

    // read count
    std::ifstream inputstream(pathname_observation.c_str(), std::ios::binary);
    LocationLookupUUID uuid;
    reader.ReadAndCheckExtendedHeader(inputstream, uuid, dailycount);
    inputstream.close();

    // return result
    return dailycount;
  }

  size_t AnalysisInputManager::TotalObservations(enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable)
  {
    // initialise to zero
    size_t total(0);

    // search for matching source info
    for (std::vector<AnalysisInput>::const_iterator inputiter(inputs.begin()); inputiter != inputs.end(); inputiter++)
    {
      if (inputiter->Source() == source && inputiter->Observable() == observable)
      {
	// source found - add obs from all days
	for (int64_t daynumber = inputiter->StartDay(); daynumber <= inputiter->EndDay(); daynumber++)
	{
	  // get count for this day
	  size_t dailycount = TotalDailyObservations(source, observable, daynumber);

	  // accumulate
	  total += dailycount;
	}
      }
    }

    // return result
    return total;
  }
 
  void AnalysisInputManager::RetrieveDay(AnalysisInputRetrieval& retrieval, 
				    enum AnalysisInput::sourcetype source,
				    enum AnalysisInput::observabletype observable,
				    int64_t daynumber)
				     
  {
    // Correlation range numbers
    std::string pathname_local_correlation_range;
    PathNameLocalCorrelationRange(pathname_local_correlation_range, source, observable);
    retrieval.local_correlation_ranges = cache_localcorrelationrange.Retrieve(pathname_local_correlation_range);

    // Use cached static location file if it's a static location or assume mobile otherwise
    if (staticlocation.find(source) != staticlocation.end())
    {
      // Pathname to fixed location lookup (for all time)
      std::string pathname_location_lookup;
      PathNameLocationLookup(pathname_location_lookup, source);

      // Set retrieval pointer to cache of fixed lookup
      retrieval.location_lookup = cache_locationlookup.Retrieve(pathname_location_lookup);
    }
    else
    {
      // make pathname to daily mobile location lookup
      std::string pathname_mobile_location_lookup;
      PathNameMobileLocationLookup(pathname_mobile_location_lookup, source, observable, daynumber);

      // Read into our local mobile lookup structure
      std::ifstream inputstream(pathname_mobile_location_lookup.c_str(), std::ios::binary);
      LocationLookupRawBinaryReader mobile_location_reader;
      mobile_location_reader.Read(retrieval.mobile_location_lookup, inputstream);
      inputstream.close();

      // Set retrieval pointer to newly read mobile locations
      retrieval.location_lookup = &retrieval.mobile_location_lookup;
    }

    // Read observations
    std::string pathname_observation;
    PathNameObservation(pathname_observation, source, observable, daynumber);
    ObservationRawBinaryReader reader(retrieval.local_correlation_ranges->size());
    std::ifstream inputstream(pathname_observation.c_str(), std::ios::binary);
    reader.Read(retrieval.observations, inputstream);
    inputstream.close();
  }

  void AnalysisInputManager::RetrieveLocations(AnalysisInputRetrievalCollection& collection, 
					       enum AnalysisInput::sourcetype source,
					       enum AnalysisInput::observabletype observable,
					       std::vector<uint64_t>& location_id_subset)
  {
    // Only use static locations here
    if (staticlocation.find(source) == staticlocation.end())
      return;

    // Correlation range numbers (lazy-loaded cache)
    std::string pathname_local_correlation_range;
    PathNameLocalCorrelationRange(pathname_local_correlation_range, source, observable);
    std::vector<double>* local_correlation_ranges = cache_localcorrelationrange.Retrieve(pathname_local_correlation_range);

    // Fixed location lookup for all time (lazy loaded cache)
    std::string pathname_location_lookup;
    PathNameLocationLookup(pathname_location_lookup, source);
    LocationLookup* location_lookup = cache_locationlookup.Retrieve(pathname_location_lookup);

    // Get file offset information for all ids
    std::vector< std::vector<LocationPointerRecord> > pointer_records;
    {
      std::string pathname_location_pointers;
      PathNameLocationPointers(pathname_location_pointers, source, observable);
      std::ifstream pointerstream(pathname_location_pointers.c_str(), std::ios::binary);
      LocationPointersRawBinaryReader().Read(pointer_records, pointerstream, location_id_subset);
    }

    // number of locally correlated uncertainty ranges
    size_t num_correlation_ranges = local_correlation_ranges->size();

    // location id size in array
    size_t location_id_element_size = sizeof(uint64_t);
    
    // observation entry size
    size_t observation_data_element_size = ObservationEntries::EntrySizeBytes(num_correlation_ranges) - location_id_element_size;

    // There is an array of pointers for each location
    std::vector<uint64_t>::const_iterator location_id_iterator(location_id_subset.begin());
    for (std::vector< std::vector<LocationPointerRecord> >::const_iterator record_iterator(pointer_records.begin());
	 record_iterator != pointer_records.end();
	 record_iterator++, location_id_iterator++)
    {
      // The array of pointers for this location - one for each day on which there is data
      const std::vector<LocationPointerRecord>& location_pointer_records = *record_iterator;

      // Location ID we are looking at
      uint64_t location_id = *location_id_iterator;

      // create a retrieval object to store results
      AnalysisInputRetrieval* retrieval = collection.NewRetrieval();
      retrieval->local_correlation_ranges = local_correlation_ranges;
      retrieval->location_lookup = location_lookup;
      retrieval->observations.SetDimensions(num_correlation_ranges, location_pointer_records.size());

      // pointer to raw data for reading
      char* observation_pointer = (char*) retrieval->observations.Bytes();

      // read data for each day
      for (std::vector<LocationPointerRecord>::const_iterator time_iterator(location_pointer_records.begin());
	   time_iterator != location_pointer_records.end();
	   time_iterator++)
      {

	// the day number is used as the identifier in this context
	*(uint64_t*)observation_pointer = time_iterator->daynumber;

	// increment output pointer ready to store observation data
	observation_pointer += location_id_element_size;

	// open file
	// NOTE: if quicker could cache input streams in a map indexed by day number
	std::string pathname_observation;
	PathNameObservation(pathname_observation, source, observable, time_iterator->daynumber);
	std::ifstream obsstream(pathname_observation.c_str(), std::ios::binary);

	// seek to offset expressed in pointers file
	obsstream.seekg(time_iterator->fileoffset, std::ifstream::beg);	

	// read information
	obsstream.read(observation_pointer, observation_data_element_size);

	// and increment output pointer ready to read next element
	observation_pointer += observation_data_element_size;
      }
    }
  }


  void AnalysisInputManager::PathNameObservation(std::string& pathname, enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable, int64_t daynumber) const
  {
    std::string datestring("");
    CalendarDay(TimeBaseDays(EPOCH).Time(daynumber)).Text(datestring);

    pathname = basepath + std::string("/") + 
      AnalysisInput::sourcename[source] + std::string("/") +
      datestring.substr(0, 4)+ std::string("/") +
      AnalysisInput::sourcename[source] + std::string("_") +
      AnalysisInput::observablename[observable] + std::string("_") +
      datestring + std::string(".bin");
  }

  void AnalysisInputManager::PathNameLocationPointers(std::string& pathname, enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable) const
  {
    pathname = basepath + std::string("/") + 
      AnalysisInput::sourcename[source] + std::string("_") +
      AnalysisInput::observablename[observable] + std::string("_locationpointers.bin");
  }

  void AnalysisInputManager::PathNameLocationLookup(std::string& pathname, enum AnalysisInput::sourcetype source) const
  {
    // Retrieve source name (static source)
    std::string sourcename = staticlocation.find(source)->second;

    // Build path name
    pathname = basepath + std::string("/locationlookup_") + sourcename  + std::string(".bin");
  }
  
  void AnalysisInputManager::PathNameMobileLocationLookup(std::string& pathname, enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable, int64_t daynumber) const
  {
    // Convert daynumber to date string
    std::string datestring("");
    CalendarDay(TimeBaseDays(EPOCH).Time(daynumber)).Text(datestring);

    // Build path name containing date for mobile location lookup
    pathname = basepath + std::string("/") +
      AnalysisInput::sourcename[source] + std::string("/") +
      datestring.substr(0, 4)+ std::string("/") +
      AnalysisInput::sourcename[source] + std::string("_") +
      AnalysisInput::observablename[observable] + std::string("_mobilelocations_") +
      datestring + std::string(".bin");
  }

  void AnalysisInputManager::PathNameLocalCorrelationRange(std::string& pathname, enum AnalysisInput::sourcetype source, enum AnalysisInput::observabletype observable) const
  {
    pathname = basepath + std::string("/") +
      AnalysisInput::sourcename[source] + std::string("_") +
      AnalysisInput::observablename[observable] + std::string("_localcorrelationranges.bin");
  }
 
  void CacheLocationLookup::Read(LocationLookup& result, std::istream& inputstream)
  {
    LocationLookupRawBinaryReader reader;
    reader.Read(result, inputstream);
  }

  void CacheLocalCorrelationRange::Read(std::vector<double>& result, std::istream& inputstream)
  {
    LocalCorrelationRangeRawBinaryReader reader;
    reader.Read(result, inputstream);
  }


}
