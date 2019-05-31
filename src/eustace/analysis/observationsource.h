#ifndef _EUSTACE_ANALYSIS_OBSERVATIONSOURCE_H_
#define _EUSTACE_ANALYSIS_OBSERVATIONSOURCE_H_

#include <stdint.h>
#include <math.h>
#include <memory.h>
#include <vector>

namespace EUSTACE
{
  class LocationLookupUUID
  {
  public:
    static const size_t sizebytes = 16;

  public:
    LocationLookupUUID()
    {
      memset(bytes, 0, LocationLookupUUID::sizebytes);
    }

    LocationLookupUUID(const uint8_t* _bytes)
    {
      memcpy(bytes, _bytes, LocationLookupUUID::sizebytes);
    }

  public:
    uint8_t* Bytes()
    { return bytes; }

    const uint8_t* Bytes() const
    { return bytes; }

  private:
    uint8_t bytes[LocationLookupUUID::sizebytes];
  };

  class LocationEntry
  {
  public:
    double latitude;
    double longitude;
  };

  class LocationLookup
  {
  public:
    LocationLookup()
    {
    }

    const LocationLookupUUID& ID() const
    { return id; }

    LocationLookupUUID& ID()
    { return id; }

    const std::vector<LocationEntry>& Lookup() const
    { return lookup; }

    std::vector<LocationEntry>& Lookup()
    { return lookup; }

  private:
    LocationLookupUUID id;
    std::vector<LocationEntry> lookup;
  };

  class ObservationEntries
  {
  public:
    ObservationEntries() :
      num_locally_correlated(0),
      num_entries(0)
    {
    }

    void SetDimensions(size_t _num_locally_correlated, size_t _num_entries)
    {
      num_entries = _num_entries;
      num_locally_correlated = _num_locally_correlated;
      entries.resize(EntrySizeBytes(_num_locally_correlated) * _num_entries);
    }

    size_t NumEntries() const
    { return num_entries; }

    static size_t EntrySizeBytes(size_t num_locally_correlated)
    { return 8 + 8 + 8 + 8*num_locally_correlated; }

    size_t TotalSizeBytes() const
    { return entries.size(); }

    uint8_t* Bytes()
    { return &entries[0]; }

    const uint8_t* Bytes() const
    { return &entries[0]; }

    LocationLookupUUID& LocationLookupID()
    { return location_lookup_uuid; }

    const LocationLookupUUID& LocationLookupID() const
    { return location_lookup_uuid; }

    size_t NumLocalCorrelationRanges() const
    { return num_locally_correlated; }

  private:
    size_t num_locally_correlated;
    size_t num_entries;
    LocationLookupUUID location_lookup_uuid;
    std::vector<uint8_t> entries;
  };

  class ObservationIterator
  {
  public:
    ObservationIterator(const ObservationEntries& _entries) :
      entries(&_entries),
      elementsize(ObservationEntries::EntrySizeBytes(_entries.NumLocalCorrelationRanges())),
      index(0),
      pointer(_entries.Bytes())
    {
    }

    uint64_t Identifier() const
    {
      return *(uint64_t*)pointer;
    }

    double Measurement() const
    {
      return *(double*)(pointer+8);
    }

    double UncertaintyUncorrelated() const
    {
      return *(double*)(pointer+16);
    }

    double UncertaintyLocallyCorrelated(size_t local_index) const
    {
      return ((double*)(pointer+24))[local_index];
    }

    bool HasElement() const
    {
      return (index < entries->NumEntries());
    }

  public:

    void Next()
    {
      index++;
      pointer += elementsize;
    }

  private:
    const ObservationEntries* entries;
    size_t elementsize;
    size_t index;
    const uint8_t* pointer;
  };

}

#endif 
