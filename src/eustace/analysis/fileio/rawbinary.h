// Reading of raw binary observations

#ifndef _EUSTACE_ANALYSIS_FILEIO_RAWBINARY_H_
#define _EUSTACE_ANALYSIS_FILEIO_RAWBINARY_H_

#include "eustace/analysis/observationsource.h"
#include <stdint.h>
#include <math.h>
#include <memory.h>
#include <fstream>
#include <vector>

namespace EUSTACE
{
  class FormatException
  {
  };

  class FormatExceptionBadHeader : public FormatException
  {
  };

  class FormatExceptionUnexpectedDataSize : public FormatException
  {
  };

  class RawBinaryReader
  {
  public:
    RawBinaryReader();

  public:

    virtual const char* FormatID() const = 0;

  protected:

    virtual void ReadAndCheckHeader(std::istream& inputstream) const throw(FormatExceptionBadHeader);

    static size_t NumEntriesRemaining(std::istream& inputstream, size_t entrysize) throw(FormatExceptionUnexpectedDataSize);

  };

  class RawBinaryLocationCategoryReader : public RawBinaryReader
  {
  public:

    RawBinaryLocationCategoryReader();

  public:

    virtual size_t EntrySize() const = 0;

    virtual void  ReadAndCheckExtendedHeader(std::istream& inputstream, LocationLookupUUID& uuid, size_t& num_entries) const throw(FormatException);    
  };

  class ObservationRawBinaryReader : public RawBinaryLocationCategoryReader
  {
  public:   
    ObservationRawBinaryReader(uint64_t _num_correlation_ranges);

    virtual const char* FormatID() const;

    virtual size_t EntrySize() const;

    virtual void Read(ObservationEntries& obs, std::istream& inputstream) const throw(FormatException);

  private:
    uint64_t num_correlation_ranges;    
  };

  class LocationLookupRawBinaryReader : public RawBinaryLocationCategoryReader
  {
  public:
    LocationLookupRawBinaryReader();

  public:
    
    virtual const char* FormatID() const;

    virtual size_t EntrySize() const;

    virtual void Read(LocationLookup& loc, std::istream& inputstream) const throw(FormatException);
    
  };

  class LocalCorrelationRangeRawBinaryReader : public RawBinaryReader
  {
  public:
    LocalCorrelationRangeRawBinaryReader();

  public:
    
    virtual const char* FormatID() const;

    virtual void Read(std::vector<double>& ranges, std::istream& inputstream) const throw(FormatException);
    
  };
}

#endif // #ifdef _EUSTACE_ANALYSIS_FILEIO_RAWBINARY_H_
