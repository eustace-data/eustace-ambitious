 // Use of location pointer files to read time sequences of data

#ifndef _EUSTACE_ANALYSIS_FILEIO_LOCATIONPOINTERS_H_
#define _EUSTACE_ANALYSIS_FILEIO_LOCATIONPOINTERS_H_

#include "rawbinary.h"

namespace EUSTACE
{
  class LocationPointerRecord
  {
  public:
    int64_t daynumber;
    uint64_t fileoffset;
  };

  class LocationPointersRawBinaryReader : public RawBinaryReader
  {
  public:   
    LocationPointersRawBinaryReader();

    virtual const char* FormatID() const;

    virtual void ReadAndCheckExtendedHeader(LocationLookupUUID& uuid, uint64_t& location_count, std::istream& inputstream) const throw(FormatException);

    virtual void Read(std::vector< std::vector<LocationPointerRecord> >& record, std::istream& inputstream, const std::vector<uint64_t>& location_id_subset) const throw(FormatException);

  };

}

#endif // #ifdef _EUSTACE_ANALYSIS_FILEIO_LOCATIONPOINTERS_H_

