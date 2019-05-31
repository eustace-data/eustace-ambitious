  // Use of location pointer files to read time sequences of data

#include "locationpointers.h"

namespace EUSTACE
{
  // Internal offset within location pointer file
  // - not to be confused with LocationPointerRecord which is a record of the offset within an observation file
  class LocalOffsetRecord
  {
  public:
    uint64_t fileoffset;
    uint64_t count;
  };


  LocationPointersRawBinaryReader::LocationPointersRawBinaryReader() :
      RawBinaryReader()
  {
  }

  const char* LocationPointersRawBinaryReader::FormatID() const
  {
    return "EUSTACELOCP00001";
  }

  void LocationPointersRawBinaryReader::ReadAndCheckExtendedHeader(LocationLookupUUID& uuid, uint64_t& location_count, std::istream& inputstream) const throw(FormatException)
  {
    // check format identifier is valid
    RawBinaryReader::ReadAndCheckHeader(inputstream);

    // get uuid of location lookup
    inputstream.read((char*)uuid.Bytes(), LocationLookupUUID::sizebytes);

    // get count of location ids
    inputstream.read((char*)&location_count, sizeof(uint64_t));
  }

  void LocationPointersRawBinaryReader::Read(std::vector< std::vector<LocationPointerRecord> >& record, std::istream& inputstream, const std::vector<uint64_t>& location_id_subset) const throw(FormatException)
  {
    // read first part of header
    LocationLookupUUID uuid;
    uint64_t location_count(0);
    ReadAndCheckExtendedHeader(uuid, location_count, inputstream);
    
    // start of local offsets table
    size_t local_offset_table_start = inputstream.tellg();

    // read daily data from the specified subset of locations
    record.resize(location_id_subset.size());
    std::vector< std::vector<LocationPointerRecord> >::iterator record_iterator(record.begin());
    for (std::vector<uint64_t>::const_iterator location_id_iterator(location_id_subset.begin()); location_id_iterator != location_id_subset.end(); location_id_iterator++, record_iterator++)
    {
      // entry in local offsets table for this id
      LocalOffsetRecord localoffset;
      memset(&localoffset, 0, sizeof(LocalOffsetRecord));
      inputstream.seekg(local_offset_table_start + size_t(*location_id_iterator)*sizeof(LocalOffsetRecord), std::istream::beg);
      inputstream.read((char*)&localoffset, sizeof(LocalOffsetRecord));

      // corresponding output element
      std::vector<LocationPointerRecord>& current_record = *record_iterator;

      // seek to specified offset that has data for this id
      inputstream.seekg(localoffset.fileoffset, std::istream::beg);

      // allocate buffer
      current_record.resize(localoffset.count);

      // read
      inputstream.read((char*)&current_record[0], localoffset.count*sizeof(LocationPointerRecord));
    }
  }
}
