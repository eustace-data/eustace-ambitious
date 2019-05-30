#include "rawbinary.h"
#include <iostream>

namespace EUSTACE
{
  RawBinaryReader::RawBinaryReader()
  {
  }

  void RawBinaryReader::ReadAndCheckHeader(std::istream& inputstream) const throw(FormatExceptionBadHeader)
  {
    // file format ID (should match the 16-byte string returned by FormatID() )
    char file_format_id[16];
    inputstream.read(file_format_id, sizeof(file_format_id));

    // read and ignore code identifier
    char code_id[16];
    inputstream.read(code_id, sizeof(code_id));

    // check format ID matches
    if (strncmp(file_format_id, FormatID(), 16) != 0)
      throw FormatExceptionBadHeader();
  }
    
  size_t RawBinaryReader::NumEntriesRemaining(std::istream& inputstream, size_t entrysize) throw(FormatExceptionUnexpectedDataSize)
  {
    // seek to end and back to determine bytes remaining to end of file
    size_t startposition = inputstream.tellg();
    inputstream.seekg(0, std::ios::end);
    size_t endposition = inputstream.tellg();
    inputstream.seekg(startposition, std::ios::beg);
    size_t databytes = endposition - startposition;

    // check is multiple of expected element size
    if (databytes % entrysize != 0)
      throw FormatExceptionUnexpectedDataSize();

    // compute num entries
    size_t num_entries = databytes / entrysize;
    return num_entries;
  }

  RawBinaryLocationCategoryReader::RawBinaryLocationCategoryReader() :
    RawBinaryReader()
  {
  }

  void RawBinaryLocationCategoryReader::ReadAndCheckExtendedHeader(std::istream& inputstream, LocationLookupUUID& uuid, size_t& num_entries) const throw(FormatException)
  {
    // make sure basic header is ok
    ReadAndCheckHeader(inputstream);

    // get uuid of location lookup
    inputstream.read((char*)uuid.Bytes(), LocationLookupUUID::sizebytes);

    // size of one entry
    size_t entrysize = EntrySize();

    // find number of entries present
    num_entries = RawBinaryReader::NumEntriesRemaining(inputstream, entrysize);
  }

  ObservationRawBinaryReader::ObservationRawBinaryReader(uint64_t _num_correlation_ranges) :
      RawBinaryLocationCategoryReader(),
      num_correlation_ranges(_num_correlation_ranges)
  {
  }

  const char* ObservationRawBinaryReader::FormatID() const
  {
    return "EUSTACEOBSN00001";
  }

  size_t ObservationRawBinaryReader::EntrySize() const
  {
    return ObservationEntries::EntrySizeBytes(num_correlation_ranges);
  }

  void ObservationRawBinaryReader::Read(ObservationEntries& obs, std::istream& inputstream) const throw(FormatException)
  {
    // make sure header ok
    size_t num_entries(0);
    ReadAndCheckExtendedHeader(inputstream, obs.LocationLookupID(), num_entries);

    // allocate and read entries
    obs.SetDimensions(num_correlation_ranges, num_entries);
    inputstream.read((char*)obs.Bytes(), num_entries * EntrySize());
  }

  LocationLookupRawBinaryReader::LocationLookupRawBinaryReader() :
    RawBinaryLocationCategoryReader()
  {
  }

  const char* LocationLookupRawBinaryReader::FormatID() const
  {
    return "EUSTACELOCN00001";
  }

  size_t LocationLookupRawBinaryReader::EntrySize() const
  {
    return sizeof(LocationEntry);
  }

  void LocationLookupRawBinaryReader::Read(LocationLookup& loc, std::istream& inputstream) const throw(FormatException)
  {
    // check header and get location UUID
    size_t num_entries(0);
    ReadAndCheckExtendedHeader(inputstream, loc.ID(), num_entries);

    // allocate and read
    loc.Lookup().resize(num_entries);
    inputstream.read((char*)&(loc.Lookup()[0]), num_entries * EntrySize());
  }

  LocalCorrelationRangeRawBinaryReader::LocalCorrelationRangeRawBinaryReader() :
    RawBinaryReader()
  {
  }

  const char* LocalCorrelationRangeRawBinaryReader::FormatID() const
  {
    return "EUSTACECORN00001";
  }

  void LocalCorrelationRangeRawBinaryReader::Read(std::vector<double>& ranges, std::istream& inputstream) const throw(FormatException)
  {
    // make sure header ok
    ReadAndCheckHeader(inputstream);

    // size of one entry
    size_t entrysize = sizeof(double);

    // find number of entries present
    size_t num_entries = RawBinaryReader::NumEntriesRemaining(inputstream, entrysize);
    
    // allocate and read
    ranges.resize(num_entries);
    inputstream.read((char *)&ranges[0], num_entries * entrysize);
  }
}
