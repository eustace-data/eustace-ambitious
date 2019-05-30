#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(ncdf_helper_cpp, "$Revision: 1293 $")

#include "ambitious/ncdf/ncdf_helper.hpp"


#ifndef USE_PNETCDF
#else
MPI_Datatype_nc nc_type_to_mpi(nc_type type) {
  switch (type) {
    case NC_BYTE: return MPI_BYTE; break;
    case NC_CHAR: return MPI_CHAR; break;
    case NC_SHORT: return MPI_SHORT; break;
    case NC_INT: return MPI_INT; break;
    case NC_FLOAT: return MPI_FLOAT; break;
    case NC_DOUBLE: return MPI_DOUBLE; break;
    case NC_UBYTE: return MPI_UINT8_T; break;
    case NC_USHORT: return MPI_UNSIGNED_SHORT; break;
    case NC_UINT: return MPI_UNSIGNED; break;
    case NC_INT64: return MPI_INT64_T; break;
    case NC_UINT64: return MPI_UINT64_T; break;
  }
  return MPI_DATATYPE_NULL;
}
#endif
