// 

#include "mesh.h"
#include <netcdf.h>
#include <iostream>

namespace EUSTACE
{
  class MeshNetCDFHelper
  {
  public:
    MeshNetCDFHelper(const char* filename) throw (MeshReaderException);

    virtual ~MeshNetCDFHelper();

  public:

    void GetDimensions(std::vector<size_t>& dimensions, const char* name) const throw (MeshReaderException);

    template<class DataT> void ReadData(std::vector<DataT>& contents, const char* name) const throw (MeshReaderException);

  private:
    int ncid;
  };

  MeshNetCDFHelper::MeshNetCDFHelper(const char* filename) throw (MeshReaderException) :
    ncid(0)
  {
    // Open it and set member file id
    int errorcode = nc_open(filename, NC_NOWRITE, &ncid);
    if (errorcode)
    {
      throw MeshReaderException();
    }
  }

  MeshNetCDFHelper::~MeshNetCDFHelper()
  {
    if (ncid)
    {
      nc_close(ncid);
      ncid = 0;
    }
  }

  void MeshNetCDFHelper::GetDimensions(std::vector<size_t>& dimensions, const char* name) const throw (MeshReaderException)
  {
    // ID of variable
    int varid(0);
    nc_inq_varid(ncid, name, &varid);

    // Number of dimensions the variable has
    int ndims(0);
    nc_inq_varndims(ncid, varid, &ndims);

    // ID for each dimension
    std::vector<int> dimension_ids(ndims);
    nc_inq_var(ncid, varid, NULL, NULL, NULL, &dimension_ids[0], NULL);

    // Count of each dimension
    for (std::vector<int>::const_iterator dimension_id_iterator(dimension_ids.begin()); dimension_id_iterator != dimension_ids.end(); dimension_id_iterator++)
    {
      size_t length(0);
      nc_inq_dimlen(ncid, *dimension_id_iterator, &length);
      dimensions.push_back(length);
    }
  }

  template<class DataT> void MeshNetCDFHelper::ReadData(std::vector<DataT>& contents, const char* name) const throw (MeshReaderException)
  {
    int varid(0);
    nc_inq_varid(ncid, name, &varid);
    nc_get_var(ncid, varid, &contents[0]);
  }

  MeshReader::MeshReader()
  {
  }

  MeshReader::~MeshReader()
  {
  }

  void MeshReader::Read(const char* filename) throw(MeshReaderException)
  {
    // Open it
    MeshNetCDFHelper netcdf(filename);

    // Get dimensions of x and y
    std::vector<size_t> dimensions_x;
    std::vector<size_t> dimensions_y;
    netcdf.GetDimensions(dimensions_x, "Mesh2_node_x");
    netcdf.GetDimensions(dimensions_y, "Mesh2_node_y");

    // Should both have one dimension
    if (dimensions_x.size() != 1 || dimensions_y.size() != 1)
      throw MeshReaderException();

    // Number of points
    size_t npoints = dimensions_x[0];

    // Number of triangles
    std::vector<size_t> dimensions_triangles;
    netcdf.GetDimensions(dimensions_triangles, "Mesh2_face_nodes");
    size_t ntriangles = dimensions_triangles[0];

    // Make output buffers
    points.resize(3*npoints);
    triangles.resize(3*ntriangles);
   
    // Get latitude and longitude
    std::vector<double> latitude(npoints);
    std::vector<double> longitude(npoints);
    netcdf.ReadData<double>(longitude, "Mesh2_node_x");
    netcdf.ReadData<double>(latitude, "Mesh2_node_y");

    // Get vertices
    netcdf.ReadData<int32_t>(triangles, "Mesh2_face_nodes");

    // Convert (latitude,longitude) to unit vectors in cartesian coordinates
    std::vector<double>::const_iterator latitude_iterator(latitude.begin());
    std::vector<double>::const_iterator longitude_iterator(longitude.begin());
    std::vector<double>::iterator point_iterator(points.begin());
    while (point_iterator != points.end())
    {
      double latitude_radians = *latitude_iterator * M_PI / 180.0;
      double longitude_radians = *longitude_iterator * M_PI / 180.0;
      double z = sin(latitude_radians);
      double minor_r = cos(latitude_radians);
      double x = cos(longitude_radians) * minor_r;
      double y = sin(longitude_radians) * minor_r;
      *(point_iterator++) = x;
      *(point_iterator++) = y;
      *(point_iterator++) = z;
      latitude_iterator++;
      longitude_iterator++;
    }
  }

}
