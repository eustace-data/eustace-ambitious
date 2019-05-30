 // Reading of mesh files

#ifndef _EUSTACE_ANALYSIS_FILEIO_MESH_H_
#define _EUSTACE_ANALYSIS_FILEIO_MESH_H_

#include <stdint.h>
#include <math.h>
#include <vector>

namespace EUSTACE
{
  class MeshReaderException
  {
  };

  class MeshReader
  {
  public:
    MeshReader();

    virtual ~MeshReader();

  public:

    void Read(const char* pathname) throw(MeshReaderException);

  public:

    const std::vector<double>& Points() const
    { return points; }

    const std::vector<int32_t>& Triangles() const
    { return triangles; }

  private:    
    std::vector<double> points;
    std::vector<int32_t> triangles;
  };
}

#endif // #ifdef _EUSTACE_ANALYSIS_FILEIO_MESH_H_
