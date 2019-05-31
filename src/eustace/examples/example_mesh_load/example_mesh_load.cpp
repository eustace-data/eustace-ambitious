/* Example to load EUSTACE mesh */

#include <eustace/analysis/fileio/mesh.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

int main(int argc, char* argv[])
{
  std::cerr << "example_mesh_load" << std::endl;

  // Get filename from command line
  if (argc != 2)
  {
    std::cerr << "Usage: example_mesh_load filename" << std::endl;
    exit(1);
  }
  const char* filename = argv[1];
  std::cerr << "Opening mesh file: " << filename << std::endl;


  // Read it
  EUSTACE::MeshReader reader;
  try
  {
    reader.Read(filename);
  }
  catch (const EUSTACE::MeshReaderException&)
  {
    std::cerr << "Error reading file";
    exit(1);
  }

  // Show data
  const std::vector<double>& points = reader.Points();
  const std::vector<int32_t>& triangles = reader.Triangles();
  size_t num_points = points.size() / 3;
  size_t num_triangles = triangles.size() / 3;
  for (size_t point_index = 0; point_index < num_points; point_index++)
  {
    std::cout << "point [" << std::setw(12) << point_index << "]: " << 
      std::setw(12) << points[3*point_index  ] << " " <<
      std::setw(12) << points[3*point_index+1] << " " <<
      std::setw(12) << points[3*point_index+2] << std::endl;
  }
  for (size_t triangle_index = 0; triangle_index < num_triangles; triangle_index++)
  {
    std::cout << "triangle [" << std::setw(12) << triangle_index << "]: " << 
      std::setw(12) << triangles[3*triangle_index  ] << " " <<
      std::setw(12) << triangles[3*triangle_index+1] << " " <<
      std::setw(12) << triangles[3*triangle_index+2] << std::endl;
  }

  return 0;
}
