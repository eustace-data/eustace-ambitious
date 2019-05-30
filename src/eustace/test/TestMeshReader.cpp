
#include <iostream>
#include <fstream>
#include <vector>
#include <memory.h>
#include <cppunit/extensions/HelperMacros.h>
#include <eustace/analysis/fileio/mesh.h>

using namespace EUSTACE;
using namespace std;

class TestMeshReader  : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( TestMeshReader );
  CPPUNIT_TEST( testReadIcosahedron );
  CPPUNIT_TEST_SUITE_END();

public:
 
  void testReadIcosahedron()
  {
    // Filename
    const char* filename = EUSTACE_TEST_DATA_DIRECTORY "/eustace_mesh_level0.nc";

    // Make reader and read data
    MeshReader reader;
    reader.Read(filename);
    
    // Check sizes
    CPPUNIT_ASSERT_EQUAL(size_t(3*12), reader.Points().size());
    CPPUNIT_ASSERT_EQUAL(size_t(3*20), reader.Triangles().size());    

    // Expected triangle contents
    int32_t expected_triangles[60] =
    {
       4,  7,  8,
       4,  9,  7,
       5, 11,  6,
       5,  6, 10,
       0,  3,  4,
       0,  5,  3,
       2,  1,  7,
       2,  6,  1,
       8, 11,  0,
       8,  1, 11,
       9,  3, 10,
       9, 10,  2,
       8,  0,  4,
      11,  5,  0,
       4,  3,  9,
       5, 10,  3,
       7,  1,  8,
       6, 11,  1,
       7,  9,  2,
       6,  2, 10
    };

    // Expected vertices contents
    // - but will be rotated to put face centre at pole
    double a = 0.850650808352;
    double b = 0.525731112119;
    double expected_points[36] = 
    {
      a, b, 0,
     -a, b, 0,
     -a,-b, 0,
      a,-b, 0,
      b, 0, a,
      b, 0,-a,
     -b, 0,-a,
     -b, 0, a,
      0, a, b,
      0,-a, b,
      0,-a,-b,
      0, a,-b
    };

    // Rotation needed to put poles at centres of faces
    double golden_ratio = (1.0 + sqrt(5.0)) * 0.5;
    double one_minus_inverse_ratio = 1.0 - (1.0 / golden_ratio);
    double s = 1.0 / sqrt(1.0 + one_minus_inverse_ratio*one_minus_inverse_ratio);
    double c = sqrt(1.0 - s*s);
    for (int i  = 0; i < 12; i++)
    {
      double x = expected_points[i*3  ];
      double z = expected_points[i*3+2];
      expected_points[i*3  ] = c*x - s*z;
      expected_points[i*3+2] = s*x + c*z;
    }
    
    // Check
    CPPUNIT_ASSERT_EQUAL(expected_triangles[ 0], reader.Triangles()[ 0]);
    CPPUNIT_ASSERT_EQUAL(expected_triangles[ 1], reader.Triangles()[ 1]);
    CPPUNIT_ASSERT_EQUAL(expected_triangles[ 2], reader.Triangles()[ 2]);
    CPPUNIT_ASSERT_EQUAL(expected_triangles[ 3], reader.Triangles()[ 3]);
    CPPUNIT_ASSERT_EQUAL(expected_triangles[ 4], reader.Triangles()[ 4]);
    CPPUNIT_ASSERT_EQUAL(expected_triangles[ 5], reader.Triangles()[ 5]);
    CPPUNIT_ASSERT_EQUAL(expected_triangles[57], reader.Triangles()[57]);
    CPPUNIT_ASSERT_EQUAL(expected_triangles[58], reader.Triangles()[58]);
    CPPUNIT_ASSERT_EQUAL(expected_triangles[59], reader.Triangles()[59]);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[ 0], reader.Points()[ 0], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[ 1], reader.Points()[ 1], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[ 2], reader.Points()[ 2], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[ 3], reader.Points()[ 3], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[ 4], reader.Points()[ 4], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[ 5], reader.Points()[ 5], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[ 6], reader.Points()[ 6], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[ 7], reader.Points()[ 7], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[ 8], reader.Points()[ 8], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[ 9], reader.Points()[ 9], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[10], reader.Points()[10], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[11], reader.Points()[11], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[12], reader.Points()[12], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[13], reader.Points()[13], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[14], reader.Points()[14], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[15], reader.Points()[15], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[16], reader.Points()[16], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[17], reader.Points()[17], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[18], reader.Points()[18], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[19], reader.Points()[19], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[20], reader.Points()[20], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[21], reader.Points()[21], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[22], reader.Points()[22], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[23], reader.Points()[23], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[24], reader.Points()[24], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[25], reader.Points()[25], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[26], reader.Points()[26], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[27], reader.Points()[27], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[28], reader.Points()[28], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[29], reader.Points()[29], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[30], reader.Points()[30], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[31], reader.Points()[31], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[32], reader.Points()[32], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[33], reader.Points()[33], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[34], reader.Points()[34], 1.0E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_points[35], reader.Points()[35], 1.0E-6);
    
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestMeshReader );
