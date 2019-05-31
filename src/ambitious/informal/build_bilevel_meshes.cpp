#include "ambitious/debuglog/debuglog.hpp"
#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(build_bilevel_meshes_cpp, "$Revision: 1278 $")

// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
 
#include <eustace/analysis/fileio/mesh.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "ambitious/timer/timer.hpp"
#include "ambitious/to_string/to_string.hpp"
#include "ambitious/meshes/meshes.hpp"

#include <mpi.h>






typedef std::vector<fmesh::Mesh > MeshVector;
typedef std::vector<fmesh::TriangleLocator* > LocatorVector;







int main(int argc, char* argv[])
{
  VersionRegistry::log();
  
  //   MPI_Init(&argc, &argv);
  int rank = 0;
  //  int sz = 1;
  //  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //  MPI_Comm_size(MPI_COMM_WORLD, &sz);
  //  LOG_("MPI size = " << sz << std::endl);
  //  LOG_("MPI rank = " << rank << std::endl);

 // Get filename from command line
  if (argc != 4) {
    std::cerr << "Usage: example_mesh_load fileroot level0 level1" << std::endl;
    VERSIONS_PURGE();
    return 1;
  }
  const char* fileroot = argv[1];
  std::vector<int64_t> levels(2);
  levels[0] = atoi(argv[2]);
  levels[1] = atoi(argv[3]);

  // Set verbosity to 1.
  TimerHierarchy timer(2);

  std::cout << "*************************** gobal_mesh *************************" << std::endl;
  MeshBilevel global_mesh(&timer, false);
  global_mesh.set_levels(levels[0], levels[1]).import_python(fileroot);
  LOG_("global_mesh.macro_boundary_.rows = " << global_mesh.macro_boundary_.rows() << std::endl);

  std::cout << "*************************** submesh *************************" << std::endl;
  MeshBilevel submesh(&timer, false);
  std::set<int64_t> macro_vertices;
  macro_vertices.clear();
  macro_vertices.insert(4);
  submesh.define_subset_from_vertices(global_mesh, macro_vertices, false);
  LOG_("submesh.macro_boundary_.rows = " << submesh.macro_boundary_.rows() << std::endl);

  std::cout << "*************************** global_mesh.macro_mesh *************************" << std::endl;
  std::cout << global_mesh.macro_mesh() << std::endl;
  std::cout << "*************************** global_mesh.micro_mesh *************************" << std::endl;
  std::cout << global_mesh.micro_mesh() << std::endl;
  std::cout << "*************************** submesh.macro_mesh *************************" << std::endl;
  std::cout << submesh.macro_mesh() << std::endl;
  std::cout << "*************************** submesh.micro_mesh *************************" << std::endl;
  std::cout << submesh.micro_mesh() << std::endl;

  //  std::cout << "*************************** submesh.micro_mesh.S *************************" << std::endl;
  //  std::cout << submesh.micro_mesh().S() << std::endl;
  //  std::cout << "*************************** submesh.micro_mesh.TV *************************" << std::endl;
  //  std::cout << submesh.micro_mesh().TV() << std::endl;


  //  LOG_("Python\n"
  //       << "submesh.vertex_index\n"
  //       << submesh.supermesh_index_->vertex()
  //       << "submesh.triangle_index0\n"
  //       << submesh.supermesh_index_->triangle_macro()
  //       << "submesh.triangle_index1\n"
  //       << submesh.supermesh_index_->triangle_micro());
  //  LOG_("Supermesh\n"
  //       << "submesh.vertex_index\n"
  //       << submesh.supermesh_index_->vertex()
  //       << "submesh.triangle_index0\n"
  //       << submesh.supermesh_index_->triangle_macro()
  //       << "submesh.triangle_index1\n"
  //       << submesh.supermesh_index_->triangle_micro());
  

  
  //  Eigen::Matrix<double, 1, 2> points(90.0, 0.0);
  Eigen::Matrix<double, Eigen::Dynamic, 3> points;
  points.resize(submesh.micro_mesh().nV(), 3);
  for (int64_t v=0; v < (int64_t)submesh.micro_mesh().nV(); ++v) {
    points.row(v) <<
      submesh.micro_mesh().S(v)[0],
      submesh.micro_mesh().S(v)[1],
      submesh.micro_mesh().S(v)[2];      
  }
  
  Eigen::Matrix<int, Eigen::Dynamic, 1> point2T;
  Eigen::Matrix<double, Eigen::Dynamic, 3> point2bary;
  Eigen::Matrix<int, Eigen::Dynamic, 1> T0;
  Eigen::Matrix<double, Eigen::Dynamic, 3> bary0;

  timer.tic("Locate points");
  MeshBilevel* themesh;
  for (int64_t switcher=0; switcher < 2; ++switcher) {
    if (switcher == 0) {
      timer.tic("Global mesh");
      std::cout << "Global mesh" << std::endl;
      themesh = &global_mesh;
    } else {
      timer.tic("Sub-mesh");
      std::cout << "Sub-mesh" << std::endl;
      themesh = &submesh;
    }
    for (int64_t level=0; level<2; ++level) {
      const fmesh::Mesh* mesh;
      if (level == 0) {
	timer.tic("Macro level");
	std::cout << "Macro level" << std::endl;
	timer.tic("macro_bary");
	themesh->macro_bary(points, point2T, point2bary);
	timer.toc();
	mesh = &themesh->macro_mesh();
      } else {
	timer.tic("Micro level");
	std::cout << "Micro level" << std::endl;
	timer.tic("micro_bary");
	themesh->micro_bary(points, point2T, point2bary, T0, bary0);
	timer.toc();
	mesh = &themesh->micro_mesh();
      }
      int64_t not_found = 0;
      for (int64_t v=0; v < point2T.rows(); ++v) {
	int64_t t = point2T(v,0);
	not_found += (t < 0);
      }
      std::cout << "not_found = " << not_found << std::endl;
      if (false) {
      if (point2T(0,0) < 0) {
	std::cout << "Point 0 not found" << std::endl;
      } else {
	int64_t t = point2T(0,0);
	std::cout << "TV(t) = "
		  << mesh->TV(t)[0] << ", "
		  << mesh->TV(t)[1] << ", "
		  << mesh->TV(t)[2]
		  << std::endl;
	std::cout << "S(TV(t)) = \t"
		  << mesh->S(mesh->TV(t)[0]) << std::endl << "\t\t"
		  << mesh->S(mesh->TV(t)[1]) << std::endl << "\t\t"
		  << mesh->S(mesh->TV(t)[2])
		  << std::endl;
	std::cout << "Bary = \t"
		  << point2bary(0,0) << ", "
		  << point2bary(0,1) << ", "
		  << point2bary(0,2) << ", "
		  << std::endl;
	if (level == 1) {
	  mesh = &themesh->macro_mesh();
	  t = T0(0,0);
	  std::cout << "t0 = " << t << std::endl;
	  std::cout << "TV0(t0) = "
		    << mesh->TV(t)[0] << ", "
		    << mesh->TV(t)[1] << ", "
		    << mesh->TV(t)[2]
		    << std::endl;
	  std::cout << "S0(TV0(t0)) = \t"
		    << mesh->S(mesh->TV(t)[0]) << std::endl << "\t\t"
		    << mesh->S(mesh->TV(t)[1]) << std::endl << "\t\t"
		    << mesh->S(mesh->TV(t)[2])
		    << std::endl;
	  std::cout << "Bary0 = \t"
		    << bary0(0,0) << ", "
		    << bary0(0,1) << ", "
		    << bary0(0,2) << ", "
		    << std::endl;
	}
      }
      }
      timer.toc();
    }
    timer.toc();
  }
  timer.toc();
  timer.toc();

  /*
  for (int64_t v=0; v < submesh.micro_mesh().nV(); ++v) {
    std::cout << "W(v=" << v << ") =\t"
	      << submesh.weights()[v] << ",\t"
	      << submesh.macro_hat()[v] << std::endl;
  }
  */
  
  Eigen::Matrix<double, Eigen::Dynamic, 1> DC0;
  Eigen::Matrix<double, Eigen::Dynamic, 1> DC0inv;
  Eigen::SparseMatrix<double> C1, G1, B1;
  Eigen::Matrix<double, Eigen::Dynamic, 1> Tareas;
  timer.tic("Calculate global mesh Q-blocks.");
  calcQblocks(global_mesh.micro_mesh(), DC0, DC0inv, C1, G1, B1, Tareas);
  timer.toc();
  timer.tic("Calculate submesh Q-blocks.");
  calcQblocks(submesh.micro_mesh(), DC0, DC0inv, C1, G1, B1, Tareas);
  timer.toc();





  //  if (rank == 0) {
    global_mesh.export_direct("global_mesh");
    submesh.export_direct("submesh_mesh");
    //  }

  
    if (rank == 0) {
      timer.print();
    }

    //  LOG_("MPI barrier rank = " << rank << std::endl);
    //  MPI_Barrier(MPI_COMM_WORLD);

    //  LOG_("MPI pre-finalize rank = " << rank << std::endl);
    //  MPI_Finalize();
    //  LOG_("MPI post-finalize rank = " << rank << std::endl);

  VERSIONS_PURGE();
  return 0;
}
