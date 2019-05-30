#ifndef MESHES_HPP
#define MESHES_HPP

#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(meshes_hpp, "$Revision: 1297 $")

#include "ambitious/debuglog/debuglog.hpp"

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <map>
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <mpi.h>

#include <eustace/analysis/fileio/mesh.h>
#include <eustace/timeutils/timebase.h>
#include <fmesher.hh>
#include "ambitious/ostreamer/ostreamer.hpp"
#include "ambitious/timer/timer.hpp"
#include "ambitious/bidirmap/bidirmap.hpp"
#include "ambitious/to_string/to_string.hpp"
#include "ambitious/meshes/span.hpp"





enum class TimeUnits { Day, Year, Season };
typedef Span<int64_t> DaySpan;


// Forward declaration
class Converter {
 public:
  static void euclidean_to_latlong(const fmesh::Point& loc, std::pair<double, double>& latlong);
  static void latlong_to_euclidean(const std::pair<double, double>& latlong, fmesh::Point& loc);
};





/** Implements the FEM constructions of Lindgren et al (2011), appendix A.2
 * 
 * Constructs the precision building blocks \f$C_0\f$, \f$C_0^{-1}\f$, \f$C_1\f$,
 * \f$G_1\f$, \f$B_1\f$, and the triangle areas (Tareas).
 *
 * The precision of \f$(\kappa^2 - \Delta) x(s) = w(s)\f$ with a space-independent
 * spatial (inverse) range parameter \f$\kappa\f$ is
 * \f$Q(\kappa) = \kappa^4 C_0 + 2\kappa^2 G_1 + G_2 \f$ where \f$G_2=G_1 C_0^{-1} G_1\f$
 */

template <class DerivedA, class DerivedB>
void calcQblocks(fmesh::Mesh M,
		 Eigen::MatrixBase<DerivedA> const & DC0_,
		 Eigen::MatrixBase<DerivedA> const & DC0inv_,
		 Eigen::SparseMatrixBase<DerivedB> const & C1_,
		 Eigen::SparseMatrixBase<DerivedB> const & G1_,
		 Eigen::SparseMatrixBase<DerivedB> const & B1_,
		 Eigen::MatrixBase<DerivedA> const & Tareas_)
{
  Eigen::MatrixBase<DerivedA>& DC0 = const_cast< Eigen::MatrixBase<DerivedA>& >(DC0_);
  Eigen::MatrixBase<DerivedA>& DC0inv = const_cast< Eigen::MatrixBase<DerivedA>& >(DC0inv_);
  Eigen::SparseMatrixBase<DerivedB>& C1 = const_cast< Eigen::SparseMatrixBase<DerivedB>& >(C1_);
  Eigen::SparseMatrixBase<DerivedB>& G1 = const_cast< Eigen::SparseMatrixBase<DerivedB>& >(G1_);
  Eigen::SparseMatrixBase<DerivedB>& B1 = const_cast< Eigen::SparseMatrixBase<DerivedB>& >(B1_);
  Eigen::MatrixBase<DerivedA>& Tareas = const_cast< Eigen::MatrixBase<DerivedA>& >(Tareas_);

  DC0.derived().resize(M.nV(), 1);
  DC0.derived().setZero();
  DC0inv.derived().resize(M.nV(), 1);
  C1.derived().resize(M.nV(), M.nV());
  G1.derived().resize(M.nV(), M.nV());
  B1.derived().resize(M.nV(), M.nV());
  Tareas.derived().resize(M.nT(), 1);

  uint64_t nnz1 = M.nV() * 6;
    
  uint64_t n1_triplet = M.nT() * 9;
  std::vector< Eigen::Triplet<double> > C1_triplet;
  std::vector< Eigen::Triplet<double> > G1_triplet;
  std::vector< Eigen::Triplet<double> > B1_triplet;
  C1_triplet.reserve(n1_triplet);
  G1_triplet.reserve(n1_triplet);
  B1_triplet.reserve(n1_triplet);

  fmesh::Point e[3];
  for (int64_t t = 0; t < (int64_t)M.nT(); t++) {
    const fmesh::Int3Raw& tv = M.TV()[t].raw();
    const fmesh::Point& s0 = M.S()[tv[0]];
    const fmesh::Point& s1 = M.S()[tv[1]];
    const fmesh::Point& s2 = M.S()[tv[2]];
    e[0].diff(s2,s1);
    e[1].diff(s0,s2);
    e[2].diff(s1,s0);
    
    fmesh::PointRaw eij[3];
    for (int i=0; i<3; i++) {
      eij[i][i] = fmesh::Vec::scalar(e[i],e[i]);
      for (int j=i+1; j<3; j++) {
	eij[i][j] = fmesh::Vec::scalar(e[i],e[j]);
	eij[j][i] = eij[i][j];
      }
    }
    
    bool b[3];
    b[0] = (M.TT()[t][0] < 0 ? true : false);
    b[1] = (M.TT()[t][1] < 0 ? true : false);
    b[2] = (M.TT()[t][2] < 0 ? true : false);
    
    double a = M.triangleArea(t);
    Tareas(t,0) = a;
    
    /* "Flat area" better approximation for use in G-calculation. */
    double fa = fmesh::Point().cross(e[0],e[1]).length()/2.0;
    
    double vij;
    for (int i=0; i<3; i++) {
      DC0(tv[i],0) += a/3.0;
      C1_triplet.push_back(Eigen::Triplet<double>(tv[i], tv[i], a/6.0));
      G1_triplet.push_back(Eigen::Triplet<double>(tv[i], tv[i], eij[i][i]/(4.*fa)));
      for (int j=i+1; j<3; j++) {
	C1_triplet.push_back(Eigen::Triplet<double>(tv[i], tv[j], a/12.0));
	C1_triplet.push_back(Eigen::Triplet<double>(tv[j], tv[i], a/12.0));
	vij = eij[i][j]/(4.*fa);
	G1_triplet.push_back(Eigen::Triplet<double>(tv[i], tv[j], vij));
	G1_triplet.push_back(Eigen::Triplet<double>(tv[j], tv[i], vij));
      }
    }
    
    if (b[0] || b[1] || b[2]) {
      vij = -1./(4.*fa);
      for (int i=0; i<3; i++) {
	for (int j=0; j<3; j++) {
	  for (int k=0; k<3; k++) {
	    if (b[k] && (i != k)) {
	      B1_triplet.push_back(Eigen::Triplet<double>(tv[i], tv[j], eij[k][j]*vij));
	    }
	  }
	}
      }
    }
  }
  
  G1.derived().reserve(nnz1);
  B1.derived().reserve(nnz1);
  C1.derived().setFromTriplets(C1_triplet.begin(), C1_triplet.end());
  G1.derived().setFromTriplets(G1_triplet.begin(), G1_triplet.end());
  B1.derived().setFromTriplets(B1_triplet.begin(), B1_triplet.end());
  
  DC0inv.derived() = 1.0 / DC0.array();
}







/**
 * If target_level < micro_level, returns the macro-triangle
 * containing the micro triangle.
 *
 * Subdivision level \f$A < B\f$, \f$T_B\f$ is the level B triangle index.
 * Then the macro-triangle \f$T_A\f$ containing \f$T_B\f$ is given by
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *   T_A = (T_B >> ((B-A) << 1))
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * If target_level > micro_level, returns the index of the 0th subdivision micro-triangle.
 */
int64_t get_macro_triangle_index(int64_t micro_index,
				 int64_t micro_level,
				 int64_t macro_level);


/** Get the micro-triangle index for one of the corners (or middle, for corner == 3) of a
 * macro triangle */
int64_t get_micro_triangle_index(int64_t macro_index,
				 int64_t macro_level,
				 int64_t micro_level);
int64_t get_micro_triangle_index(int64_t macro_index,
				 int64_t macro_level,
				 int64_t micro_level,
				 int64_t corner);

/**
 * Which micro-triangle is it, with respect to the macro-triangle?
 * The remainder of division by 4 at the level above the target, if
 * recursive==false.
 *
 * If recursive==true, the remainder of division by 4^(B-A)-1.
 */
int64_t get_micro_in_macro_triangle_index(int64_t micro_index,
					  int64_t micro_level,
					  int64_t macro_level,
					  bool recursive = false);


/** Allocate and initialise a fmesh::TriangleLocator object for a Mesh
 */
fmesh::TriangleLocator* initialise_locator(const fmesh::Mesh& M);

/** Helper function for allocating and initialising locators for a vector of meshes
 */
void initialise_locators_helper(std::vector<fmesh::TriangleLocator* >& locators,
				std::vector<fmesh::Mesh >& meshes,
				TimerHierarchy* timer);


/** Map between CW/CCW triangle ordering in the subdivision scheme.
 *
 * For each triangle subdivision, the local indices are mapped according to
 * `0 -> 0`, `1 -> 2`, `2 -> 1`, and `3 -> 3`.
 * This reordering is done within each of the 20 level-0 macro-triangles.
 */
int64_t triangle_index_swapwise(int64_t triangle_index, int64_t level);



/** Build a subdivided icosahedron from stored files generated by a Python script
 *
 * The Python version stores the triangle-vertex information and sub-triangles clockwise:
 *
 *                1
 *                /\
 *               /  \
 *             e0    e1
 *             /  t1  \
 *            /        \
 *          “2”---e2---“0”
 *           /\         /\
 *          /  \  t3   /  \
 *        e0    e1   e0    e1
 *        /  t0  \   /  t2  \
 *       /        \ /        \
 *     0/---------“1”---------\2
 *
 * The C++ code stores counterclockwise,
 * which would be the following:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~
 *     Macro-triangle: 0 2 1
 *     t0:  0  "1" "2"
 *     t1: "2" "0"  1
 *     t2: "1"  2  "0"
 *     t3: "0" "2" "1"
 * ~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Renumbering the local node references,
 * as well as the triangle orders:
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~
 *                2
 *                /\
 *               /  \
 *              /    \
 *             /  t2  \
 *            /        \
 *          “1”--------“0”
 *           /\         /\
 *          /  \  t3   /  \
 *         /    \     /    \
 *        /  t0  \   /  t1  \
 *       /        \ /        \
 *     0/---------“2”---------\1
 * ~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * The C++ counterclockwise storage becomes:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~
 * Macro-triangle: 0 1 2
 * t0:  0  "2" "1"
 * t1: "2"  1  "0"
 * t2: "1" "0"  2
 * t3: "0" "1" "2"
 * ~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Subdivision level A < B, T_B is the level B triangle index.
 * Then the macro-triangle T_A containing T_B is given by
 * ~~~~~~~~~~~~~~~~~~~~~~~~~
 *   T_A = (T_B >> ((B-A) << 1))
 * ~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~
 * T_A == 0 is subdivided into T_{A+1} = 0,1,2,3
 * T_A == 1 is subdivided into T_{A+1} = 4,5,6,7
 * etc
 * ~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * It's not clear what the global order of the vertices are, since
 * that emerges from the edge construction in the Python code, but
 * they do at least fulfil V_0 < V_1 < V_2 < ..., where V_. are the
 * extra vertices added at each subdivision level.  The order of the
 * vertices is kept intact; only the triangle and per-triangle
 * orderings are changed to CCW order in the C++ code.
 *
 */
void initialise_mesh_from_python(fmesh::Mesh& mesh, int64_t level, const std::string& fileroot);

/** Initialise multiple meshes from Python-generated files.
 */
void initialise_meshes_from_python(std::vector<fmesh::Mesh >& meshes,
				   const std::vector<int64_t>& levels,
				   const std::string& fileroot,
				   TimerHierarchy* timer);






/** Calculate a modified halo hat functions that decreases to zero at a specified halo depth.
 *
 * - The halo hat functions for root nodes in a macro triangle add up to one.
 * - The halo hat functions for non-root nodes are defined to be zero.
 * - The halo hat functions are zero for triangles not connected to a root node.
 *
 * Necessary and sufficient conditions for valid root node sets:
 * -# Necessary: Interior nodes must be root nodes.
 * -# Necessary: Boundary nodes must be root nodes if they have more than one non-boundary edge connection.
 * -# Necessary: At least one node of each interior edge must be a root node. (Boundary-to-boundary connections.)
 * -# Necessary: At least one node of each triangle must be a root node. (Isolated triangles.)
 * -# Together, the four conditions above are Sufficient.
 * -# Conditions 1 and 2 can be combined: Any node must be a root node if they are attached to more than 2 triangles.
 *
 * @param qbary_ A row vector of triangle barycentric coordinates, with element 0 denoting the root node.
 * @param root_nodes_ A row vector of booleans, with true denoting that a macro node is a root node.
 *   Halo hat functions for non-root nodes are zero. 
 * @param boundary_edges_ A row vector of booleans, with true denoting that a opposing macro node is on the boundary.
 * @param halo_depth The relative halo depth. the minimum absolute halo depth is 0
 *   (meaning that only the opposing edge is in the halo).
 *   A macro triangle with 8 subdivision intervals along each edge,
 *   and absolute halo depth 3, has a relative halo depth of \f$3/8\f$.
 */
template <typename DerivedA, typename DerivedB, typename DerivedC, typename DerivedD>
void calc_macro_halo_hat(const Eigen::DenseBase<DerivedA> & qbary,
                         const Eigen::DenseBase<DerivedB> & root_nodes,
                         const Eigen::DenseBase<DerivedC> & boundary_edges,
                         double halo_depth,
                         Eigen::DenseBase<DerivedD> const & hh_)
{
  Eigen::DenseBase<DerivedD>& hh = const_cast< Eigen::DenseBase<DerivedD>& >(hh_);

  LOG("" << std::endl);
  hh.derived().resize(1,3);
  LOG("" << std::endl);
  hh.derived().setZero();
  LOG("" << std::endl);
  if (!root_nodes(0,0) && !root_nodes(0,1) && !root_nodes(0,2)) {
    return;
  }
  LOG("" << std::endl);
  for (int i = 0; i < 3; ++i) {
  LOG("" << std::endl);
    if (!root_nodes(0, (i+1)%3) && !root_nodes(0,(i+2)%3)) {
  LOG("" << std::endl);
      hh(0,i) = 1.0;
      return;
    }
  }
  LOG("" << std::endl);
  // At least two nodes are now root nodes.

  double hhsum = 0.0;
  for (int i = 0; i < 3; ++i) {
  LOG("" << std::endl);
    if (root_nodes(0,i)) {
  LOG("" << std::endl);
      hh(0, i) = std::max(0.0, qbary(0, i) - halo_depth);
      hhsum += hh(0,i);
    }
  }
  LOG("" << std::endl);
  hh.derived() /= hhsum;
  LOG("" << std::endl);

  if (root_nodes(0,0) && root_nodes(0,1) && root_nodes(0,2)) {
  LOG("" << std::endl);
    return;
  }
  LOG("" << std::endl);
  // Precisely one node is now a non-root node.

  // Non-root node with two adjacent boundary edges: share
  // Non-root node with one adjacent boundary edges: opposing root node takes all
  // Non-root node with no boundary edges in the triangle: not allowed
  for (int i = 0; i < 3; ++i) {
  LOG("" << std::endl);
    if (!root_nodes(0,i)) {
  LOG("" << std::endl);
      if (boundary_edges(0,(i+1)%3)) {
  LOG("" << std::endl);
        hh(0,(i+1)%3) += hh(0,i);
      }
  LOG("" << std::endl);
      if (boundary_edges(0,(i+2)%3)) {
  LOG("" << std::endl);
        hh(0,(i+2)%3) += hh(0,i);
      }
  LOG("" << std::endl);
      hh(0,i) = 0.0;
  LOG("" << std::endl);
    }
  LOG("" << std::endl);
  }
  LOG("" << std::endl);
  // Normalise
  hhsum = hh(0,0) + hh(0,1) + hh(0,2);
  LOG("" << std::endl);
  hh.derived() /= hhsum;
  LOG("" << std::endl);
}



/** Calculate the macro-triangle weight function for given quasi-Barycentric coordinates.
 *
 * The square of the weight functions for the three corners of a macro triangle add up to one.
 *
 * To handle (super)mesh boundaries, the input must be computed with calc_macro_halo_hat.
 *
 * 1   1 1/2 0 0
 * 1   1   0 0
 * 1/2 0   0
 * 0   0
 * 0
 *
 * 1 1
 * 1
 *
 * 1 2 1
 * 2 2
 * 1
 *
 * 1 3 3 1
 * 3 6 3
 * 3 3
 * 1
 *
 * 1  4  6 4 1
 * 4 12 12 4
 * 6 12 6
 * 4  4
 * 1
 *
 * \f{align*}{b0^4 + 4 b0^3 (b1+b2) + b0^2 (6 b1^2 /2 + 12 b1 b2 + 6 b2^2 /2)
 *  \\&= b0^2 (b0^2 + 4 b0 (b1+b2) + (3 b1^2 + 12 b1 b2 + 3 b2^2)
 * \\&= b0^2 (b0^2 + b1 (4 b0 + 3 b1 + 6 b2) + b2 (4 b0 + 6 b1 + 3 b2)
 * \\&= b0^2 (b0 (b0 + b1 + b2)  + b1 (3 b0 + 3 b1 + 6 b2) + b2 (3 b0 + 6 b1 + 3 b2)
 * \\&= b0^2 (b0 + b1 (3 b0 + 3 b1 + 3 b2 + 3 b2) + b2 (3 b0 + 3 b1 + 3 b1 + 3 b2)
 * \\&= b0^2 (b0 + 3 b1 (1 + b2) + 3 b2 (1 + b1))
 * \\&= b0^2 (b0 + 3 b1 + 3 b2 + 6 b1 b2)
 * \\&= b0^2 (1 + 2 b1 + 2 b2 + 6 b1 b2)
 * \f}
 */
template <typename DerivedA>
double calc_macro_weights(const Eigen::DenseBase<DerivedA> & qbary, int64_t corner)
{
  return qbary(0,corner) * std::sqrt(1 + 2*qbary(0,(corner+1)%3) + 2*qbary(0,(corner+2)%3) + 6*qbary(0,(corner+1)%3)*qbary(0,(corner+2)%3));
}




/**
 * Quasi-Barycentric coordinates of sub-triangle nodes with respect to a
 * macro triangle one level above:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * macro_bary(0,...,5):
 *    (1,0,0)
 *    (0,1,0)
 *    (0,0,1)
 *    (0,1/2,1/2)
 *    (1/2,0,1/2)
 *    (1/2,1/2,0)
 * macro_tri: (0,1,2)
 * sub_tri(0,1,2,3):
 *    (0,5,4)
 *    (5,1,3)
 *    (4,3,2)
 *    (3,4,5)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Can operate in-place
 */
template <typename DerivedB, typename DerivedC, typename DerivedD, typename DerivedE>
void map_tbary(int64_t sub_level,
	       int64_t target_level,
	       const Eigen::DenseBase<DerivedB> & T,
	       const Eigen::DenseBase<DerivedC> & bary,
	       Eigen::DenseBase<DerivedD> const & T_target_,
	       Eigen::DenseBase<DerivedE> const & bary_target_)
{
  Eigen::DenseBase<DerivedD>& T_target = const_cast< Eigen::DenseBase<DerivedD>& >(T_target_);
  Eigen::DenseBase<DerivedE>& bary_target = const_cast< Eigen::DenseBase<DerivedE>& >(bary_target_);
  assert(target_level <= sub_level);
  assert(T.rows() == bary.rows());

  T_target.derived().resize(T.rows(), T.cols());
  bary_target.derived().resize(bary.rows(), bary.cols());
  
  Eigen::Matrix<double, 6, 3> macro_bary;
  macro_bary << 1,0,0, 0,1,0, 0,0,1, 0,0.5,0.5, 0.5,0,0.5, 0.5,0.5,0;
  Eigen::Matrix<int64_t, 4, 3> sub_triangles;
  sub_triangles << 0,5,4, 5,1,3, 4,3,2, 3,4,5;

  Eigen::Matrix<double, 1, 3> target_bary;
  for (Eigen::Index point_index=0; point_index < T.rows(); ++point_index) {
    int64_t target_triangle = T(point_index);
    if (target_triangle < 0) { // point not found / doesn't exist
      T_target(point_index) = -1;
      bary_target.row(point_index) << 0.0, 0.0, 0.0;
    } else {
      target_bary = bary.row(point_index);
      for (int64_t level=sub_level; level > target_level; --level) {
	int64_t sub_triangle_index = get_micro_in_macro_triangle_index(target_triangle,
								       level,
								       level-1,
								       false);
	target_triangle = get_macro_triangle_index(target_triangle, level, level-1);
	target_bary =
	  macro_bary.row(sub_triangles(sub_triangle_index, 0)) * target_bary(0) +
	  macro_bary.row(sub_triangles(sub_triangle_index, 1)) * target_bary(1) +
	  macro_bary.row(sub_triangles(sub_triangle_index, 2)) * target_bary(2);
      }
      T_target(point_index) = target_triangle;
      bary_target.row(point_index) = target_bary;
    }
  }
}




template <typename DerivedA, typename DerivedB, typename DerivedC>
void map_points_to_mesh(const fmesh::Mesh& M,
			const fmesh::TriangleLocator* locator,
			const Eigen::DenseBase<DerivedA>& points,
			Eigen::DenseBase<DerivedB> const & T_,
			Eigen::DenseBase<DerivedC> const & bary_)
{
  Eigen::DenseBase<DerivedB>& T = const_cast< Eigen::DenseBase<DerivedB>& >(T_);
  Eigen::DenseBase<DerivedC>& bary = const_cast< Eigen::DenseBase<DerivedC>& >(bary_);
  T.derived().resize(points.rows(), 1);
  bary.derived().resize(points.rows(), 3);

  int t;
  fmesh::Point s;
  fmesh::Point b;

  for (int i=0; i < points.rows(); ++i) {
    if (points.cols() == 3) {
      s[0] = points(i,0);
      s[1] = points(i,1);
      s[2] = points(i,2);
    } else { // Convert from lat/lon
      Converter::latlong_to_euclidean(std::pair<double, double>(points(i,0), points(i,1)), s);
    }
    if (locator != NULL) {
      LOG("Use locator\n");
      t = locator->locate(s);
    } else {
      fmesh::Dart dh = fmesh::Dart();
      dh = M.locate_point(dh, s);
      if (dh.isnull()) {
	t = -1;
      } else {
	t = dh.t();
      }
    }
    if (t>=0) { /* Point located. */
      M.barycentric(fmesh::Dart(M,t),s,b); /* Coordinates relative to
					      canonical vertex
					      ordering. */
      T(i,0) = t;
      bary.row(i) << b[0], b[1], b[2];
    } else { /* Point not found. */
      T(i,0) = -1;
      bary.row(i) << 0.0, 0.0, 0.0;
    }
  }
}





/**
 * Input: mesh, lat&long, time_mesh, time, is_cell (T/F)
 * Output: Evaluation/integration weights.
 * 
 */






class MeshBilevelIndexMap {
 public:
  typedef BidirectionalMap<int64_t,int64_t> IndexMap;
 protected:
  int64_t macro_basis_;
  int64_t micro_basis_;
  /** Map between originator indices and actual indices.
   * AtoB: Indexed by Python or supermesh indices to return corresponding C mesh/submesh indices.
   * BtoA: Indexed by C mesh/submesh indices to return corresponding Python or supermesh indices.
  */
  IndexMap vertex_index_;
  IndexMap triangle_index0_;
  IndexMap triangle_index1_;
 public:
  MeshBilevelIndexMap()
      : macro_basis_(0),
        micro_basis_(0)
  {
    clear();
  }
  void clear() {
    vertex_index_.clear();
    triangle_index0_.clear();
    triangle_index1_.clear();
  }
  void num_basis(int64_t macro_basis, int64_t micro_basis) {
    macro_basis_ = macro_basis;
    micro_basis_ = micro_basis;
  }
  int64_t macro_num_basis() const {
    return macro_basis_;
  }
  int64_t micro_num_basis() const {
    return macro_basis_;
  }
  /** Vertex IndexMap */
  IndexMap& vertex() {
    return vertex_index_;
  }
  /** Triangle IndexMap for a macro mesh */
  IndexMap& triangle_macro() {
    return triangle_index0_;
  }
  /** Triangle IndexMap for a micro mesh */
  IndexMap& triangle_micro() {
    return triangle_index1_;
  }

  void log_dump() {
    PLOG_("MeshBilevelIndexMap" << std::endl);
    PLOG_("macro_basis_ = " << macro_basis_ << std::endl);
    PLOG_("micro_basis_ = " << micro_basis_ << std::endl);
    PLOG_("vertex_index_ = " << std::endl << vertex_index_);
    PLOG_("triangle_index0_ = " << std::endl << triangle_index0_);
    PLOG_("triangle_index1_ = " << std::endl << triangle_index1_);
  }
};















class MeshBilevel; // Forward declaration.
class SpaceTimeMapper; // Forward declaration.





/** Iterator template for the individual elements of a vector of vectors
 *
 * Use this to iterate over all the entries.
 */
template <class Value>
class dualiterator {
 public:
  typedef std::forward_iterator_tag iterator_category;
  typedef Value value_type;
  typedef std::ptrdiff_t difference_type;
  typedef value_type* pointer;
  typedef value_type& reference;

  typedef dualiterator<Value> self_type;
  typedef std::vector<value_type> singlevector_type;
  typedef std::vector<singlevector_type> dualvector_type;
 protected:
  dualvector_type& vec_;
  typename dualvector_type::iterator outer_;
  typename singlevector_type::iterator inner_;
 public:
  dualiterator(dualvector_type& vec, bool at_beginning)
      : vec_(vec)
  {
    if (at_beginning) {
      to_begin();
    } else {
      to_end();
    }
  }
  self_type& to_begin() {
    outer_ = vec_.begin();
    if (outer_ != vec_.end()) {
      inner_ = (*outer_).begin();
      while ((outer_ != vec_.end()) &&
             (inner_ == (*outer_).end())) {
        ++outer_;
        inner_ = (*outer_).begin();
      }
    }
    return *this;
  }
  self_type& to_end() {
    outer_ = vec_.end();
    return *this;
  }
  bool at_end() const {
    return outer_ == vec_.end();
  }
  reference operator*() {
    assert(!at_end());
    return *inner_;
  }
  pointer operator->() {
    assert(!at_end());
    return &(*inner_);
  }
  bool operator==(const self_type& rhs) const {
    return !operator!=(rhs);
  }
  bool operator!=(const self_type& rhs) const {
    return ((outer_ != rhs.outer_) ||
            (!at_end() && (inner_ != inner_)));
  }
  self_type& operator++() {
    assert(!at_end());
    ++inner_;
    while (!at_end() &&
           (inner_ == (*outer_).end())) {
      ++outer_;
      inner_ = (*outer_).begin();
    }
    return *this;
  }
  self_type operator++(int) {
    assert(!at_end());
    self_type tmp(*this);
    operator++();
    return tmp;
  }
};


/** Iterator template for the individual elements of a vector of vectors
 *
 * Use this to iterate over all the entries.
 */
template <class Value>
class const_dualiterator {
 public:
  typedef std::forward_iterator_tag iterator_category;
  typedef Value value_type;
  typedef std::ptrdiff_t difference_type;
  typedef const value_type* pointer;
  typedef const value_type& reference;

  typedef const_dualiterator<Value> self_type;
  typedef std::vector<value_type> singlevector_type;
  typedef std::vector<singlevector_type> dualvector_type;
 protected:
  const dualvector_type& vec_;
  typename dualvector_type::const_iterator outer_;
  typename singlevector_type::const_iterator inner_;
 public:
  const_dualiterator(const dualvector_type& vec, bool at_beginning)
      : vec_(vec)
  {
    if (at_beginning) {
      to_begin();
    } else {
      to_end();
    }
  }
  self_type& to_begin() {
    outer_ = vec_.cbegin();
    if (outer_ != vec_.cend()) {
      inner_ = (*outer_).cbegin();
      while ((outer_ != vec_.cend()) &&
             (inner_ == (*outer_).cend())) {
        ++outer_;
        inner_ = (*outer_).cbegin();
      }
    }
    return *this;
  }
  self_type& to_end() {
    outer_ = vec_.cend();
    return *this;
  }
  bool at_end() const {
    return outer_ == vec_.cend();
  }
  reference operator*() const {
    assert(!at_end());
    return *inner_;
  }
  pointer operator->() const {
    assert(!at_end());
    return &(*inner_);
  }
  bool operator==(const self_type& rhs) const {
    return !operator!=(rhs);
  }
  bool operator!=(const self_type& rhs) const {
    return ((outer_ != rhs.outer_) ||
            (!at_end() && (inner_ != inner_)));
  }
  self_type& operator++() {
    assert(!at_end());
    ++inner_;
    while (!at_end() &&
           (inner_ == (*outer_).cend())) {
      ++outer_;
      inner_ = (*outer_).cbegin();
    }
    return *this;
  }
  self_type operator++(int) {
    assert(!at_end());
    self_type tmp(*this);
    operator++();
    return tmp;
  }
};








  /** Class for handling proto-information for submesh/supermesh indices */
  class MPIIOBimeshProtoInfo {
   public:
    typedef Eigen::Matrix<int64_t,Eigen::Dynamic,1> triangle_storage_type;
    typedef Eigen::Matrix<int64_t,Eigen::Dynamic,3> edge_storage_type;
    typedef triangle_storage_type vertex_storage_type;

   protected:
    int my_rank_;
    /// The number of priorities used for triangles (normally 3)
    int64_t num_triangle_;
    /// The number of priorities used for edges (normally 4)
    int64_t num_edge_;
    /// The number of priorities used for vertices
    int64_t num_vertex_;
    /// The number of priorities (max of the other three)
    int64_t num_priorities_;
    triangle_storage_type triangle_;
    edge_storage_type edge_;
    vertex_storage_type vertex_;
   protected:
    /** Construct MPI I/O priorities for the micro level of a bilevel mesh
     *
     * @param comm The MPI communicator that will be used for communication
     * @param bimesh The bilevel mesh
     *
     */
    void construct_priorities(MPI_Comm comm, MeshBilevel& bimesh);
   public:
    /** Construct MPI I/O priorities for the micro level of a bilevel mesh
     *
     * @param comm The MPI communicator that will be used for communication
     * @param bimesh The bilevel mesh
     *
     */ 
    MPIIOBimeshProtoInfo(MPI_Comm comm, MeshBilevel& bimesh);

    int64_t num_priorities() const {
      return num_priorities_;
    }

    size_t triangle_size() const {
      return triangle_.rows();
    }
    
    size_t edge_size() const {
      return edge_.rows();
    }
    
    size_t vertex_size() const {
      return vertex_.rows();
    }

    triangle_storage_type& triangle() {
      return triangle_;
    }
    
    edge_storage_type& edge() {
      return edge_;
    }
    
    vertex_storage_type& vertex() {
      return vertex_;
    }
    
    int64_t triangle(int64_t tri) const {
      return triangle_(tri);
    }
    
    int64_t edge(int64_t tri, int64_t vtx) const {
      return edge_(tri, vtx);
    }
    
    int64_t vertex(int64_t vtx) const {
      return vertex_(vtx);
    }
    
  };





class MPIIOPriorities {
 public:

  class index_type {
   public:
    int64_t submesh_start;
    int64_t supermesh_start;
    int64_t count;
    index_type(int64_t thesubmesh_start, int64_t thesupermesh_start, int64_t thecount)
        : submesh_start(thesubmesh_start),
          supermesh_start(thesupermesh_start),
          count(thecount)
    {
    }
  };
  typedef std::vector<index_type> index_collection_type;
  typedef std::vector<index_collection_type> priority_index_type;

  typedef dualiterator<index_type> iterator;
  typedef const_dualiterator<index_type> const_iterator;

 protected:
  int my_rank_;
  /// Vector of length num_priorities_, each element a collection of
  /// index specifications
  priority_index_type indices_;

  /** Construct MPI I/O priorities for the micro level of a bilevel mesh
   *
   * @param comm The MPI communicator that will be used for communication
   * @param bimesh The bilevel mesh
   *
   */
  void construct_indices(MeshBilevel& bimesh, MPIIOBimeshProtoInfo& info);
 public:
  /** Construct MPI I/O priorities for the micro level of a bilevel mesh
   *
   * @param comm The MPI communicator that will be used for communication
   * @param bimesh The bilevel mesh
   *
   */ 
  MPIIOPriorities(MPI_Comm comm, MeshBilevel& bimesh);

  uint64_t size() const {
    return indices_.size();
  }

  uint64_t size(int64_t k) const {
    return indices_[k].size();
  }

  const priority_index_type& indices() const {
    return indices_;
  }

  priority_index_type& indices() {
    return indices_;
  }

  const index_collection_type& indices(int64_t k) const {
    return indices_[k];
  }

  index_collection_type& indices(int64_t k) {
    return indices_[k];
  }

  const index_type& indices(int64_t i, int64_t j) const {
    return indices_[i][j];
  }

  index_type& indices(int64_t i, int64_t j) {
    return indices_[i][j];
  }

  iterator begin() {
    return iterator(indices_, true);
  }
  iterator end() {
    return iterator(indices_, false);
  }

  const_iterator cbegin() const {
    return const_iterator(indices_, true);
  }
  const_iterator cend() const {
    return const_iterator(indices_, false);
  }


  void log_dump(std::ostream& os);
};




/** Indexing information for get/put_varm pnetcdf calls
 *
 * ncmpi_put_varm(ncid, varid,
 *                start.data(),
 *                count.data(),
 *                stride.data(),
 *                imap.data(),
 *                &buffer.data()[buffer_offset],
 *                bufcount,
 *                buftype)
 */
class MPIIOPrioritiesVarm {
 public:
  
  class index_type {
   public:
    typedef std::vector<MPI_Offset> vector_int;
    vector_int start;
    vector_int count;
    vector_int stride;
    vector_int imap;
    MPI_Offset buffer_offset; /// Index to the first element in the local buffer
    MPI_Offset bufcount;
    index_type(int n_dim)
        : start(n_dim),
          count(n_dim),
          stride(n_dim),
          imap(n_dim),
          buffer_offset(0),
          bufcount(0)
    {
    }
    int ndim() const {
      return start.size();
    }

    void construct(const MPIIOPriorities::index_type& blob,
                   SpaceTimeMapper& mapper,
                   const std::vector<bool>& active_dim,
                   const std::vector<int>& dim_idx,
                   int64_t num_fields);

    /** Iterate over the local buffer elements */
    template <class Value>
    class iterator {
     public:
      typedef std::forward_iterator_tag iterator_category;
      typedef Value value_type;
      typedef std::ptrdiff_t difference_type;
      typedef value_type* pointer;
      typedef value_type& reference;

      typedef iterator self_type;
     protected:
      index_type* info_;
      pointer ptr_;
      std::vector<int64_t> idx_;
     protected:
      void check_end() {
        if ((ptr_ != NULL) &&
            (idx_[0] >= info_->count[0])) {
          ptr_ = NULL;
        }
      }
      bool at_end() const {
        return (ptr_ == NULL);
      }
     public:
      iterator()
          : info_(),
            ptr_(NULL),
            idx_()
      {
      }
      iterator(index_type& info, pointer buf)
          : info_(&info),
            ptr_(NULL),
            idx_()
      {
        if (buf) {
          idx_.resize(info_->ndim(), 0);
          ptr_ = buf + info_->buffer_offset;
          check_end();
        }
      }
      reference operator*() {
        assert(!at_end());
        return *ptr_;
      }
      pointer operator->() {
        assert(!at_end());
        return &(*ptr_);
      }
      bool operator==(const self_type& rhs) const {
        if (at_end()) {
          return rhs.at_end();
        } else if (rhs.at_end()) {
          return false;
        }
        bool is_equal = true;
        for (int k = 0; is_equal && (k < info_->ndim()); k++) {
          is_equal = idx_[k] == rhs.idx_[k];
        }
        return is_equal;
      }
      bool operator!=(const self_type& rhs) const {
        if (at_end()) {
          return !rhs.at_end();
        } else if (rhs.at_end()) {
          return true;
        }
        bool is_equal = true;
        for (int k = 0; is_equal && (k < info_->ndim()); k++) {
          is_equal = idx_[k] == rhs.idx_[k];
        }
        return !is_equal;
      }
      self_type& operator++() {
        assert(!at_end());
        int k = info_->ndim() - 1;
        ++idx_[k];
        ptr_ += info_->imap[k];
        while ((k > 0) && (idx_[k] >= info_->count[k])) {
          ptr_ -= info_->imap[k] * idx_[k];
          idx_[k] = 0;
          --k;
          ++idx_[k];
          ptr_ += info_->imap[k];
        }
        check_end();
        return *this;
      }
      self_type operator++(int) {
        assert(!at_end());
        self_type tmp(*this);
        operator++();
        return tmp;
      }
    };

    /** Iterate over the local buffer elements */
    template <class Value>
    class const_iterator {
     public:
      typedef std::forward_iterator_tag iterator_category;
      typedef Value value_type;
      typedef std::ptrdiff_t difference_type;
      typedef const value_type* pointer;
      typedef const value_type& reference;

      typedef const_iterator self_type;
     protected:
      const index_type* info_;
      pointer ptr_;
      std::vector<int64_t> idx_;
     protected:
      void check_end() {
        if ((ptr_ != NULL) &&
            (idx_[0] >= info_->count[0])) {
          ptr_ = NULL;
        }
      }
      bool at_end() const {
        return (ptr_ == NULL);
      }
     public:
      const_iterator()
          : info_(),
            ptr_(NULL),
            idx_()
      {
      }
      const_iterator(const index_type& info, pointer buf)
          : info_(&info),
            ptr_(NULL),
            idx_()
      {
        if (buf) {
          idx_.resize(info_->ndim(), 0);
          ptr_ = buf + info_->buffer_offset;
          check_end();
        }
      }
      reference operator*() const {
        assert(!at_end());
        return *ptr_;
      }
      pointer operator->() const {
        assert(!at_end());
        return &(*ptr_);
      }
      bool operator==(const self_type& rhs) const {
        if (at_end()) {
          return rhs.at_end();
        } else if (rhs.at_end()) {
          return false;
        }
        bool is_equal = true;
        for (int k = 0; is_equal && (k < info_->ndim()); k++) {
          is_equal = idx_[k] == rhs.idx_[k];
        }
        return is_equal;
      }
      bool operator!=(const self_type& rhs) const {
        if (at_end()) {
          return !rhs.at_end();
        } else if (rhs.at_end()) {
          return true;
        }
        bool is_equal = true;
        for (int k = 0; is_equal && (k < info_->ndim()); k++) {
          is_equal = idx_[k] == rhs.idx_[k];
        }
        return !is_equal;
      }
      self_type& operator++() {
        assert(!at_end());
        int k = info_->ndim() - 1;
        ++idx_[k];
        ptr_ += info_->imap[k];
        while ((k > 0) && (idx_[k] >= info_->count[k])) {
          ptr_ -= info_->imap[k] * idx_[k];
          idx_[k] = 0;
          --k;
          ++idx_[k];
          ptr_ += info_->imap[k];
        }
        check_end();
        return *this;
      }
      self_type operator++(int) {
        assert(!at_end());
        self_type tmp(*this);
        operator++();
        return tmp;
      }
    };

    template <class Value>
    iterator<Value> begin(typename iterator<Value>::pointer buf) {
      return iterator<Value>(*this, buf);
    }
    template <class Value>
    iterator<Value> end(typename iterator<Value>::pointer buf) {
      return iterator<Value>(*this, NULL);
    }

    template <class Value>
    const_iterator<Value> cbegin(typename const_iterator<Value>::pointer buf) const {
      return const_iterator<Value>(*this, buf);
    }
    template <class Value>
    const_iterator<Value> cend(typename const_iterator<Value>::pointer buf) const {
      return const_iterator<Value>(*this, NULL);
    }

  };

  typedef std::vector<index_type> index_collection_type;
  typedef std::vector<index_collection_type> priority_index_type;


  /** Iterator template for the individual elements within one priority level
   *
   * Use this to iterate over all the buffer entries.
   */
  template <class Value>
  class valueiterator {
   public:
    typedef std::forward_iterator_tag iterator_category;
    typedef Value value_type;
    typedef std::ptrdiff_t difference_type;
    typedef value_type* pointer;
    typedef value_type& reference;

    typedef valueiterator<Value> self_type;
    typedef typename index_collection_type::iterator blob_iterator;
    typedef typename index_type::iterator<Value> value_iterator;
   protected:
    index_collection_type* blobs_;
    pointer buf_;
    blob_iterator blob_;
    value_iterator value_;
    value_iterator value_end_;
   public:
    valueiterator(index_collection_type& blobs, pointer buf)
        : blobs_(&blobs),
          buf_(buf)
    {
      if (buf) {
        to_begin();
      } else {
        to_end();
      }
    }
    self_type& next_blob() {
      ++blob_;
      if (!at_end()) {
        value_ = (*blob_).begin<Value>(buf_);
        value_end_ = (*blob_).end<Value>(buf_);
      }
      return *this;
    }
    self_type& to_begin() {
      blob_ = blobs_->begin();
      if (!at_end()) {
        value_ = (*blob_).begin<Value>(buf_);
        value_end_ = (*blob_).end<Value>(buf_);
        while (!at_end() &&
               (value_ == value_end_)) {
          next_blob();
        }
      }
      return *this;
    }
    self_type& to_end() {
      blob_ = blobs_->end();
      return *this;
    }
    bool at_end() const {
      return blob_ == blobs_->end();
    }
    reference operator*() {
      assert(!at_end());
      return *value_;
    }
    pointer operator->() {
      assert(!at_end());
      return &(*value_);
    }
    bool operator==(const self_type& rhs) const {
      return !operator!=(rhs);
    }
    bool operator!=(const self_type& rhs) const {
      return ((blob_ != rhs.blob_) ||
              (!at_end() && (value_ != value_)));
    }
    self_type& operator++() {
      assert(!at_end());
      ++value_;
      while (!at_end() &&
             (value_ == value_end_)) {
        next_blob();
      }
      return *this;
    }
    self_type operator++(int) {
      assert(!at_end());
      self_type tmp(*this);
      operator++();
      return tmp;
    }
  };


  /** Iterator template for the individual elements within one priority level
   *
   * Use this to iterate over all the buffer entries.
   */
  template <class Value>
  class const_valueiterator {
   public:
    typedef std::forward_iterator_tag iterator_category;
    typedef Value value_type;
    typedef std::ptrdiff_t difference_type;
    typedef const value_type* pointer;
    typedef const value_type& reference;

    typedef const_valueiterator<Value> self_type;
    typedef typename index_collection_type::const_iterator blob_iterator;
    typedef typename index_type::const_iterator<Value> value_iterator;
   protected:
    const index_collection_type* blobs_;
    pointer buf_;
    blob_iterator blob_;
    value_iterator value_;
    value_iterator value_end_;
   public:
    const_valueiterator(const index_collection_type& blobs, pointer buf)
        : blobs_(&blobs),
          buf_(buf)
    {
      if (buf) {
        to_begin();
      } else {
        to_end();
      }
    }
    self_type& next_blob() {
      ++blob_;
      if (!at_end()) {
        value_ = (*blob_).cbegin<Value>(buf_);
        value_end_ = (*blob_).cend<Value>(buf_);
      }
      return *this;
    }
    self_type& to_begin() {
      blob_ = blobs_->cbegin();
      if (!at_end()) {
        value_ = (*blob_).cbegin<Value>(buf_);
        value_end_ = (*blob_).cend<Value>(buf_);
        while (!at_end() &&
               (value_ == value_end_)) {
          next_blob();
        }
      }
      return *this;
    }
    self_type& to_end() {
      blob_ = blobs_->cend();
      return *this;
    }
    bool at_end() const {
      return blob_ == blobs_->cend();
    }
    reference operator*() const {
      assert(!at_end());
      return *value_;
    }
    pointer operator->() const {
      assert(!at_end());
      return &(*value_);
    }
    bool operator==(const self_type& rhs) const {
      return !operator!=(rhs);
    }
    bool operator!=(const self_type& rhs) const {
      return ((blob_ != rhs.blob_) ||
              (!at_end() && (value_ != value_)));
    }
    self_type& operator++() {
      assert(!at_end());
      ++value_;
      while (!at_end() &&
             (value_ == value_end_)) {
        next_blob();
      }
      return *this;
    }
    self_type operator++(int) {
      assert(!at_end());
      self_type tmp(*this);
      operator++();
      return tmp;
    }
  };


 protected:
  int my_rank_;
  priority_index_type indices_;
 public:
  MPIIOPrioritiesVarm() {};

  uint64_t size() const {
    return indices_.size();
  }

  void initialise(MPIIOPriorities& prio) {
    indices_.clear();
    indices_.resize(prio.size());
  }

  void construct(int64_t num_fields,
                 const std::vector<bool>& active_dim,
                 MPIIOPriorities& priorities,
                 SpaceTimeMapper& mapper);

  void log_dump(std::ostream& os);


  class iterator : public priority_index_type::iterator {
   public:
    iterator(priority_index_type::iterator iter) : priority_index_type::iterator(iter) {
    }
    template <class Value>
    valueiterator<Value> begin(typename valueiterator<Value>::pointer buf) {
      return valueiterator<Value>(operator*(), buf);
    }
    template <class Value>
    valueiterator<Value> end(typename valueiterator<Value>::pointer buf = NULL) {
      return valueiterator<Value>(operator*(), NULL);
    }
    template <class Value>
    const_valueiterator<Value> cbegin(typename const_valueiterator<Value>::pointer buf) const {
      return const_valueiterator<Value>(operator*(), buf);
    }
    template <class Value>
    const_valueiterator<Value> cend(typename const_valueiterator<Value>::pointer buf = NULL) const {
      return const_valueiterator<Value>(operator*(), NULL);
    }
  };

  class const_iterator : public priority_index_type::const_iterator {
   public:
    const_iterator(priority_index_type::iterator iter) : priority_index_type::const_iterator(iter) {
    }
    const_iterator(priority_index_type::const_iterator iter) : priority_index_type::const_iterator(iter) {
    }
    template <class Value>
    const_valueiterator<Value> cbegin(typename const_valueiterator<Value>::pointer buf) const {
      return const_valueiterator<Value>(operator*(), buf);
    }
    template <class Value>
    const_valueiterator<Value> cend(typename const_valueiterator<Value>::pointer buf = NULL) const {
      return const_valueiterator<Value>(operator*(), NULL);
    }
  };

  iterator begin() {
    return indices_.begin();
  }
  iterator end() {
    return indices_.end();
  }
  template <class Value>
  valueiterator<Value> begin(index_collection_type& blobs,
                             typename valueiterator<Value>::pointer buf) {
    return valueiterator<Value>(blobs, buf);
  }
  template <class Value>
  valueiterator<Value> end(index_collection_type& blobs,
                           typename valueiterator<Value>::pointer buf = NULL) {
    return valueiterator<Value>(blobs, NULL);
  }

  const_iterator cbegin() const {
    return indices_.begin();
  }
  const_iterator cend() const {
    return indices_.end();
  }
  template <class Value>
  const_valueiterator<Value> cbegin(const index_collection_type& blobs,
                                    typename const_valueiterator<Value>::pointer buf) const {
    return const_valueiterator<Value>(blobs, buf);
  }
  template <class Value>
  const_valueiterator<Value> cend(const index_collection_type& blobs,
                                  typename const_valueiterator<Value>::pointer buf = NULL) const {
    return const_valueiterator<Value>(blobs, NULL);
  }


};
































class MeshBilevel
{
  std::vector<int64_t> levels_;
  std::vector<fmesh::Mesh > meshes_;
  std::vector<fmesh::TriangleLocator* > locators_;
  TimerHierarchy* timer_;
  bool use_locators_;
  bool initialised_meshes_;
  bool initialised_locators_;

public:
  typedef typename MeshBilevelIndexMap::IndexMap IndexMap;
  MeshBilevelIndexMap* python_index_;
  MeshBilevelIndexMap* supermesh_index_;
 public:
  /** Number of interior points along a macro edge */
  int64_t macro_edge_size_;
  /** Number of interior points inside a macro triangle */
  int64_t macro_interior_size_;

 public:
  /**  Is the edge on the macro mesh boundary?
   *
   * TRUE if the macro edge is on the boundary.
   * The column index refers to opposing edge (edge 0 is the 1-2 edge, 1 is the 2-0 edge,
   * 2 is the 0-1 edge. */
  Eigen::Matrix<bool, Eigen::Dynamic, 3> macro_boundary_;
  /** Is the edge on a supermesh boundary?
   *
   * If mesh is associated with a supermesh, then TRUE if the macro
   * edge is the edge of the macro mesh of a super-mesh. Otherwise
   * equal to macro_boundary_.
   *
   * The column index refers to opposing
   * edge (edge 0 is the 1-2 edge, 1 is the 2-0 edge, 2 is the 0-1
   * edge. */
  Eigen::Matrix<bool, Eigen::Dynamic, 3> supermesh_macro_boundary_;

  /** Index of the micro vertex of each macro vertex */
  Eigen::Matrix<int64_t, Eigen::Dynamic, 1> macro_vertex_begin_;
  /** Index of the first micro vertex for each edge of a macro triangle */
  Eigen::Matrix<int64_t, Eigen::Dynamic, 3> macro_edge_begin_;
  /** Index of the first micro vertex for the interior of each macro triangle */
  Eigen::Matrix<int64_t, Eigen::Dynamic, 1> macro_interior_begin_;

  /** macro mesh vertex index collection of domain decomposition root nodes of some supermesh */
  std::set<int64_t> macro_domain_root_set_;
  /** If non-negative, identifies a macro domain decomposition root node */
  int64_t macro_domain_root_vertex_;

  Eigen::Matrix<double,Eigen::Dynamic,1> macro_weights_;
  Eigen::Matrix<double,Eigen::Dynamic,1> macro_halo_hat_;
  Eigen::Matrix<double,Eigen::Dynamic,1> macro_hat_;
  bool initialised_weights_;

  MPIIOPriorities* mpiio_priorities_;

public:
  MeshBilevel(TimerHierarchy* timer = NULL,
	      bool use_locators = false) :
    levels_(2), meshes_(2), locators_(2),
    timer_(timer),
    use_locators_(use_locators),
    initialised_meshes_(false), initialised_locators_(false),
    python_index_(NULL), supermesh_index_(NULL),
    macro_edge_size_(0), macro_interior_size_(0),
    macro_boundary_(0,3), supermesh_macro_boundary_(0,3),
    macro_edge_begin_(0,3), macro_interior_begin_(0),
    macro_domain_root_set_(),
    macro_domain_root_vertex_(-1),
    macro_weights_(0,1), macro_halo_hat_(0,1), macro_hat_(0,1),
    initialised_weights_(false),
    mpiio_priorities_(NULL)
  {
    meshes_[0].useVT(true);
    meshes_[1].useVT(true);
    meshes_[0].useTTi(true);
    meshes_[1].useTTi(true);
    locators_[0] = NULL;
    locators_[1] = NULL;
    set_levels(0, 1);
  }
  ~MeshBilevel() {
    clear_mpiio_priorities();
    clear_meshes();
  }

  MeshBilevel& set_levels(int64_t level0, int64_t level1)
  {
    assert(level0 <= level1);
    clear_meshes();
    levels_[0] = level0;
    levels_[1] = level1;
    macro_edge_size_ = (1 << (levels_[1]-levels_[0])) - 1;
    macro_interior_size_ = (macro_edge_size_ * (macro_edge_size_ - 1)) / 2;
    return *this;
  }

  void log_python_index() {
    if (python_index_) {
      PLOG_("PYTHON INDEX" << std::endl);
      python_index_->log_dump();
    }
  }
  void log_supermesh_index() {
    if (supermesh_index_) {
      PLOG_("SUPERMESH INDEX" << std::endl);
      supermesh_index_->log_dump();
    }
  }

  void log_dump() {
    log_python_index();
    log_supermesh_index();

    std::cout << std::boolalpha;
    PLOG_("LEVELS\t" << levels_[0] << ", " << levels_[1] << std::endl);
    //   std::vector<fmesh::Mesh > meshes_;
    //  std::vector<fmesh::TriangleLocator* > locators_;
    //  TimerHierarchy* timer_;
    PLOG_("USE LOCATORS\n" << use_locators_ << std::endl);
    PLOG_("INITIALISED MESHES\t" << initialised_meshes_ << std::endl);
    PLOG_("INITIALISED LOCATORS\t" << initialised_locators_ << std::endl);

    PLOG_("#POINTS OF MACRO EDGE\t" << macro_edge_size_ << std::endl);
    PLOG_("#POINTS IN TRIANGLE INTERIOR\t" << macro_interior_size_ << std::endl);

    PLOG_("IS OPPOSING MACRO EDGE A BOUNDARY?\n" << macro_boundary_ <<std::endl);
    PLOG_("IS OPPOSING SUPERMESH MACRO EDGE A BOUNDARY?\n" << supermesh_macro_boundary_ <<std::endl);
    PLOG_("MACRO VERTEX IS HERE\n" << macro_vertex_begin_ << std::endl);
    PLOG_("MACRO EDGE STARTS HERE\n" << macro_edge_begin_ << std::endl);
    PLOG_("MACRO TRIANGLE INTERIOR STARTS HERE\n" << macro_interior_begin_ << std::endl);

    PLOG_("MACRO VERTEX SET FOR DOMAIN ROOTS OF A SUPERMESH?\n" <<
          macro_domain_root_set_ << std::endl);

    /** If non-negative, identifies a macro domain decomposition root node */
    PLOG_("THE MACRO DOMAIN ROOT VERTEX\t" << macro_domain_root_vertex_ << std::endl);

    PLOG_("INITIALISED WEIGHTS\t" << initialised_weights_ << std::endl);
    //  Eigen::Matrix<double,Eigen::Dynamic,1> macro_weights_;
    //  Eigen::Matrix<double,Eigen::Dynamic,1> macro_halo_hat_;
    //  Eigen::Matrix<double,Eigen::Dynamic,1> macro_hat_;
    PLOG_("END OF LOG DUMP\n");
}
  

  const fmesh::Mesh& macro_mesh() {
    assert(initialised_meshes_);
    return meshes_[0];
  }
  const fmesh::Mesh& micro_mesh() {
    assert(initialised_meshes_);
    return meshes_[1];
  }
  int64_t macro_num_basis() const {
    assert(initialised_meshes_);
    return meshes_[0].nV();
  }
  int64_t micro_num_basis() const {
    assert(initialised_meshes_);
    return meshes_[1].nV();
  }
  int64_t super_macro_num_basis() const {
    assert(supermesh_index_);
    return supermesh_index_->macro_num_basis();
  }
  int64_t super_micro_num_basis() const {
    assert(supermesh_index_);
    return supermesh_index_->micro_num_basis();
  }
  Eigen::Matrix<double,Eigen::Dynamic,1>& macro_hat() {
    return macro_hat_;
  }
  Eigen::Matrix<double,Eigen::Dynamic,1>& macro_weights() {
    return macro_weights_;
  }
  Eigen::Matrix<double,Eigen::Dynamic,1>& macro_halo_hat() {
    return macro_halo_hat_;
  }

  /** Number of interior points along a macro triangle edge */
  int64_t macro_edge_size() const {
    return macro_edge_size_;
  }
  /** Number of interior points inside a macro triangle */
  int64_t macro_interior_size() const {
    return macro_interior_size_;
  }

  /** fmesh::Dart for edge opposing subvertex vtx of macro triangle tri. */
  fmesh::Dart macro_opposing_edge(int64_t tri, int64_t vtx) const {
    return fmesh::Dart(meshes_[0], tri, 1, (vtx + 1) % 3); // Opposing edge
  }
  /** fmesh::Dart for edge opposing subvertex vtx of micro triangle tri. */
  fmesh::Dart micro_opposing_edge(int64_t tri, int64_t vtx) const {
    return fmesh::Dart(meshes_[1], tri, 1, (vtx + 1) % 3); // Opposing edge
  }
  /** Check if an edge is "owned" by a triangle in a mesh,
   * disregarding any possible supermesh status. */
  bool triangle_is_edge_owner(const fmesh::Dart& dh) {
    return dh.onBoundary() || (dh.v() < dh.vo());
  }
  /** Check if the vtx-opposing edge is "owned" by the triangle tri in
   * the macro mesh, disregarding any possible supermesh status. */
  bool macro_triangle_is_edge_owner(int64_t tri, int64_t vtx) {
    return triangle_is_edge_owner(macro_opposing_edge(tri, vtx));
  }

  template <typename DerivedA, typename DerivedB, typename DerivedC>
  void macro_bary(const Eigen::DenseBase<DerivedA>& points,
		  Eigen::DenseBase<DerivedB> const & T_,
		  Eigen::DenseBase<DerivedC> const & bary_) {
    assert(initialised_meshes_);
    timer_tic("Map points directly to macro mesh");
    map_points_to_mesh(meshes_[0], locators_[0], points, T_, bary_);
    timer_toc();
  }
  template <typename DerivedA, typename DerivedB, typename DerivedC>
  void micro_bary(const Eigen::DenseBase<DerivedA>& points,
		  Eigen::DenseBase<DerivedB> const & T_,
		  Eigen::DenseBase<DerivedC> const & bary_) {
    assert(initialised_meshes_);
    timer_tic("Map points to micro mesh");
    map_points_to_mesh(meshes_[1], locators_[1], points, T_, bary_);
    timer_toc();
  }
  template <typename DerivedA, typename DerivedB, typename DerivedC,
	    typename DerivedD, typename DerivedE>
  void micro_bary(const Eigen::DenseBase<DerivedA>& points,
		  Eigen::DenseBase<DerivedB> const & T_,
		  Eigen::DenseBase<DerivedC> const & bary_,
		  Eigen::DenseBase<DerivedD> const & T0_,
		  Eigen::DenseBase<DerivedE> const & bary0_) {
    assert(initialised_meshes_);
    timer_tic("Map points to micro mesh");
    timer_tic("Find subtriangle and barycentric coordinates");
    map_points_to_mesh(meshes_[1], locators_[1], points, T_, bary_);
    timer_toc_tic("Calculate macro quasi-barycentric coordinates");
    map_tbary(levels_[1], levels_[0], T_, bary_, T0_, bary0_);
    timer_toc();
    timer_toc();
  }

  void import_python(const std::string& fileroot);
  void import_direct(const std::string& fileroot);
  void export_direct(const std::string& fileroot);

  void copy_python_index_to_submesh(MeshBilevel& super_mesh);

  std::set<int64_t> find_mesh_subset(std::set<int64_t>& macro_vertices);
  void define_subset_from_vertices(MeshBilevel& super_mesh,
                                   std::set<int64_t>& super_macro_vertices,
                                   bool claim_boundary);
  void define_subset_from_triangles(MeshBilevel& super_mesh,
                                    std::set<int64_t>& super_macro_triangles,
                                    bool claim_boundary);

  /** Internal method to be called when constructing meshes.
   *
   * Needs to be called after setting up macro_boundary etc, and before initialise_weights.
   *
   * Should also be called when constructing submeshes, but only when claiming boundaries.
   */
  void initialise_macro_domain_root_set();
  
  /** Internal method to be called when constructing submeshes without claiming boundaries.
   *
   * Copies the subdomain root vertex information from the super_mesh.
   *
   * Needs to be called after initialising supermesh_index_, and before initialise_weights.
   */
  void initialise_macro_domain_root_set(MeshBilevel& super_mesh);

  void clear_macro_domain_root_set() {
    macro_domain_root_set_.clear();
  }
  bool in_macro_domain_root_set(int64_t vertex) const {
    std::set<int64_t>::const_iterator it = macro_domain_root_set_.find(vertex);
    return it != macro_domain_root_set_.end();
  }
  std::set<int64_t>& macro_domain_root_set() {
    return macro_domain_root_set_;
  }
  
  
  template <typename DerivedA, typename DerivedB>
  void set_weights(int64_t vertex,
                   int64_t macro_tri,
                   int64_t corner,
                   double relative_halo_depth,
                   const Eigen::DenseBase< DerivedA > & triangle_macro_domain_root,
                   const Eigen::DenseBase< DerivedB > & qbary) {
    
    Eigen::Matrix<double,1,3> halo_hat;
    LOG("macro_boundary_.rows = " << macro_boundary_.rows() << std::endl);
    LOG("supermesh_macro_boundary_.rows = " << supermesh_macro_boundary_.rows() << std::endl);
    LOG("macro_tri = " << macro_tri << std::endl);
    Eigen::Matrix<bool,1,3> supermesh_macro_boundary_row = supermesh_macro_boundary_.row(macro_tri);
    LOG("calc_macro_halo_hat(relative_depth)" << std::endl);
    calc_macro_halo_hat(qbary, triangle_macro_domain_root,
                        supermesh_macro_boundary_row,
                        relative_halo_depth,
                        halo_hat);
    LOG("macro_halo_hat assignment" << std::endl);
    macro_halo_hat_(vertex) = halo_hat(corner);
    LOG("calc_macro_halo_hat(0.0)" << std::endl);
    calc_macro_halo_hat(qbary, triangle_macro_domain_root,
                        supermesh_macro_boundary_row,
                        0.0,
                        halo_hat);
    LOG("macro_weights calculation and assignment" << std::endl);
    macro_weights_(vertex) = calc_macro_weights(halo_hat, corner);
    LOG("macro_hat assignment" << std::endl);
    macro_hat_(vertex) = qbary(corner);
  }

  void initialise_weights(int64_t macro_domain_root_vertex = -1, double relative_halo_depth = 0.0) {
    assert(initialised_meshes_);
    LOG_("macro_edge_size_ = " << macro_edge_size_ << std::endl);
    LOG_("macro_interior_size_ = " << macro_interior_size_ << std::endl);
    LOG("macro_domain_root_vertex = " << macro_domain_root_vertex << std::endl);
    LOG("macro_domain_root_vertex_ = " << macro_domain_root_vertex_ << std::endl);
    if (macro_domain_root_vertex < 0) {
      macro_domain_root_vertex = macro_domain_root_vertex_;
    }
    LOG("macro_domain_root_vertex = " << macro_domain_root_vertex << std::endl);

    macro_weights_.resize(meshes_[1].nV(),1);
    macro_halo_hat_.resize(meshes_[1].nV(),1);
    macro_hat_.resize(meshes_[1].nV(),1);
    LOG("macro_weights size = " << macro_weights_.size() << std::endl);
    LOG("macro_halo_hat size = " << macro_halo_hat_.size() << std::endl);
    LOG("macro_hat size = " << macro_hat_.size() << std::endl);

    for (int64_t tri=0; tri < (int64_t)meshes_[0].nT(); ++tri) {
      int64_t corner0 = 0;
      if (macro_domain_root_vertex < 0) {
        corner0 = 3;
      } else {
        while (meshes_[0].TV(tri)[corner0] != macro_domain_root_vertex) {
          ++corner0;
          if (corner0 > 2) {
            break;
          }
        }
      }
      LOG("corner0 = " << corner0 << std::endl);

      if (corner0 > 2) {
	// Triangle not connected to the root node; fill with zeros.
	for (int64_t v=0; v < 3; ++v) {
	  macro_weights_(meshes_[0].TV(tri)[v]) = 0.0;
	  macro_halo_hat_(meshes_[0].TV(tri)[v]) = 0.0;
	  macro_hat_(meshes_[0].TV(tri)[v]) = 0.0;
	}
	for (int64_t e=0; e < 3; ++e) {
	  fmesh::Dart dh(meshes_[0], tri, 1, (e+1)%3);
	  for (int64_t v=macro_edge_begin_(tri,e);
	       v < macro_edge_begin_(tri,e) + macro_edge_size_;
	       ++v) {
	    macro_weights_(v) = 0.0;
	    macro_halo_hat_(v) = 0.0;
	    macro_hat_(v) = 0.0;
	  }
	}
	for (int64_t v = macro_interior_begin_(tri);
	     v < macro_interior_begin_(tri) + macro_interior_size_;
	     ++v) {
	  macro_weights_(v) = 0.0;
	  macro_halo_hat_(v) = 0.0;
	  macro_hat_(v) = 0.0;
	}
	continue;
      }

      LOG("Extract macro domain root information" << std::endl);
      Eigen::Matrix<double,1,3> triangle_macro_domain_root;
      triangle_macro_domain_root(0) = in_macro_domain_root_set(meshes_[0].TV(tri)[0]);
      triangle_macro_domain_root(1) = in_macro_domain_root_set(meshes_[0].TV(tri)[1]);
      triangle_macro_domain_root(2) = in_macro_domain_root_set(meshes_[0].TV(tri)[2]);

      Eigen::Matrix<double,1,3> qbary;
      fmesh::Dart dh0(meshes_[1],
		      get_micro_triangle_index(tri, levels_[0], levels_[1], corner0),
		      1,
		      corner0);
      qbary << 0.0, 0.0, 0.0;
      qbary(corner0) = 1.0;
      LOG("Construct weight for root node" << std::endl);
      set_weights(dh0.v(), tri, corner0, relative_halo_depth, triangle_macro_domain_root, qbary);
        
      LOG("Construct edge and inner weights" << std::endl);
      for (int64_t idx0=macro_edge_size_; idx0 >= 0; --idx0) {
	dh0.orbit2();
	fmesh::Dart dh = dh0;
	for (int64_t idx1=macro_edge_size_ + 1 - idx0; idx1 >= 0; --idx1) {
	  int64_t idx2 = macro_edge_size_ + 1 - idx0 - idx1;
	  qbary(corner0) = double(idx0) / double(macro_edge_size_ + 1);
          qbary((corner0+1)%3) = double(idx1) / double(macro_edge_size_ + 1);
          qbary((corner0+2)%3) = double(idx2) / double(macro_edge_size_ + 1);
          set_weights(dh.v(), tri, corner0, relative_halo_depth, triangle_macro_domain_root, qbary);
	  if (idx1 > 0) {
	    dh.orbit2().orbit0rev().orbit0rev();
	  }
	}
	if (idx0 > 0) {
	  dh0.orbit0rev().orbit0rev();
	}
      }
    }
    initialised_weights_ = true;
  }

  void clear_mpiio_priorities() {
    if (mpiio_priorities_) {
      delete mpiio_priorities_;
      mpiio_priorities_ = NULL;
    }
  }
/** Construct MPI I/O priorities for the micro level of a bilevel mesh
 *
 * Contructs an MPIIOPriorities object and stores it internally.
 * Access through mpiio_priorities().  Check for existance with has_mpiio_priorities().
 *
 * @param comm The MPI communicator that will be used for communication
 *
 * This is a collective MPI operation. (I think; all nodes in the communicator must call the method)
 */
  void construct_mpiio_priorities(MPI_Comm comm);
  bool has_mpiio_priorities() const {
    return mpiio_priorities_ != NULL;
  }
  MPIIOPriorities& mpiio_priorities() {
    assert(mpiio_priorities_ != NULL);
    return *mpiio_priorities_;
  }
  
protected:
  TIMER_HELPERS

public:
  void clear_weights() {
    macro_weights_.resize(0,1);
    macro_halo_hat_.resize(0,1);
    macro_hat_.resize(0,1);
    initialised_weights_ = false;
  }

protected:
  void initialise_python_index() {
    if (python_index_ != NULL) {
      python_index_->clear();
    } else {
      python_index_ = new MeshBilevelIndexMap;
    }
  }

  void initialise_supermesh_index() {
    if (supermesh_index_ != NULL) {
      supermesh_index_->clear();
    } else {
      supermesh_index_ = new MeshBilevelIndexMap;
    }
  }

  void clear_python_index() {
    if (python_index_ != NULL) {
      python_index_->clear();
      delete python_index_;
      python_index_ = NULL;
    }
  }
 
  void clear_supermesh_index() {
    if (supermesh_index_ != NULL) {
      supermesh_index_->clear();
      delete supermesh_index_;
      supermesh_index_ = NULL;
    }
  }
 
  void clear_meshes() {
    clear_locators();
    clear_weights();
    if (initialised_meshes_) {
      clear_python_index();
      clear_supermesh_index();
      meshes_[0].clear();
      meshes_[1].clear();
      initialised_meshes_ = false;
    }
    macro_domain_root_vertex_ = -1;
  }
 
  void clear_locators() {
    if (initialised_locators_) {
      delete locators_[0];
      locators_[0] = NULL;
      delete locators_[1];
      locators_[1] = NULL;
      initialised_locators_ = false;
    }
  }
  void initialise_locators() {
    assert(initialised_meshes_);
    clear_locators();
    if (use_locators_) {
      if (timer_ != NULL) timer_->tic("Initialising locators");
      initialise_locators_helper(locators_, meshes_, timer_);
      if (timer_ != NULL) timer_->toc();
      initialised_locators_ = true;
    }
  }


  /** Reorganise the finer scale mesh into more convenient node order:
   * Macro nodes
   * Nodes along macro edge
   * Nodes strictly inside a macro triangle
   */
  void reorganise_from_python(const fmesh::Mesh& mesh_py);

};



class DateTimeHelper {
  EUSTACE::CalendarDay epoch_;
  EUSTACE::TimeBaseDays calendar_;
 public:
  DateTimeHelper(int64_t epoch_calyear) : epoch_(epoch_calyear, 1, 1), calendar_(epoch_) {
  }

  int64_t day_from_other(int64_t day_other, const DateTimeHelper& dt_other) {
    return day_other - dt_other.get_day(epoch_calyear(), 1, 1);
  }

  int64_t whole_of_frac(double frac) const {
    return int64_t(std::floor(frac));
  }
  double fraction_of_frac(double frac) const {
    return (frac - std::floor(frac));
  }

  int64_t epoch_calyear() const {
    return (epoch_.tm_year + 1900);
  }
  std::tm get_tm(int64_t day_number) const {
    return calendar_.Time(day_number);
  }
  int64_t get_calyear(int64_t day_number) const {
    std::tm now = calendar_.Time(day_number);
    return now.tm_year + 1900;
  }
  int64_t get_calyear(double frac_day) const {
    return get_calyear(whole_of_frac(frac_day));
  }
  int64_t get_year(double frac_day) const {
    return get_calyear(whole_of_frac(frac_day)) - epoch_calyear();
  }
  int64_t get_day(int64_t year, int64_t month, int64_t day) const {
    EUSTACE::CalendarDay the_day = EUSTACE::CalendarDay(year, month, day);
    return calendar_.Number(the_day);
  }
  int64_t days_in_calyear(int64_t calyear) const {
    return get_day(calyear + 1, 1, 1) - get_day(calyear, 1, 1);
  }
  int64_t get_yday(int64_t day_number) const {
    std::tm now = calendar_.Time(day_number);
    return now.tm_yday;
  }
  double get_yday(double frac_day) const {
    return get_yday(whole_of_frac(frac_day)) + fraction_of_frac(frac_day);
  }
  double fraction_of_year(double frac_day) const {
    return get_yday(frac_day) / double(days_in_calyear(get_calyear(frac_day)));
  }
  double fraction_of_year(int64_t day_number, double fraction_of_day) const {
    return (double(get_yday(day_number)) + fraction_of_day) /
        double(days_in_calyear(get_calyear(day_number)));
  }

  double fractional_day(double frac_year) const {
    int64_t calyear = whole_of_frac(frac_year) + epoch_calyear();
    return double(get_day(calyear, 1, 1)) +
        double(days_in_calyear(calyear)) * fraction_of_frac(frac_year);
  }
  int64_t whole_day(double frac_year) const {
    return whole_of_frac(fractional_day(frac_year));
  }

  double fractional_day(int64_t day_number, double fraction_of_day) const {
    return double(day_number) + fraction_of_day;
  }
  double fractional_year(int64_t day_number, double fraction_of_day) const {
    return double(get_year(day_number)) + fraction_of_year(day_number, fraction_of_day);
  }
  double fractional_year(double frac_day) const {
    return double(get_year(frac_day)) + fraction_of_year(frac_day);
  }
};



enum class TimeBasis { Bspline2, Harmonic };

class TimeMeshUnitless {
 public:
  class SuperMeshInfo {
    friend class TimeMeshUnitless;
   protected:
    int64_t subdomain_intervals_;
    int64_t shift_size_;
    int64_t shift_;
   public:
    SuperMeshInfo()
        : subdomain_intervals_(0),
          shift_size_(0),
          shift_(0)
    {
    }
    SuperMeshInfo(int64_t subdomain_intervals,
                  int64_t shift_size,
                  int64_t shift)
        : subdomain_intervals_(subdomain_intervals),
          shift_size_(shift_size),
          shift_(shift)
    {
    }
    int64_t subdomain_intervals() const {
      return subdomain_intervals_;
    }
    int64_t shift_size() const {
      return shift_size_;
    }
    int64_t shift() const {
      return shift_;
    }
    SuperMeshInfo& step() {
      shift_++;
      return *this;
    }
    int64_t offset() const {
      return shift_size_ * shift_;
    }
    int64_t next_offset() const {
      return offset() + shift_size_;
    }
    
  };
 protected:
  double domain_begin_;
  double domain_end_;
  double step_size_;
  int64_t num_intervals_;
  int64_t num_basis_;
  TimeBasis basis_;
  bool is_cyclic_;
  TimeMeshUnitless* super_mesh_;
  SuperMeshInfo super_mesh_info_;
 public:
  /** Create a temporal mesh
   * @param num_intervals Number of intervals for Bspline2, maximum number of orders for Harmonic (not counting the constant), see below.
   *
   * For Bspline2, the number of basis functions is num_intervals + 2*(!is_cyclic)
   *
   * For Harmonic, the number of basis functions is 1 + num_intervals*2, so that the cyclic
   * basis on a \f$[0,2\pi)\f$ interval is \f$\{1,\cos(t),\sin(t),\dots\}\f$, and the non-cyclic
   * basis is \f$\{1,\cos(t/2),\sin(t/2),\dots\}\f$
   */
  TimeMeshUnitless(double domain_begin, double domain_end, int64_t num_intervals,
                   TimeBasis basis, bool is_cyclic);
  TimeMeshUnitless(TimeMeshUnitless* the_super_mesh,
                   const SuperMeshInfo& info,
                   int64_t extra_shift);


  double domain_begin() const
  {
    return domain_begin_;
  }
  double domain_end() const
  {
    return domain_end_;
  }
  double step_size() const
  {
    return step_size_;
  }
  int64_t num_intervals() const
  {
    return num_intervals_;
  }
  int64_t num_basis() const 
  {
    return num_basis_;
  }
  TimeBasis get_basis() const 
  {
    return basis_;
  }
  bool is_cyclic() const 
  {
    return is_cyclic_;
  }

  TimeMeshUnitless* super_mesh_unitless() {
    return super_mesh_;
  }
  const SuperMeshInfo& super_info() const {
    assert(super_mesh_);
    return super_mesh_info_;
  }
  int64_t super_offset() const {
    if (super_mesh_) {
      return super_mesh_info_.offset();
    }
    return 0;
  }

  /** Compute the index of the interval in which the time falls */
  int64_t index(double the_time) {
    return int64_t(std::floor((the_time - domain_begin_) / step_size_));
  }

  typedef std::vector<std::pair<int64_t, double>> eval_type;
  void eval(double the_time, eval_type& output);

  void get_weights(int64_t halo_depth,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>& halo_hat,
                   Eigen::Matrix<double, Eigen::Dynamic, 1>& weights);
};



class TimeMesh : public TimeMeshUnitless {
  TimeUnits time_units_;
  DateTimeHelper dt_helper_;
 public:
  TimeMesh(double domain_begin, double domain_end, int64_t num_intervals,
           int64_t epoch_year, TimeUnits time_units, TimeBasis basis, bool is_cyclic) :
      TimeMeshUnitless(domain_begin, domain_end, num_intervals, basis, is_cyclic),
      time_units_(time_units),
      dt_helper_(epoch_year)
  {
  }
  TimeMesh(TimeMesh* the_super_mesh,
           const TimeMeshUnitless::SuperMeshInfo& info,
           int64_t extra_shift)
      : TimeMeshUnitless(the_super_mesh, info, extra_shift),
        time_units_(the_super_mesh->time_units_),
        dt_helper_(the_super_mesh->epoch_calyear())
  {
  }
        
  TimeMesh* super_mesh() {
    return static_cast<TimeMesh*>(super_mesh_);
  }

  TimeUnits time_units() {
    return time_units_;
  }
  int64_t epoch_calyear() {
    return dt_helper_.epoch_calyear();
  }

  int64_t index_day(double frac_day) {
    switch (time_units_) {
      case TimeUnits::Day : return index(frac_day); break;
      case TimeUnits::Year : return index(dt_helper_.fractional_year(frac_day)); break;
      case TimeUnits::Season : return index(dt_helper_.fraction_of_year(frac_day)); break;
    }
    return 0;
  }
  int64_t index_year(double frac_year) {
    switch (time_units_) {
      case TimeUnits::Day : return index(dt_helper_.fractional_day(frac_year)); break;
      case TimeUnits::Year : return index(frac_year); break;
      case TimeUnits::Season : return index(dt_helper_.fraction_of_frac(frac_year)); break;
    }
    return 0;
  }
    
  void eval_day(double frac_day, eval_type& output) {
    switch (time_units_) {
      case TimeUnits::Day : return eval(frac_day, output); break;
      case TimeUnits::Year : return eval(dt_helper_.fractional_year(frac_day), output); break;
      case TimeUnits::Season : return eval(dt_helper_.fraction_of_year(frac_day), output); break;
    }
  }
  void eval_year(double frac_year, eval_type& output) {
    switch (time_units_) {
      case TimeUnits::Day : return eval(dt_helper_.fractional_day(frac_year), output); break;
      case TimeUnits::Year : return eval(frac_year, output); break;
      case TimeUnits::Season : return eval(dt_helper_.fraction_of_frac(frac_year), output); break;
    }
  }

  /** Return the first and last day properly covered by the domain */
  DaySpan day_span() {
    switch (time_units_) {
      case TimeUnits::Day : return DaySpan(std::floor(domain_begin()),
                                           std::ceil(domain_end()) - 1); break;
      case TimeUnits::Year : return DaySpan(dt_helper_.whole_day(domain_begin()),
                                            dt_helper_.whole_day(domain_end())); break;
      case TimeUnits::Season : return DaySpan(true); break;
    }
    return DaySpan(false);
  }
  /** Return the first and last day covered by the complete basis function support */
  DaySpan day_span_full() {
    switch (time_units_) {
      case TimeUnits::Day : return DaySpan(std::floor(domain_begin() - step_size_),
                                           std::ceil(domain_end() + step_size_) - 1); break;
      case TimeUnits::Year : return DaySpan(dt_helper_.whole_day(domain_begin() - step_size_),
                                            dt_helper_.whole_day(domain_end() + step_size_)); break;
      case TimeUnits::Season : return DaySpan(true); break;
    }
    return DaySpan(false);
  }


};







/** Map points in space&time to nodes in a time-space-season block */
class TimeMapper {
 public:
  typedef double key_type;
  typedef TimeMesh::eval_type info_type;
  typedef std::unordered_map<key_type, info_type* > map_type;
 private:
  /** Temporal (daily or seasonal) mesh. May be NULL. */
  TimeMesh* mesh_;
  DateTimeHelper dt_helper_;
  bool is_seasonal_;

  map_type map_;
 public:
  TimeMapper(TimeMesh* mesh, int64_t epoch_calyear) :
      mesh_(mesh),
      dt_helper_(epoch_calyear),
      is_seasonal_(false),
      map_()
  {
    if (mesh_) {
      is_seasonal_ = (mesh_->time_units() == TimeUnits::Season);
    }
  }

  TimeMesh* mesh() {
    return mesh_;
  }
  bool is_seasonal() {
    return is_seasonal_;
  }
  int64_t size() {
    return map_.size();
  }

  void clear() {
    for (map_type::iterator iter = map_.begin();
         iter != map_.end();
         ++iter) {
      delete iter->second;
    }
    map_.clear();
  }
  
  ~TimeMapper() {
    clear();
  }

  /** Return the number of basis functions */
  int64_t num_basis() const {
    if (mesh_) {
      return mesh_->num_basis();
    }
    return 1;
  }

  int64_t super_offset() const {
    if (mesh_) {
      return mesh_->super_offset();
    }
    return 0;
  }


  key_type key_day(double frac_day);
  key_type key_year(double frac_year);

  /** Get info for point evaluation
   * @param stored true if the info object is stored in the map,
   *   false if the caller is responsible for deallocating it.
   */
  info_type* mapping(key_type key,
                     bool remember,
                     bool& stored);
  
};



    


/** Map spatial points and cells to mesh nodes */
class SpaceMapper {
 public:
  typedef std::pair<double, double> key_type;

  struct key_type_hasher {
    size_t operator()(const key_type& key) const {
      return (std::hash<double>()(key.first)
              ^ (std::hash<double>()(key.second) << 1));
    }
  };
  
  typedef std::map<std::int64_t, double> info_type;
  typedef std::unordered_map<key_type, info_type*, key_type_hasher > map_type;
 private:
  /** Spatial mesh. May be NULL. */
  MeshBilevel* mesh_;
  bool is_cellular_;
  double cell_width_;
  std::vector<key_type> bbox_latlong_; // [0] = min, [1] = max
  double bbox_long_centre_;
  
  map_type map_;

 public:
  SpaceMapper(MeshBilevel* mesh, double cell_width) :
      mesh_(mesh),
      is_cellular_(cell_width > 0.0),
      cell_width_(cell_width),
      bbox_latlong_(2)
  {
    make_bbox();
  }

  void make_bbox();

  bool in_bbox(const key_type& latlong) {
    double longitude = latlong.second - bbox_long_centre_;
    while (longitude < 180.0) {
      longitude += 360.0;
    }
    while (longitude > 180.0) {
      longitude -= 360.0;
    }
    return ((bbox_latlong_[0].first <= latlong.first) &&
            (bbox_latlong_[0].second <= longitude) &&
            (bbox_latlong_[1].first >= latlong.first) &&
            (bbox_latlong_[1].second >= longitude));
  }

  MeshBilevel* mesh() {
    return mesh_;
  }
  bool is_cellular() {
    return is_cellular_;
  }
  bool cell_width() {
    return cell_width_;
  }
  int64_t size() {
    return map_.size();
  }

  void clear() {
    for (map_type::iterator iter = map_.begin();
         iter != map_.end();
         ++iter) {
      delete iter->second;
    }
    map_.clear();
  }
  
  ~SpaceMapper() {
    clear();
  }

  /** Return the number of basis functions */
  int64_t num_basis() {
    if (mesh_) {
      return mesh_->micro_mesh().nV();
    }
    return 1;
  }
  /** Return the number of basis functions in the supermesh */
  int64_t super_num_basis() {
    if (mesh_) {
      return mesh_->micro_mesh().nV();
    }
    return 1;
  }


  key_type key_space(double lat, double lon) {
    if (!mesh_) {
      return key_type(0.0, 0.0);
    } else {
      return key_type(lat, lon);
    }
  }

  void eval(key_type key, info_type& info);

  /** Get info for point evaluation
   * @stored returns true if the info object is stored in the map,
   *   false if the caller is responsible for deallocating it.
   */
  info_type* mapping(key_type key,
                     bool remember,
                     bool& stored);
  
};


/** Map points in space&time to nodes in a time-space-season block */
class SpaceTimeMapper {
 public:
  class IsLocal {
   public:
    bool time_;
    bool space_point_;
    bool space_cell_;
    bool season_;
    IsLocal(bool time, bool space_point, bool space_cell, bool season) :
        time_(time), space_point_(space_point), space_cell_(space_cell), season_(season) {
    }
    IsLocal() :
        time_(false), space_point_(false), space_cell_(false), season_(false) {
    }
  };
  class Maps {
   public:
    /** Temporal basis mapping. Mesh may be NULL. */
    TimeMapper* time_;
    /** Spatial basis mapping for points. Mesh may be NULL. */
    SpaceMapper* space_point_;
    /** Spatial basis mapping for cells. Mesh may be NULL. */
    SpaceMapper* space_cell_;
    /** Seasonal basis mapping. Mesh may be NULL. */
    TimeMapper* season_;
    /** Keeps track of locally allocated mappers. */
    IsLocal is_local_;
    Maps(TimeMapper* time, SpaceMapper* space_point, SpaceMapper* space_cell, TimeMapper* season) :
        time_(time), space_point_(space_point), space_cell_(space_cell), season_(season),
        is_local_()
    {
    }
    Maps() :
        time_(NULL), space_point_(NULL), space_cell_(NULL), season_(NULL), is_local_()
    {
    }

    /** Construct locally handled mappers for unused dimensions */
    void make_local(int64_t epoch_calyear);

    void clear();
    
    ~Maps() {
      clear();
    }
  };

  DateTimeHelper dt_helper_;
  Maps maps_;

 public:
  /** Verifies the basic internal state consistency.
   *
   * - Presence of mappers
   * - time_units of the time mapper is TimeUnits::Day or TimeUnits::Year
   * - The number of basis functions of the point and cell spatial mappers match eachother.
   * - time_units of the season mapper is TimeUnits::Season
   */
  void verify_basics();
  const DateTimeHelper& dt_helper() const {
    return dt_helper_;
  }
  int64_t epoch_calyear() {
    return dt_helper_.epoch_calyear();
  }
  void clear() {
    maps_.clear();
  }
  void set(TimeMapper* time_map,
           SpaceMapper* space_point_map,
           SpaceMapper* space_cell_map,
           TimeMapper* season_map,
           int64_t epoch_calyear);
  SpaceTimeMapper(TimeMapper* time_map,
                  SpaceMapper* space_point_map,
                  SpaceMapper* space_cell_map,
                  TimeMapper* season_map,
                  int64_t epoch_calyear) :
      dt_helper_(epoch_calyear),
      maps_(time_map, space_point_map, space_cell_map, season_map)
  {
    maps_.make_local(dt_helper_.epoch_calyear());
    verify_basics();
  }
  ~SpaceTimeMapper() {
    clear();
  }

  Maps& maps() {
    return maps_;
  }

  int64_t num_basis_time() const {
    assert(maps_.time_ != NULL);
    return maps_.time_->num_basis();
  }
  int64_t num_basis_space() const {
    assert(maps_.space_point_ != NULL);
    return maps_.space_point_->num_basis();
  }
  int64_t num_basis_season() const {
    assert(maps_.season_ != NULL);
    return maps_.season_->num_basis();
  }

  int64_t super_offset_time() const {
    return maps_.time_->super_offset();
  }

  int64_t super_num_basis_time() const {
    assert(maps_.time_ != NULL);
    assert(maps_.time_->mesh() != NULL);
    assert(maps_.time_->mesh()->super_mesh() != NULL);
    return maps_.time_->mesh()->super_mesh()->num_basis();
  }
  int64_t super_num_basis_space() const {
    assert(maps_.space_point_ != NULL);
    return maps_.space_point_->num_basis();
  }
  /** Seasons can't have super meshes.
   * @return The regular num_basis_season output
   */
  int64_t super_num_basis_season() const {
    return num_basis_season();;
  }

  /** Return the aggregated number of basis functions
   * (i.e. the number of space-time-season basis nodes)
   */
  int64_t num_basis() const {
    return num_basis_time() * num_basis_space() * num_basis_season();
  }
  /** Return the node index of a field/time/space/season point */
  int64_t index(int64_t index_field, int64_t index_time, int64_t index_space, int64_t index_season) const;
  
  /** Return the aggregate supermesh node index of a field/time/space/season point */
  int64_t super_index(int64_t index_field, int64_t index_time, int64_t index_space, int64_t index_season) const;

  /** Get info for a time-space-season point evaluation */
  class Point {
   public:
    double lat_;
    double lon_;
    int64_t day_;
    bool is_cell_;
    Point() :
        lat_(0.0), lon_(0.0), day_(0), is_cell_(false) {
    }
    Point(double lat,
          double lon) :
        lat_(lat), lon_(lon), day_(0), is_cell_(false) {
    }
    Point(double lat,
          double lon,
          int64_t day) :
        lat_(lat), lon_(lon), day_(day), is_cell_(false) {
    }
    Point(double lat,
          double lon,
          int64_t day,
          bool is_cell) :
        lat_(lat), lon_(lon), day_(day), is_cell_(is_cell) {
    }
    void set(double lat,
             double lon,
             int64_t day,
             bool is_cell) {
      lat_ = lat;
      lon_ = lon;
      day_ = day;
      is_cell_ = is_cell;
    }
  };

  typedef std::pair<int64_t, int64_t> offset_type;
  /** Get info for point evaluation or cell average evaluation.
   * @param point The index to use as the triplet output row index
   * @param use_solar_time true for local solar time, false for UTC=0 time
   * @param row_index The index to use as the triplet output row index
   * @param remember true if the mappers should remember the point
   * @param field_weights Vector of weights for each field.
   *    If NULL, {1.0} is used instead.
   * @return The number of weights added to the output
   *
   * The dimension order is field, time, space, season, where the season index varies the fastest.
   */
  int64_t mapping(const Point& point, bool use_solar_time,
                  const offset_type& offset, bool remember,
                  const std::vector<double>* field_weights,
                  std::vector<Eigen::Triplet<double>>& output);
};












#include "ambitious/ncdf/ncdf.hpp"


/** Initiate iput_varm for all blobs in a collection
 */
template <class Value>
NCHelper& iput_varm(NCHelper& nch,
                    int varid,
                    MPIIOPrioritiesVarm& object,
                    const Value* buf) {
  for (auto iter = object.begin();
       iter != object.end();
       ++iter) {
    iput_varm(nch, varid, *iter, buf);
  }
  return nch;
}

/** Initiate iget_varm for all blobs in a collection
 */
template <class Value>
NCHelper& iget_varm(NCHelper& nch,
                    int varid,
                    MPIIOPrioritiesVarm& object,
                    Value* buf) {
  for (auto iter = object.begin();
       iter != object.end();
       ++iter) {
    iget_varm(nch, varid, *iter, buf);
  }
  return nch;
}

/** Initiate iput_varm for all blobs in a collection
 */
template <class Value>
NCHelper& iput_varm(NCHelper& nch,
                    int varid,
                    typename MPIIOPrioritiesVarm::index_collection_type& object,
                    const Value* buf) {
  for (auto iter = object.begin();
       iter != object.end();
       ++iter) {
    iput_varm(nch, varid, *iter, buf);
  }
  return nch;
}

/** Initiate iget_varm for all blobs in a collection
 */
template <class Value>
NCHelper& iget_varm(NCHelper& nch,
                    int varid,
                    typename MPIIOPrioritiesVarm::index_collection_type& object,
                    Value* buf) {
  for (auto iter = object.begin();
       iter != object.end();
       ++iter) {
    iget_varm(nch, varid, *iter, buf);
  }
  return nch;
}

/** Initiate iput_varm for a blob */
template <class Value>
NCHelper& iput_varm(NCHelper& nch,
                    int varid,
                    typename MPIIOPrioritiesVarm::index_type& blob,
                    const Value* buf) {
  nch.iput_varm(varid,
                blob.start,
                blob.count,
                blob.stride,
                blob.imap,
                buf + blob.buffer_offset);
  return nch;
}

/** Initiate iget_varm for a blob */
template <class Value>
NCHelper& iget_varm(NCHelper& nch,
                    int varid,
                    typename MPIIOPrioritiesVarm::index_type& blob,
                    Value* buf) {
  nch.iget_varm(varid,
                blob.start,
                blob.count,
                blob.stride,
                blob.imap,
                buf + blob.buffer_offset);
  return nch;
}


template <class Value, class Object>
NCHelper& iput_varm(NCHelper& nch,
               const std::string& varid,
               Object& object,
               const Value* buf) {
  return iput_varm(nch, nch.varid(varid), object, buf);
}
template <class Value, class Object>
NCHelper& iget_varm(NCHelper& nch,
                    const std::string& varid,
                    Object& object,
                    Value* buf) {
  return iget_varm(nch, nch.varid(varid), object, buf);
}













/**
Blocks: Daily, Slow, Season, Systematic
Operations: Apply Q, Apply M, ...

For each day:
  Create all new blocks that start on this day
  For each data Source:
    Load day of data
    Accumulate data info for each Block
  For each Block:     
    If all data info has been accumulated, perform Operation

*/








#endif
