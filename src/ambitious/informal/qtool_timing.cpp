#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <string>
#include <iostream>
#include <boost/random.hpp>
#include <boost/math/distributions/normal.hpp>

#include "ambitious/timer/timer.hpp"
#include "ambitious/poq/poq.hpp"
#include "ambitious/qtool/qtool.hpp"
#include <fmesher.hh>

#include "ambitious/informal/qtool_timing_Q_blocks.cpp"

typedef double MyScalar;
typedef Eigen::Matrix<MyScalar, Eigen::Dynamic, 1> Vector;
typedef Eigen::Matrix<MyScalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef fmesh::Matrix<MyScalar> FMatrix;

int pow2(unsigned int exponent)
{
  if (exponent == 0) {
    return 1;
  } else {
    return 2 << (exponent-1);
  }
}


fmesh::TriangleLocator* initialise_locator_local(const fmesh::Mesh& M)
{
  int the_dimensions[] = {0,1};
  std::vector<int> dimensions(the_dimensions,
			      the_dimensions +
			      sizeof(the_dimensions) / sizeof(int) );
  return new fmesh::TriangleLocator(&M, dimensions, true);
}


template <typename DerivedA, typename DerivedB, typename DerivedC>
void map_points_to_mesh(const fmesh::Mesh& M,
			const fmesh::TriangleLocator& locator,
			const Eigen::DenseBase<DerivedA>& points,
			Eigen::DenseBase<DerivedB> const & point2T_,
			Eigen::DenseBase<DerivedC> const & point2bary_)
{
  Eigen::DenseBase<DerivedB>& point2T = const_cast< Eigen::DenseBase<DerivedB>& >(point2T_);
  Eigen::DenseBase<DerivedC>& point2bary = const_cast< Eigen::DenseBase<DerivedC>& >(point2bary_);
  point2T.derived().resize(points.rows(), 1);
  point2bary.derived().resize(points.rows(), 3);

  int t;
  fmesh::Point s;
  fmesh::Point b;

  for (int i=0; i<points.rows(); i++) {
    s[0] = points(i,0);
    s[1] = points(i,1);
    s[2] = points(i,2);
    t = locator.locate(s);
    if (t>=0) { /* Point located. */
      M.barycentric(fmesh::Dart(M,t),s,b); /* Coordinates relative to
					      canonical vertex
					      ordering. */
      point2T(i,0) = t;
      point2bary(i,0) = b[0];
      point2bary(i,1) = b[1];
      point2bary(i,2) = b[2];
    } else { /* Point not found. */
      point2T(i,0) = -1;
    }
  }
}


typedef Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> MatrixX3Rd;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXXRd;

template <int DimCols=Eigen::Dynamic, typename T, typename Derived>
void map_fmesh_to_Eigen(fmesh::Matrix<T>& A,
			Eigen::MatrixBase<Derived> const & B_)
{
  Eigen::MatrixBase<Derived>& B = const_cast< Eigen::MatrixBase<Derived>& >(B_);
  //  B.derived().resize(A.rows(), A.rows());

  typedef Eigen::Matrix<T, Eigen::Dynamic, DimCols, Eigen::RowMajor> MatrixRT;
  B = Eigen::Map<const MatrixRT>(A.raw(), A.rows(), A.cols());
}

template <int DimCols=Eigen::Dynamic, typename T>
void map_fmesh_to_Eigen(fmesh::Matrix<T>& A,
			Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, DimCols, Eigen::RowMajor> >
			& B)
{
  typedef Eigen::Matrix<T, Eigen::Dynamic, DimCols, Eigen::RowMajor> MatrixRT;
  new (&B) Eigen::Map<const MatrixRT>(A.raw(), A.rows());
}


template <typename T, typename Derived>
void map_fmesh_to_Eigen(fmesh::SparseMatrix<T>& A,
			Eigen::SparseMatrixBase<Derived> const & B_)
{
  Eigen::SparseMatrixBase<Derived>& B = const_cast< Eigen::SparseMatrixBase<Derived>& >(B_);
  B.derived().resize(A.rows(), A.rows());
  B.derived().setZero();

  fmesh::Matrix1< fmesh::SparseMatrixTriplet<T> > Ttrip;
  int nnz = A.tolist(Ttrip);

  std::vector< Eigen::Triplet<double> > tripletList;
  tripletList.reserve(nnz);
  for (int i=0; i < nnz; ++i) {
    tripletList.push_back(Eigen::Triplet<double>(Ttrip(i).r, Ttrip(i).c, Ttrip(i).value));
  }
  B.derived().reserve(nnz);
  B.derived().setFromTriplets(tripletList.begin(), tripletList.end());
}





int main(int argc, char *argv[])
{
  Timer timertotal;
  TimerHierarchy timer(1);
  timertotal.tic();

  /*
   * Handle options
   */

  if (argc == 1) {
    std::cout << "Syntax:\t" << argv[0] << " nbThreads [subdivision_level=0]" << std::endl;
    exit(0);
  }

  if (argc > 1) {
    Eigen::setNbThreads(atol(argv[1]));
  }
  std::cout << "nbThreads:\t" << Eigen::nbThreads() << std::endl;

  int global_subdivision_level = 0;
  if (argc > 2) {
    global_subdivision_level = atol(argv[2]);
  }
  int globe_points_param = pow2(global_subdivision_level);
  std::cout << "global_subdivision_level:\t" << global_subdivision_level << std::endl;
  std::cout << "globe_points_param:\t" << globe_points_param << std::endl;

  /*
   * Construct global mesh
   */

  timer.tic("Construct global mesh");

  FMatrix* S0;
  fmesh::Mesh M(fmesh::Mesh::Mtype_sphere, 0, true, true);
  fmesh::MeshC MC(&M);
  MC.setOptions(MC.getOptions() | fmesh::MeshC::Option_offcenter_steiner);
  
  timer.tic("Make globe points");
  S0 = (FMatrix*)fmesh::make_globe_points(globe_points_param);
  timer.toc().tic("Attach points to mesh");
  M.S_append(*S0);
  delete S0;
  timer.toc().tic("Construct covering proto mesh");
  MC.CET(3, 1.0);
  timer.toc().tic("Add vertices");
  fmesh::vertexListT vertices;
  for (size_t v=0; v < M.S().rows(); v++) {
    vertices.push_back(v);
  }
  MC.DT(vertices);
  timer.toc();

  std::cout << "M.nV:\t" << M.nV() << ",\t" << "M.nT:\t" << M.nT() << std::endl;
  //  for (size_t v=0; v < 12; v++) {
  //    std::cout << "M.S(" << v << "):\t" << M.S(v) << std::endl;
  //  }

  timer.toc(-1);

  /*
   * Locate points
   */

  timer.tic("Locate points");

  Eigen::Matrix<double, Eigen::Dynamic, 3> points(M.nV(), 3);
  Eigen::Matrix<int, Eigen::Dynamic, 1> triangles(points.rows());
  Eigen::Matrix<double, Eigen::Dynamic, 3> barycentric(points.rows(), 3);

  timer.tic("Make points");
  for (int v=0; v < points.rows(); v++) {
    points(v,0) = M.S(v)[0];
    points(v,1) = M.S(v)[1];
    points(v,2) = M.S(v)[2];
  }
  timer.toc();
    
  timer.tic("Build locator");
  fmesh::TriangleLocator* locator = initialise_locator_local(M);
  timer.toc().tic("Map points");
  map_points_to_mesh(M, *locator, points, triangles, barycentric);
  timer.toc();

  //  std::cout << "Triangles:\n" << triangles << std::endl;
  //  std::cout << "Barycentric:\n" << barycentric << std::endl;

  timer.tic("Build A matrix");
  std::vector< Eigen::Triplet<double> > tripletList;
  tripletList.reserve(points.rows() * 3);
  for (int i=0; i < points.rows(); ++i) {
    for (int j=0; j < 3; ++j) {
      tripletList.push_back(Eigen::Triplet<double>(i,
						   M.TV(triangles(i))[j],
						   barycentric(i, j)));
    }
  }
  Eigen::SparseMatrix<double> Amatrix(points.rows(), M.nV());
  Amatrix.setFromTriplets(tripletList.begin(), tripletList.end());
  timer.toc();

  timer.toc(-1);
  
  /*
   * Build FEM matrices
   */

  timer.tic("Construct FEM matrices in fmesher and convert to Eigen");
  timer.tic("Construct FEM matrices");

  fmesh::SparseMatrix<double> C0;
  fmesh::SparseMatrix<double> C1;
  fmesh::SparseMatrix<double> G1;
  fmesh::SparseMatrix<double> B1;
  fmesh::Matrix<double> Tareas;

  timer.tic("Basic matrices");
  M.calcQblocks(C0,C1,G1,B1,Tareas);
  timer.toc().tic("Higher order matrices");
  timer.tic("C0inv");
  fmesh::SparseMatrix<double> C0inv = inverse(C0,true);
  timer.toc().tic("tmp = C0inv * G1");
  fmesh::SparseMatrix<double> tmp = C0inv * G1;
  timer.toc().tic("G2 = G1 * tmp");
  fmesh::SparseMatrix<double> G2 = G1 * tmp;
  timer.toc();
  timer.toc();

  std::cout << "nnz(G1):\t" << G1.nnz() << ",\t" << "nnz(G2):\t" << G2.nnz() << std::endl;

  timer.toc();

  /*
   * Map matrices to Eigen
   */

  timer.tic("Map matrices to Eigen");

  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> MS;
  map_fmesh_to_Eigen(M.S(), MS);
  //  Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> > MapMS(NULL, M.nV());
  //  map_fmesh_to_Eigen(M.S(), MapMS);
  Eigen::Map<const MatrixX3Rd> MS2(M.S().raw(), M.nV(), 3);

  //  std::cout << "M.S()\n" << M.S() << std::endl;
  //  std::cout << "MS\n" << MS << std::endl;
  //  std::cout << "MS2\n" << MS2 << std::endl;
  std::cout << "Changing M.S(2,1)" << std::endl;
  double savevalue = M.S()(2,1);
  M.S()(2,1) = 10.0;
  //  std::cout << "M.S()\n" << M.S() << std::endl;
  //  std::cout << "MS\n" << MS << std::endl;
  //  std::cout << "MS2\n" << MS2 << std::endl;
  M.S()(2,1) = savevalue;

  timer.toc();

  /*
   * Sparse Eigen matrices
   */
/*
  timer.tic("Convert fmesh Eigen");

  Eigen::SparseMatrix<double> EC0;// = C0;
  Eigen::SparseMatrix<double> EC0inv;// = C0inv;
  Eigen::SparseMatrix<double> EG1;// = G1;
  Eigen::SparseMatrix<double> EG2;// = G2;

  timer.tic("Direct conversion");

  timer.tic("EC0");
  map_fmesh_to_Eigen(C0, EC0);
  timer.toc().tic("EC0inv");
  map_fmesh_to_Eigen(C0inv, EC0inv);
  timer.toc().tic("EG1");
  map_fmesh_to_Eigen(G1, EG1);
  timer.toc().tic("EG2");
  map_fmesh_to_Eigen(G2, EG2);
  timer.toc();

  timer.toc();

  std::cout << " G1(0,0):\t" <<  G1(0,0) << std::endl;
  std::cout << "EG1(0,0):\t" << EG1.coeffRef(0,0) << std::endl;
  //  std::cout << " G1:\n" << G1 << std::endl;
  //  std::cout << "EG1:\n" << EG1.topLeftCorner(12,12) << std::endl;

  timer.tic("Higher order matrices");

  timer.tic("Etmp = EC0inv * EG1");
  Eigen::SparseMatrix<double> Etmp = EC0inv * EG1;
  timer.toc().tic("EG2b = EG1 * Etmp");
  Eigen::SparseMatrix<double> EG2b = EG1 * Etmp;
  timer.toc();

  timer.toc();
  timer.toc(-1);
*/

  /*
   * Direct construction of FEM matrices with Eigen
   */

  
  timer.tic("Direct construction of FEM matrices with Eigen");
  
  Eigen::Matrix<double, Eigen::Dynamic, 1> EEDC0;
  Eigen::Matrix<double, Eigen::Dynamic, 1> EEDC0inv;
  Eigen::SparseMatrix<double> EEC1;
  Eigen::SparseMatrix<double> EEG1;
  Eigen::SparseMatrix<double> EEB1;
  Eigen::Matrix<double, Eigen::Dynamic, 1> EEVareas;
  Eigen::Matrix<double, Eigen::Dynamic, 1> EETareas;
  
  timer.tic("Basic matrices");
  calcQblocks(M, EEDC0, EEDC0inv, EEC1, EEG1, EEB1, EETareas);
  timer.toc();

  timer.tic("Higher order matrices");
  timer.toc().tic("EEtmp = EEDC0inv.asDiag * EEG1");
  Eigen::SparseMatrix<double> EEtmp = EEDC0inv.asDiagonal() * EEG1;
  timer.toc().tic("EEG2b = EEG1 * EEtmp");
  Eigen::SparseMatrix<double> EEG2b = EEG1 * EEtmp;
  timer.toc();

  timer.toc(-1);


  /*
   * Check that the indirect and direct Eigen matrices match.
   */
/*
  std::cout << "|EC0.diag-EEDC0|:\t" << (EC0.diagonal()-EEDC0).norm() << std::endl;
  std::cout << "|EG1-EEG1|:\t" << (EG1-EEG1).norm() << std::endl;
  std::cout << "|EG2-EEG2b|:\t" << (EG2-EEG2b).norm() << std::endl;
*/




  Eigen::SparseMatrix<double> Q = EEG2b;
  Q.diagonal() += EEDC0;

  /*
   * Sample from a model
   */

  /*
  timer.tic("Sample from a model");

  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > Q_LLT;
  timer.tic("analyzePattern(Q)");
  Q_LLT.analyzePattern(Q);
  timer.toc().tic("factorize(Q)");
  Q_LLT.factorize(Q);
  timer.toc();
  //  std::cout << "Q_LLT:\n" << Q_LLT << std::endl;

  timer.tic("Extract Cholesky factor");
  Eigen::SparseMatrix<double> Q_L = Q_LLT.matrixL();
  timer.toc();
  std::cout << "Q_L:\n" << Q_L.topLeftCorner(12,12) << std::endl;

  timer.tic("makeCompressed(Q)");
  Q.makeCompressed();
  timer.toc().tic("analyzePattern(Q)");
  Q_LLT.analyzePattern(Q);
  std::cout << "Success = " << Q_LLT.info() << std::endl;
  timer.toc().tic("factorize(Q)");
  Q_LLT.factorize(Q);
  std::cout << "Success = " << Q_LLT.info() << std::endl;
  timer.toc();
  */
  
  Eigen::Matrix<double, Eigen::Dynamic, 1> u1, u2;
  Eigen::Matrix<double, Eigen::Dynamic, 1> b(Q.rows());
  timer.tic("b = Random");
  b.Random(Q.rows());
  timer.toc();
  /*
  timer.tic("u1 = LLT.solve(b)");
  u1 = Q_LLT.solve(b);
  timer.toc();

  // P*Q*P' = LL'
  // LL'*P'*u = P*Q*P'*u = P*b
  // u = P'*solve(L', solve(L, P*b)
  timer.tic("u2 = P' solve(L', solve(L, P b))");
  u2 = Q_LLT.permutationP().transpose() *
    Q_LLT.matrixU().solve(Q_LLT.matrixL().solve(Q_LLT.permutationP() * b));
  timer.toc();
  std::cout << "|u2-u1|:\t" << (u2-u1).norm() << std::endl;

  timer.toc();
  */


  timer.tic("Qtool (AMD)");
  QTool<double, Eigen::AMDOrdering> Q_tool;
  timer.tic("Assign Q");
  Q_tool.Q(Q);
  timer.toc().print().tic("analyzePattern");
  Q_tool.analyzePattern();
  timer.toc().print().tic("factorize");
  Q_tool.factorize();
  timer.toc().print().tic("partialInversion");
  Q_tool.partialInversion();
  timer.toc().print();

  timer.tic("solveQ(b)");
  Q_tool.solve(b, u1);
  timer.toc();

  timer.tic("solveLtP(solvePtL(b)");
  Eigen::Matrix<double, Eigen::Dynamic, 1> u3;
  Q_tool.solvePtL(b, u3);
  Q_tool.solveLtP(u3, u2);
  timer.toc();
  std::cout << "|u2-u1|:\t" << (u2-u1).norm() << std::endl;

  timer.toc();

  for (int looping=0; looping < 2; ++looping) {
    if (looping % 2 == 0) {
      timer.tic("Rtest totals");
    } else {
      timer.tic("Rtest parts");
    }
    Q.makeCompressed();
    Q_tool.Q(Q);
    timer.tic("reordering");
    Q_tool.analyzePattern();
    timer.toc();
    if (looping % 2 == 0) Q_tool.Q(Q);
    timer.tic("factorize");
    Q_tool.factorize();
    std::cout << "nonZeros(L)/n = "
	      << double(Q_tool.LLt().matrixL().eval().nonZeros()) / Q.rows()
	      << std::endl;
    timer.toc();
    if (looping % 2 == 0) Q_tool.Q(Q);
    timer.tic("partialInversion");
    Q_tool.partialInversion();
    timer.toc();
    if (looping % 2 == 0) Q_tool.Q(Q);
    timer.tic("solve");
    Q_tool.solve(b, u1);
    timer.toc();
    timer.toc();
  }
  
  timer.toc(-1);

  /*
   * End
   */

  timer.toc(-1);
  timertotal.toc();

  std::cout << "\nTimings:" << std::endl;
  timer.print();
  timertotal.print("Total:\t");

  time_t tt;
  tt = std::chrono::high_resolution_clock::to_time_t ( timertotal.start() );
  std::cout << "Start:  " << ctime(&tt);
  tt = std::chrono::high_resolution_clock::to_time_t ( timertotal.ending() );
  std::cout << "Ending: " << ctime(&tt);

  exit(0);
}
