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
  for (int t = 0; t < (int)M.nT(); t++) {
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
