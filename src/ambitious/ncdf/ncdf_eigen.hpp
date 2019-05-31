/** Read/write/update netcdf files from/to Eigen objects
 * @file
 *
 * * Construct scaling information for compressed storage
 * * Convert to/from compressed storage
 * * Read/write header information
 * * 
 *
 */

#ifndef NCDF_EIGEN_HPP
#define NCDF_EIGEN_HPP

#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(ncdf_eigen_hpp, "$Revision: 1297 $")

#include "ambitious/debuglog/debuglog.hpp"
#include "ambitious/ncdf/ncdf_helper.hpp"
#include <Eigen/Dense>

#include <map>
#include <vector>
#include <string>

/** Compressed vector for netcdf
 *
 */
template <class OT, class CT>
class CompressedVector {
 public:
  typedef OT value_type;
  typedef CT compressed_type;
  typedef Eigen::Matrix<OT, Eigen::Dynamic, 1> vector_type;
  typedef Eigen::Matrix<CT, Eigen::Dynamic, 1> compressed_vector_type;

  class meta_type {
   public:
    OT offset_;
    OT scale_;
    bool allow_missing_;
    CT missing_value_;

    meta_type()
        : offset_(),
          scale_(),
          allow_missing_(false),
          missing_value_(std::numeric_limits<CT>::max())
    {
    }
    meta_type(OT offset, OT scale, bool allow_missing, CT missing)
        : offset_(offset),
          scale_(scale),
          allow_missing_(allow_missing),
          missing_value_(missing)
    {
    }
  };

 protected:
  meta_type meta_;
  compressed_vector_type comp_vec_;
 public:
  CompressedVector()
      : meta_()
  {
  }

  const meta_type& meta() {
    return meta_;
  }

  /** Estimate the safe value limits on the original scale */
  std::pair<OT, OT> limits_estimate() const {
    OT scale = meta_.scale_ - meta_.scale_ * (std::numeric_limits<OT>::epsilon() * 8);
    if (std::numeric_limits<CT>::is_signed) { // Signed storage
      return std::pair<OT, OT>(meta_.offset_ + scale * (-std::numeric_limits<CT>::max()),
                               meta_.offset_ + scale * (+std::numeric_limits<CT>::max()));
    } else if (meta_.allow_missing_) {
      return std::pair<OT, OT>(meta_.offset_ + scale * (0),
                               meta_.offset_ + scale * (std::numeric_limits<CT>::max()-1));
    } else {
      return std::pair<OT, OT>(meta_.offset_ + scale * (0),
                               meta_.offset_ + scale * (+std::numeric_limits<CT>::max()));
    }
  }

  /** Estimate the safe value limits on the compressed scale */
  std::pair<CT, CT> comp_limits_estimate() const {
    std::pair<OT, OT> lim = limits_estimate();
    return std::pair<CT, CT>(compress(lim.first), compress(lim.second));
  }

  compressed_vector_type& comp_vec() {
    return comp_vec_;
  }

  CT operator[](size_t index) const {
    return comp_vec_[index];
  }

  CompressedVector& set_scaling(OT offset, OT scale, bool allow_missing, CT missing) {
    meta_ = meta_type(offset, scale, allow_missing, missing);
    return *this;
  }

  /**
   *
   * https://www.unidata.ucar.edu/blogs/developer/entry/compression_by_scaling_and_offfset
   *
   * Reserve the maximal CT value for missing values.
   *
   * @param lim The minimum and maximum storable data values
   */
  CompressedVector& set_scaling(const std::pair<OT, OT>& lim, bool allow_missing) {
    uint64_t span;
    if (std::numeric_limits<CT>::is_signed) { // Signed storage
      span = 1 + 2 * (uint64_t)std::numeric_limits<CT>::max();
    } else { // Unsigned storage
      span = std::numeric_limits<CT>::max();
    }
    meta_.allow_missing_ = allow_missing;
    if (meta_.allow_missing_ || std::numeric_limits<CT>::is_signed) {
      meta_.scale_ = (lim.second - lim.first) / (span - 1);
    } else {
      meta_.scale_ = (lim.second - lim.first) / span;
    }
    if (std::numeric_limits<CT>::is_signed) { // Signed storage
      meta_.missing_value_ = std::numeric_limits<CT>::min();
      // (OT)midpoint should end up at (CT)0
      // (CT)0 = round[(midpoint-offset)/scale]
      // offset = midpoint
      meta_.offset_ = (lim.first + lim.second) / 2;
    } else {
      meta_.missing_value_ = std::numeric_limits<CT>::max();
      meta_.offset_ = lim.first;
    }
    /*
     * The 8*epsilon scaling is motivated by an error analysis bounding the
     * relative approximation error of (val-offset)/scale to 4*espilon in
     * the signed case and 2*epsilon in the unsigned case.
     */
    meta_.scale_ += meta_.scale_ * (std::numeric_limits<OT>::epsilon() * 8);
    return *this;
  }

  CT compress(OT val) const {
    return std::round((val - meta_.offset_) / meta_.scale_);
  }
  CompressedVector& compress(const vector_type& vec) {
    comp_vec_ = ((vec.array() - meta_.offset_) / meta_.scale_).round().template cast<CT>();
    return *this;
  }
  OT decompress(CT comp) const {
    return meta_.scale_ * (OT)comp + meta_.offset_;
  }
  CompressedVector& decompress(vector_type& vec) {
    vec = meta_.scale_ * comp_vec_.template cast<OT>().array() + meta_.offset_;
    return *this;
  }
  
};





class NCEigenHelper : public NCHelper {
 public:
  NCEigenHelper()
      : NCHelper()
  {
  }
  NCEigenHelper(int ncid)
      : NCHelper(ncid)
  {
  }
  NCEigenHelper(NCHelper& nch, bool clear_requests = true)
      : NCHelper(nch.id())
  {
    ret_ = nch.ret();
    requests_ = nch.requests();
    statuses_ = nch.statuses();
    dimid_map_ = nch.dimid_map();
    varid_map_ = nch.varid_map();
    varinfo_map_ = nch.varinfo_map();
    mpi_ = nch.mpi();
#ifdef USE_PNETCDF
    mpi_comm_ = nch.mpi_comm();
    mpi_info_ = nch.mpi_info();
#endif

    if (clear_requests) {
      nch.requests_clear();
      nch.statuses_clear();
    }
  }
  
  template <class Varid, class DerivedA>
  NCEigenHelper& put_var(const Varid& varid,
                          const Eigen::DenseBase<DerivedA>& var) {
    put_var(varid, var.data());
    return *this;
  }

  template <class Varid, class DerivedA>
  NCEigenHelper& get_var(const Varid& varid,
                         Eigen::DenseBase<DerivedA> const & var) {
    const VarInfo* info = var_info(varid);
    if (ret_ != NC_NOERR) { return *this; }
    var.derived().resize(info->size);
    get_var(varid, var.data());
    return *this;
  }

  template <class Varid, class DerivedA>
  NCEigenHelper& iput_var(const Varid& varid,
                          const Eigen::DenseBase<DerivedA>& var) {
    iput_var(varid, var.data());
    return *this;
  }

  template <class Varid, class DerivedA>
  NCEigenHelper& iput_vara(const Varid& varid,
                           const MPI_Offset_vector& start,
                           const MPI_Offset_vector& count,
                           const Eigen::DenseBase<DerivedA>& var) {
    iput_vara(varid, start, count, var.data());
    return *this;
  }

  template <class Varid, class DerivedA>
  NCEigenHelper& iget_vara(const Varid& varid,
                           const MPI_Offset_vector& start,
                           const MPI_Offset_vector& count,
                           Eigen::DenseBase<DerivedA> const & var) {
    var.derived().resize(count_size(count));
    iget_vara(varid, start, count, var.data());
    return *this;
  }

  template <class Varid, class DerivedA>
  NCEigenHelper& iput_varm(const Varid& varid,
                           const MPI_Offset_vector& start,
                           const MPI_Offset_vector& count,
                           const MPI_Offset_vector& stride,
                           const MPI_Offset_vector& imap,
                           const Eigen::DenseBase<DerivedA>& var) {
    iput_varm(varid, start, count, stride, imap, var.data());
    return *this;
  }

  template <class Varid, class DerivedA>
  NCEigenHelper& iget_varm(const Varid& varid,
                           const MPI_Offset_vector& start,
                           const MPI_Offset_vector& count,
                           const MPI_Offset_vector& stride,
                           const MPI_Offset_vector& imap,
                           Eigen::DenseBase<DerivedA> const & var) {
    iget_varm(varid, start, count, stride, imap, var.data());
    return *this;
  }

};

#endif
