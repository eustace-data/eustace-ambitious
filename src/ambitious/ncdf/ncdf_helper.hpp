#ifndef NCDF_HELPER_HPP
#define NCDF_HELPER_HPP

#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(ncdf_helper_hpp, "$Revision: 1297 $")

#include "ambitious/debuglog/debuglog.hpp"

#include <netcdf.h>
#ifdef USE_PNETCDF
#include <pnetcdf.h>
#endif
#include <map>
#include <vector>
#include <string>

/*
 * http://cucis.ece.northwestern.edu/projects/PnetCDF/doc/pnetcdf-c/ncmpi_005fwait_002fwait_005fall.html#ncmpi_005fwait_002fwait_005fall
 */

// Make a typedef available even when MPI isn't available.
#ifdef USE_PNETCDF
typedef MPI_Datatype MPI_Datatype_nc;
typedef MPI_Offset MPI_Offset_nc;
typedef MPI_Comm MPI_Comm_nc;
typedef MPI_Info MPI_Info_nc;
#else
typedef nc_type MPI_Datatype_nc;
typedef size_t MPI_Offset_nc;
typedef void* MPI_Comm_nc;
typedef void* MPI_Info_nc;
#endif
typedef std::vector<MPI_Offset_nc> MPI_Offset_vector;

#define LOG_NC_ERR(e) ELOG_("Error: " << std::string(nc_strerror(e)) << std::endl)
#ifdef USE_PNETCDF
#define LOG_PNC_ERR(e) ELOG_("Error: " << std::string(ncmpi_strerror(e)) << std::endl)
#else
#define LOG_PNC_ERR(e) ELOG_("Error: pnetcdf is not enabled." << std::endl)
#endif
/** Handle netcdf errors by printing an error message and exiting with a
  * non-zero status. For internal use only. */
#define NC_ERR(e) { if (e != NC_NOERR) { LOG_NC_ERR(e); exit(2); } }
/** Handle pnetcdf errors by printing an error message and exiting with a
  * non-zero status.  For internal use only. */
#define PNC_ERR(e) { if (e != NC_NOERR) { LOG_PNC_ERR(e); exit(2); } }
/** Handle netcdf and pnetcdf errors by printing an error message and exiting with a
  * non-zero status.  For internal use only. */
#define NC_ERROR(mpi,e) { if (mpi) { PNC_ERR(e); } else { NC_ERR(e); } }

/** Handle NCHelper errors by printing an error message and exiting with a
  * non-zero status. For internal use only. */
#define NCERR(h) { h.err(); }
/** Handle NCHelper errors from NCHelper::wait and NCHelper::wait_all
 * by printing an error message and exiting with a non-zero status.
 * For internal use only.
 */
#define NCSTAT(h) { h.stat(); }

/** Run netcdf or pnetcdf commands with automatic switching.
 *
 * Use NCDF(-1, ...) or NCDF(..., -1) to set error status when no
 * command is available and this is an error. Replace -1 with NC_NOERR
 * if this is not an error.
 *
 * For internal use only.
 * @internal
 */
#ifdef USE_PNETCDF
#define NCDF(CCC, CCCMPI) { if (mpi_) { ret_ = CCCMPI; } else { ret_ = CCC; } }
#else
#define NCDF(CCC, CCCMPI) { if (mpi_) { ELOG_("Error: Can't run ncmpi commands without pnetcdf."); ret_ = -1; } else { ret_ = CCC; } }
#endif

/** For internal use only. */
#define NCIFERR( action ) { if (ret_ != NC_NOERR) { action; } }
/** For internal use only. */
#define NCNOTERR( action ) { if (ret_ == NC_NOERR) { action; } }



#ifndef USE_PNETCDF
#else
MPI_Datatype_nc nc_type_to_mpi(nc_type type);
#endif

class NCHelper {
 public:
  class VarInfo {
   public:
    nc_type type;                   /* variable type */
    MPI_Datatype_nc type_mpi;       /* variable type */
    int ndims;                      /* number of dims */
    int dimids[NC_MAX_VAR_DIMS];    /* dimension ids */
    MPI_Offset_vector dimlen;       /* dimension lengths */
    int natts;                      /* number of attributes */
    MPI_Offset_nc size;             /* total number of elements */
  };
 protected:
  int ncid_;
  int ret_;
  std::vector<int> requests_;
  std::vector<int> statuses_;
  typedef std::map<std::string, int> id_map_t;
  typedef typename id_map_t::value_type map_t;
  typedef std::map<int, VarInfo> VarInfo_map_t;
  id_map_t dimid_map_;
  id_map_t varid_map_;
  VarInfo_map_t varinfo_map_;

  bool mpi_;
#ifdef USE_PNETCDF
  MPI_Comm mpi_comm_;
  MPI_Info mpi_info_;
#endif

 public:
  NCHelper() : ncid_(0), ret_(NC_NOERR), requests_(0), statuses_(0),
               dimid_map_(), varid_map_(),
               mpi_(false)
#ifdef USE_PNETCDF
             , mpi_comm_(NULL), mpi_info_(NULL)
#endif
  {
  }
  NCHelper(int ncid) : ncid_(0), ret_(NC_NOERR), requests_(0), statuses_(0),
                       dimid_map_(), varid_map_(),
                       mpi_(false)
#ifdef USE_PNETCDF
                     , mpi_comm_(NULL), mpi_info_(NULL)
#endif
  {
    ncid_ = ncid;
  }
  virtual ~NCHelper() {
    if (ncid_ != 0) {
      close();
    }
  }
  int id() const {
    return ncid_;
  }
  int ret() const {
    return ret_;
  }
  std::vector<int>& requests() {
    return requests_;
  }
  const std::vector<int>& requests() const {
    return requests_;
  }
  NCHelper& requests_clear() {
    requests_.clear();
    return *this;
  }
  /** Return the current error statuses from calls to wait() or wait_all()
   *
   * The #NCSTAT macro can also be used to check for error statuses.
   */
  std::vector<int>& statuses() {
    return statuses_;
  }
  const std::vector<int>& statuses() const {
    return statuses_;
  }
  NCHelper& statuses_clear() {
    statuses_.clear();
    return *this;
  }

  bool mpi() const {
    return mpi_;
  }

  const id_map_t& dimid_map() const {
    return dimid_map_;
  }
  const id_map_t& varid_map() const {
    return varid_map_;
  }
  const VarInfo_map_t& varinfo_map() const {
    return varinfo_map_;
  }
#ifdef USE_PNETCDF
  const MPI_Comm_nc& mpi_comm() const {
    return mpi_comm_;
  }
  const MPI_Info_nc& mpi_info() const {
    return mpi_info_;
  }
  MPI_Comm_nc& mpi_comm() {
    return mpi_comm_;
  }
  MPI_Info_nc& mpi_info() {
    return mpi_info_;
  }
#else
  const MPI_Comm_nc& mpi_comm() const {
    return NULL;
  }
  const MPI_Info_nc& mpi_info() const {
    return NULL;
  }
  MPI_Comm_nc& mpi_comm() {
    return NULL;
  }
  MPI_Info_nc& mpi_info() {
    return NULL;
  }
#endif



  /** Handle NCHelper errors by printing an error message and exiting with a
   * non-zero status.
   */
  void err() {
    NC_ERROR(mpi_, ret_);
  }
  /** Handle NCHelper errors from NCHelper::wait and NCHelper::wait_all
   * by printing an error message and exiting with a non-zero status.
   */
  void stat() {
    for (std::vector<int>::iterator it = statuses_.begin();
         it != statuses_.end();
         ++it) {
      NC_ERROR(mpi_, (*it));
    }
  }


  
#ifdef USE_PNETCDF
  NCHelper& set_mpi(MPI_Comm_nc comm, MPI_Info_nc info) {
    if (mpi_) {
      unset_mpi();
    }
    mpi_ = true;
    mpi_comm_ = comm;
    mpi_info_ = info;
    return *this;
  }
#else
  NCHelper& set_mpi(MPI_Comm_nc comm, MPI_Info_nc info) {
    mpi_ = false;
    return *this;
  }
#endif
  NCHelper& unset_mpi() {
    mpi_ = false;
    return *this;
  }
  int mpi_rank() {
    int rank = 0;
#ifdef USE_PNETCDF
    if (mpi_) {
      MPI_Comm_rank(mpi_comm_, &rank);
    }
#endif
    return rank;
  }
  int mpi_size() {
    int size = 1;
#ifdef USE_PNETCDF
    if (mpi_) {
      MPI_Comm_size(mpi_comm_, &size);
    }
#endif
    return size;
  }
  
  NCHelper& set_id(int ncid) {
    ncid_ = ncid;
    ret_ = NC_NOERR;
    return *this;
  }
  NCHelper& unset_id() {
    ncid_ = 0;
    ret_ = NC_NOERR;
    return *this;
  }
  NCHelper& create(const std::string& filename,
                   int cmode) {
    NCDF(nc_create(filename.c_str(), cmode, &ncid_),
         ncmpi_create(mpi_comm_, filename.c_str(), cmode, mpi_info_, &ncid_));
    LOG_("NCDF create: ncid = " << ncid_ << std::endl);
    return *this;
  }
  NCHelper& open(const std::string& filename,
                 int mode) {
    NCDF(nc_open(filename.c_str(), mode, &ncid_),
         ncmpi_open(mpi_comm_, filename.c_str(), mode, mpi_info_, &ncid_));
    return *this;
  }

  NCHelper& close() {
    LOG_("NCDF close: ncid = " << ncid_ << std::endl);
    NCDF(nc_close(ncid_),
         ncmpi_close(ncid_));
    ncid_ = 0;
    return *this;
  }

  NCHelper& clear() {
    if (!requests_.empty()) {
      ELOG_("Warning: There are " << requests_.size() << " nonblocking requests with unknown status." << std::endl);
    }
    if (ncid_) {
      close();
    }
    return *this;
  }

  int dimid(const std::string& name) {
    id_map_t::iterator i = dimid_map_.find(name);
    if (i != dimid_map_.end()) {
      LOG_("NCDF dimid: " << name << " = " << i->second << std::endl);
      return i->second;
    }
    int id = -1;
    NCDF(nc_inq_dimid(ncid_, name.c_str(), &id),
         ncmpi_inq_dimid(ncid_, name.c_str(), &id));
    NCNOTERR(dimid_map_.insert(map_t(name, id)));
    LOG_("NCDF dimid: " << name << " = " << id << std::endl);
    return id;
  }
  int varid(const std::string& name) {
    id_map_t::iterator i = varid_map_.find(name);
    if (i != varid_map_.end()) {
    LOG_("NCDF varid: " << name << " = " << i->second << std::endl);
      return i->second;
    }
    int id = -1;
    NCDF(nc_inq_varid(ncid_, name.c_str(), &id),
         ncmpi_inq_varid(ncid_, name.c_str(), &id));
    NCNOTERR(varid_map_.insert(map_t(name, id)));
    LOG_("NCDF varid: " << name << " = " << id << std::endl);
    return id;
  }

  const VarInfo* var_info(int varid) {
    VarInfo_map_t::iterator it = varinfo_map_.find(varid);
    if (it != varinfo_map_.end()) {
      return &(it->second);
    }
    // Not cached.
    
    VarInfo info;
    NCDF(nc_inq_var(ncid_, varid, NULL,
                    &info.type, &info.ndims, info.dimids, &info.natts),
         ncmpi_inq_var(ncid_, varid, NULL,
                       &info.type, &info.ndims, info.dimids, &info.natts));
    NCIFERR(return NULL);
#ifndef USE_PNETCDF
    info.type_mpi = 0;
#else
    info.type_mpi = MPI_DATATYPE_NULL;
    if (mpi_) {
      // Translate nc_type to MPI_Datatype
      info.type_mpi = nc_type_to_mpi(info.type);
    }
#endif

    info.dimlen.resize(info.ndims);
    size_t dimlen;
    MPI_Offset_nc dimlen_mpi; // netcdf version is size_t
    for (int i = 0; i < info.ndims; ++i) {
      NCDF(nc_inq_dimlen(ncid_, info.dimids[i], &dimlen),
           ncmpi_inq_dimlen(ncid_, info.dimids[i], &dimlen_mpi));
      NCIFERR(return NULL);
      if (mpi_) {
        info.dimlen[i] = dimlen_mpi;
      } else {
        info.dimlen[i] = dimlen;
      }
    }
    info.size = count_size(info.dimlen);

    std::pair<VarInfo_map_t::iterator, bool> result =
        varinfo_map_.insert(VarInfo_map_t::value_type(varid, info));

    return &(result.first->second);
  }

  std::vector<int> dimid(const std::vector<std::string>& name) {
    std::vector<int> ids;
    ids.resize(name.size());
    for (size_t i = 0; i < name.size(); ++i) {
      ids[i] = dimid(name[i]);
    }
    return ids;
  }
  std::vector<int> varid(const std::vector<std::string>& name) {
    std::vector<int> ids;
    ids.resize(name.size());
    for (size_t i = 0; i < name.size(); ++i) {
      ids[i] = varid(name[i]);
    }
    return ids;
  }

  int def_dim(const std::string& name, size_t len) {
    int id;
    NCDF(nc_def_dim(ncid_, name.c_str(), len, &id),
         ncmpi_def_dim(ncid_, name.c_str(), len, &id));
    NCNOTERR(dimid_map_.insert(map_t(name, id)));
    return id;
  }
  int def_var(const std::string& name,
              nc_type xtype,
              int dim) {
    int id;
    NCDF(nc_def_var(ncid_, name.c_str(), xtype, 1, &dim, &id),
         ncmpi_def_var(ncid_, name.c_str(), xtype, 1, &dim, &id));
    NCNOTERR(varid_map_.insert(map_t(name, id)));
    return id;
  }
  int def_var(const std::string& name,
              nc_type xtype,
              std::vector<int> dims) {
    int id;
    NCDF(nc_def_var(ncid_, name.c_str(), xtype, dims.size(), dims.data(), &id),
         ncmpi_def_var(ncid_, name.c_str(), xtype, dims.size(), dims.data(), &id));
    NCNOTERR(varid_map_.insert(map_t(name, id)));
    return id;
  }
  int def_var(const std::string& name,
              nc_type xtype,
              const std::string& dim) {
    return def_var(name, xtype, dimid(dim));
  }
  int def_var(const std::string& name,
              nc_type xtype,
              const std::vector<std::string>& dims) {
    return def_var(name, xtype, dimid(dims));
  }
  NCHelper& put_att_text(int varid,
                         const std::string& name,
                         const std::string& value) {
    NCDF(nc_put_att_text(ncid_, varid, name.c_str(), 
                         value.length(),
                         value.c_str()),
         ncmpi_put_att_text(ncid_, varid, name.c_str(), 
                            value.length(),
                            value.c_str()));
    return *this;
  }
  template <typename ValueT>
  NCHelper& put_att(int varid,
                    const std::string& name,
                    nc_type xtype,
                    const std::vector<ValueT>& value)
  {
    LOG_("value.size = " << value.size() << std::endl);
    NCDF(nc_put_att(ncid_, varid, name.c_str(), xtype,
                    value.size(), value.data()),
         ncmpi_put_att(ncid_, varid, name.c_str(), xtype,
                       value.size(), value.data()));
    return *this;
  }
  NCHelper& put_att_text(const std::string& varname,
                         const std::string& name,
                         const std::string& value) {
    return put_att_text(varid(varname), name, value);
  }
  template <typename ValueT>
  NCHelper& put_att(const std::string& varname,
                    const std::string& name,
                    nc_type xtype,
                    const std::vector<ValueT>& value) {
    return put_att(varid(varname), name, xtype, value);
  }

  NCHelper& enddef() {
    NCDF(nc_enddef(ncid_),
         ncmpi_enddef(ncid_));
    return *this;
  }


  uint64_t count_size(const MPI_Offset_vector& count) const {
    uint64_t countprod = 1;
    for (size_t k = 0; k < count.size(); ++k) {
      countprod *= count[k];
    }
    return countprod;
  }


  template <typename ValueT>
  NCHelper& put_var(int varid, const ValueT* op) {
    if (!mpi_) {
      NCDF(nc_put_var(ncid_, varid, op),
           -1);
      return *this;
    }
#ifndef USE_PNETCDF
    NCDF(-1, -1);
#else
    /* TODO in future:
     * From pnetcdf 1.6.0, setting bufcount,buftype to {anynumber},MPI_DATATYPE_NULL
     * ignores bufcount and requires matching data types between code and file,
     * just as for netcdf.
     */

    const VarInfo* varinfo = var_info(varid);
    NCIFERR(return *this);

    NCDF(-1,
         ncmpi_put_var(ncid_, varid, op,
                       varinfo->size,
                       varinfo->type_mpi));
    NCIFERR(return *this);
#endif
    return *this;
  }

  template <typename ValueT>
  NCHelper& get_var(int varid, ValueT* op) {
    if (!mpi_) {
      NCDF(nc_get_var(ncid_, varid, op),
           -1);
      return *this;
    }
#ifndef USE_PNETCDF
    NCDF(-1, -1);
#else
    /* TODO in future:
     * From pnetcdf 1.6.0, setting bufcount,buftype to {anynumber},MPI_DATATYPE_NULL
     * ignores bufcount and requires matching data types between code and file,
     * just as for netcdf.
     */

    const VarInfo* varinfo = var_info(varid);
    NCIFERR(return *this);

    NCDF(-1,
         ncmpi_get_var(ncid_, varid, op,
                       varinfo->size,
                       varinfo->type_mpi));
    NCIFERR(return *this);
#endif
    return *this;
  }

  template <typename ValueT>
  NCHelper& iput_var(int varid, const ValueT* op) {
    if (!mpi_) {
      NCDF(nc_put_var(ncid_, varid, op),
           -1);
      return *this;
    }
#ifndef USE_PNETCDF
    NCDF(-1, -1);
#else
    int new_request;
    /* TODO in future:
     * From pnetcdf 1.6.0, setting bufcount,buftype to {anynumber},MPI_DATATYPE_NULL
     * ignores bufcount and requires matching data types between code and file,
     * just as for netcdf.
     */

    const VarInfo* varinfo = var_info(varid);
    NCIFERR(return *this);

    NCDF(-1,
         ncmpi_iput_var(ncid_, varid, op,
                        varinfo->size,
                        varinfo->type_mpi,
                        &new_request));
    NCNOTERR(requests_.push_back(new_request));
#endif
    return *this;
  }

  template <typename ValueT>
  NCHelper& iget_var(int varid, ValueT* op) {
    if (!mpi_) {
      NCDF(nc_get_var(ncid_, varid, op),
           -1);
      return *this;
    }
#ifndef USE_PNETCDF
    NCDF(-1, -1);
#else
    int new_request;
    /* TODO in future:
     * From pnetcdf 1.6.0, setting bufcount,buftype to {anynumber},MPI_DATATYPE_NULL
     * ignores bufcount and requires matching data types between code and file,
     * just as for netcdf.
     */

    const VarInfo* varinfo = var_info(varid);
    NCIFERR(return *this);

    NCDF(-1,
         ncmpi_iget_var(ncid_, varid, op,
                        varinfo->size,
                        varinfo->type_mpi,
                        &new_request));
    NCNOTERR(requests_.push_back(new_request));
#endif
    return *this;
  }

  template <typename ValueT>
  NCHelper& iput_vara(int varid,
                      const MPI_Offset_vector& start,
                      const MPI_Offset_vector& count,
                      const ValueT* op) {
    // TODO: check that the size of start and count matches the #dim.
    if (!mpi_) {
      NCDF(nc_put_vara(ncid_, varid, start.data(), count.data(), op),
           -1);
      return *this;
    }
#ifndef USE_PNETCDF
    NCDF(-1, -1);
#else
    int new_request;
    /* TODO in future:
     * From pnetcdf 1.6.0, setting bufcount,buftype to {anynumber},MPI_DATATYPE_NULL
     * ignores bufcount and requires matching data types between code and file,
     * just as for netcdf.
     */

    const VarInfo* varinfo = var_info(varid);
    NCIFERR(return *this);

    // Check size
    assert((int64_t)count_size(count) <= varinfo->size);

    NCDF(-1,
         ncmpi_iput_vara(ncid_, varid,
                         start.data(),
                         count.data(),
                         op,
                         count_size(count),
                         varinfo->type_mpi,
                         &new_request));
    NCNOTERR(requests_.push_back(new_request));
#endif
    return *this;
  }

  template <typename ValueT>
  NCHelper& iget_vara(int varid,
                      const MPI_Offset_vector& start,
                      const MPI_Offset_vector& count,
                      ValueT* op) {
    // TODO: check that the size of start and count matches the #dim.
    if (!mpi_) {
      NCDF(nc_get_vara(ncid_, varid, start.data(), count.data(), op),
           -1);
      return *this;
    }
#ifndef USE_PNETCDF
    NCDF(-1, -1);
#else
    int new_request;
    /* TODO in future:
     * From pnetcdf 1.6.0, setting bufcount,buftype to {anynumber},MPI_DATATYPE_NULL
     * ignores bufcount and requires matching data types between code and file,
     * just as for netcdf.
     */

    const VarInfo* varinfo = var_info(varid);
    NCIFERR(return *this);

    // Check size
    assert((int64_t)count_size(count) <= varinfo->size);

    NCDF(-1,
         ncmpi_iget_vara(ncid_, varid,
                         start.data(),
                         count.data(),
                         op,
                         count_size(count),
                         varinfo->type_mpi,
                         &new_request));
    NCNOTERR(requests_.push_back(new_request));
#endif
    return *this;
  }

  template <typename ValueT>
  NCHelper& iput_varm(int varid,
                      const MPI_Offset_vector& start,
                      const MPI_Offset_vector& count,
                      const MPI_Offset_vector& stride,
                      const MPI_Offset_vector& imap,
                      const ValueT* op) {
    // TODO: check that the size of start and count matches the #dim.
    if (!mpi_) {
      std::vector<size_t> start_converted(start.begin(), start.end());
      std::vector<size_t> count_converted(count.begin(), count.end());
      std::vector<ptrdiff_t> stride_converted(stride.begin(), stride.end());
      std::vector<ptrdiff_t> imap_converted(imap.begin(), imap.end());
      NCDF(nc_put_varm(ncid_,
                       varid,
                       start_converted.data(),
                       count_converted.data(),
                       stride_converted.data(),
                       imap_converted.data(),
                       op),
           -1);
      return *this;
    }
#ifndef USE_PNETCDF
    NCDF(-1, -1);
#else
    int new_request;
    /* TODO in future:
     * From pnetcdf 1.6.0, setting bufcount,buftype to {anynumber},MPI_DATATYPE_NULL
     * ignores bufcount and requires matching data types between code and file,
     * just as for netcdf.
     */

    const VarInfo* varinfo = var_info(varid);
    NCIFERR(return *this);

    // Check size
    assert((int64_t)count_size(count) <= varinfo->size);

    NCDF(-1,
         ncmpi_iput_varm(ncid_, varid,
                         start.data(),
                         count.data(),
                         stride.data(),
                         imap.data(),
                         op,
                         count_size(count),
                         varinfo->type_mpi,
                         &new_request));
    NCNOTERR(requests_.push_back(new_request));
#endif
    return *this;
  }

  template <typename ValueT>
  NCHelper& iget_varm(int varid,
                      const MPI_Offset_vector& start,
                      const MPI_Offset_vector& count,
                      const MPI_Offset_vector& stride,
                      const MPI_Offset_vector& imap,
                      ValueT* op) {
    // TODO: check that the size of start and count matches the #dim.
    if (!mpi_) {
      std::vector<size_t> start_converted(start.begin(), start.end());
      std::vector<size_t> count_converted(count.begin(), count.end());
      std::vector<ptrdiff_t> stride_converted(stride.begin(), stride.end());
      std::vector<ptrdiff_t> imap_converted(imap.begin(), imap.end());
      NCDF(nc_get_varm(ncid_,
                       varid,
                       start_converted.data(),
                       count_converted.data(),
                       stride_converted.data(),
                       imap_converted.data(),
                       op),
           -1);
      return *this;
    }
#ifndef USE_PNETCDF
    NCDF(-1, -1);
#else
    int new_request;
    /* TODO in future:
     * From pnetcdf 1.6.0, setting bufcount,buftype to {anynumber},MPI_DATATYPE_NULL
     * ignores bufcount and requires matching data types between code and file,
     * just as for netcdf.
     */

    const VarInfo* varinfo = var_info(varid);
    NCIFERR(return *this);

    // Check size
    assert((int64_t)count_size(count) <= varinfo->size);

    NCDF(-1,
         ncmpi_iget_varm(ncid_, varid,
                         start.data(),
                         count.data(),
                         stride.data(),
                         imap.data(),
                         op,
                         count_size(count),
                         varinfo->type_mpi,
                         &new_request));
    NCNOTERR(requests_.push_back(new_request));
#endif
    return *this;
  }

  /**
   * The stat() function can be used to check for error statuses. Direct access
   * to the statuses is through the statuses() method.
   */
  NCHelper& wait() {
    statuses_.resize(requests_.size());
    NCDF(NC_NOERR,
         ncmpi_wait(ncid_, requests_.size(), requests_.data(), statuses_.data()));
    requests_.clear();
    return *this;
  }
  /**
   * The stat() function can be used to check for error statuses. Direct access
   * to the statuses is through the statuses() method.
   */
  NCHelper& wait_all() {
    statuses_.resize(requests_.size());
    if (mpi_) {
      LOG_("MPI rank = " << mpi_rank() << ", waiting for " << requests_.size() << " requests." << std::endl);
    }
    NCDF(NC_NOERR,
         ncmpi_wait_all(ncid_, requests_.size(), requests_.data(), statuses_.data()));
    requests_.clear();
    return *this;
  }
  /* Available from 1.7.0 */
  /*
  NCHelper& wait(int num) {
    std::assert((num == NC_GET_REQ_ALL) || (num == NC_PUT_REQ_ALL));
    NCDF(NC_NOERR,
         ncmpi_wait(ncid_, num, NULL, NULL));
    return *this;
  }
  */

  

  const VarInfo* var_info(const std::string& name) {
    return var_info(varid(name));
  }

  template <typename ValueT>
  NCHelper& put_var(const std::string& name, const ValueT* op) {
    return put_var(varid(name), op);
  }
  template <typename ValueT>
  NCHelper& get_var(const std::string& name, ValueT* op) {
    return put_var(varid(name), op);
  }
  template <typename ValueT>
  NCHelper& iput_var(const std::string& name, const ValueT* op) {
    return iput_var(varid(name), op);
  }
  template <typename ValueT>
  NCHelper& iget_var(const std::string& name, ValueT* op) {
    return iget_var(varid(name), op);
  }
  template <typename ValueT>
  NCHelper& iput_vara(const std::string& name,
                      const MPI_Offset_vector& start,
                      const MPI_Offset_vector& count,
                      const ValueT* op) {
    return iput_vara(varid(name), start, count, op);
  }
  template <typename ValueT>
  NCHelper& iget_vara(const std::string& name,
                      const MPI_Offset_vector& start,
                      const MPI_Offset_vector& count,
                      ValueT* op) {
    return iget_vara(varid(name), start, count, op);
  }
  template <typename ValueT>
  NCHelper& iput_varm(const std::string& name,
                      const MPI_Offset_vector& start,
                      const MPI_Offset_vector& count,
                      const MPI_Offset_vector& stride,
                      const MPI_Offset_vector& imap,
                      const ValueT* op) {
    return iput_varm(varid(name), start, count, stride, imap, op);
  }
  template <typename ValueT>
  NCHelper& iget_varm(const std::string& name,
                      const MPI_Offset_vector& start,
                      const MPI_Offset_vector& count,
                      const MPI_Offset_vector& stride,
                      const MPI_Offset_vector& imap,
                      ValueT* op) {
    return iget_varm(varid(name), start, count, stride, imap, op);
  }
};




#endif
