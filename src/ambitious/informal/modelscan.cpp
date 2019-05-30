#include "ambitious/debuglog/debuglog.hpp"
#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(modelscan_cpp, "$Revision: 1299 $")

#include <iostream>
#include "ambitious/timer/timer.hpp"
#include "ambitious/to_string/to_string.hpp"
#include "ambitious/ostreamer/ostreamer.hpp"
#include <string.h>
#include "ambitious/po/po.hpp"
#include "ambitious/observe/observe.hpp"
#include <eustace/analysis/inputmanager.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include "ambitious/blocks/blocks.hpp"

using namespace EUSTACE;


/** Container for meshes and mesh mappers.
 */
class ModelMeshes {
 public:
  MeshBilevel* space_domain = NULL; /// Memory not handled by ModelMeshes
  TimeMesh* daily = NULL;
  TimeMesh* climate = NULL;
  MeshBilevel* space = NULL;
  TimeMesh* season = NULL;
  TimeMapper* daily_map = NULL;
  TimeMapper* climate_map = NULL;
  SpaceTimeMapper* subdomain_map = NULL; /// Mapper to the spatial subdomain.
  SpaceMapper* space_point_map = NULL;
  SpaceMapper* space_cell_map = NULL;
  TimeMapper* season_map = NULL;
  MPI_Group subdomain_group;
  MPI_Comm subdomain_comm;
 protected:
  int64_t epoch_calyear_;
  TimerHierarchy* timer_;
  TIMER_HELPERS
 public:

  ~ModelMeshes() {
    clear();
  }
  void clear() {
    if (daily_map) { delete daily_map; daily_map = NULL; }
    if (climate_map) { delete climate_map; climate_map = NULL; }
    if (subdomain_map) { delete subdomain_map; subdomain_map = NULL; }
    if (space_point_map) { delete space_point_map; space_point_map = NULL; }
    if (space_cell_map) { delete space_cell_map; space_cell_map = NULL; }
    if (season_map) { delete season_map; season_map = NULL; }
    if (daily) { delete daily; daily = NULL; }
    if (climate) { delete climate; climate = NULL; }
    if (space) { delete space; space = NULL; }
    if (season) { delete season; season = NULL; }
    LOG_(std::endl);
    if (subdomain_comm != MPI_COMM_NULL)
      MPI_Comm_free(&subdomain_comm);
    LOG_(std::endl);
    if (subdomain_group != MPI_GROUP_NULL)
      MPI_Group_free(&subdomain_group);
    LOG_(std::endl);
  }
  ModelMeshes(MeshBilevel* the_domain,
              const DaySpan& model_span,
              int64_t epoch_calyear,
              MPI_Comm& world_comm,
              TimerHierarchy* timer)
      : epoch_calyear_(epoch_calyear),
        timer_(timer)
  {
    space_domain = the_domain;
    LOG_("space_domain->macro_boundary_.rows = " << space_domain->macro_boundary_.rows() << std::endl);
    
    space = new MeshBilevel(timer_, true);
    std::set<int64_t> macro_vertices;

    LOG_(std::endl);

    LOG_("Number of domain roots: " << space_domain->macro_domain_root_set().size() << std::endl);
    int mpi_size;
    MPI_Comm_size(world_comm, &mpi_size);
    LOG_("Number of mpi processes: " << mpi_size << std::endl);
    assert(space_domain->macro_domain_root_set().size() <= (uint64_t)mpi_size);

    std::vector<int> subdomain_group_ranks(space_domain->macro_domain_root_set().size());
    for (uint k = 0; k < subdomain_group_ranks.size(); ++k) {
      subdomain_group_ranks[k] = k;
    }
    MPI_Group world_group;
    MPI_Comm_group(world_comm, &world_group);
    MPI_Group_incl(world_group, subdomain_group_ranks.size(), subdomain_group_ranks.data(), &subdomain_group);
    MPI_Comm_create_group(world_comm, subdomain_group, 0, &subdomain_comm);
    MPI_Group_free(&world_group);

    if (MPI_COMM_NULL == subdomain_comm) {
      // This process will not be involved in the computations.
      clear();
      return;
    }

    int my_mpi_rank;
    MPI_Comm_rank(subdomain_comm, &my_mpi_rank);
    int mpi_rank = 0;
    for (auto it = space_domain->macro_domain_root_set().begin();
         it != space_domain->macro_domain_root_set().end();
         ++it) {
      if (my_mpi_rank == mpi_rank) {
        macro_vertices.insert(*it);
      }
      ++mpi_rank;
    }

    LOG_(std::endl);
    space->define_subset_from_vertices(*space_domain, macro_vertices, false);
    LOG_("space->macro_boundary_.rows = " << space->macro_boundary_.rows() << std::endl);

    space->construct_mpiio_priorities(subdomain_comm);
    space->mpiio_priorities().log_dump(std::cout);
    
    daily = new TimeMesh((double)model_span.first(),
                         (double)model_span.second() + 1.0,
                         model_span.second() - model_span.first(),
                         epoch_calyear_,
                         TimeUnits::Day,
                         TimeBasis::Bspline2,
                         false);

    DateTimeHelper dth(epoch_calyear);
    int64_t climate_interval_length = 5;
    // Span should cover the model_span, in multiples of climate_interval_length,
    // starting at an even multiple
    DaySpan climate_span(
        (dth.get_year(model_span.first()) / climate_interval_length) * climate_interval_length,
        (dth.get_year(model_span.second()) / climate_interval_length) * climate_interval_length
        + climate_interval_length);
    int64_t num_intervals = (climate_span.second() - climate_span.first()) / climate_interval_length;
    climate = new TimeMesh((double)climate_span.first(),
                           (double)climate_span.second(),
                           num_intervals,
                           epoch_calyear_,
                           TimeUnits::Year,
                           TimeBasis::Bspline2,
                           false);

    int64_t order = 2;
    season = new TimeMesh(
        0.0, 1.0, order, epoch_calyear_,
        TimeUnits::Season, TimeBasis::Harmonic, true);

    init_mappers();
  }
  void init_mappers() {
    daily_map = new TimeMapper(daily, epoch_calyear_);
    climate_map = new TimeMapper(climate, epoch_calyear_);
    season_map = new TimeMapper(season, epoch_calyear_);
    space_point_map = new SpaceMapper(space, 0.0);
    space_cell_map = new SpaceMapper(space, 0.25);
    subdomain_map = new SpaceTimeMapper(NULL, space_point_map, space_cell_map, NULL, epoch_calyear_);
  }

  bool is_active() const {
    return (subdomain_comm != MPI_COMM_NULL);
  }

};

MeshBilevel* make_mesh_domain(const std::string& fileroot,
                              const std::vector<int64_t>& levels,
                              TimerHierarchy* timer)
{
  LOG_("Make mesh domain: start\n");
  MeshBilevel global(timer, false);
  LOG_(std::endl);
  global.set_levels(levels[0], levels[1]).import_python(fileroot);
  LOG_(std::endl);
  LOG_("global.macro_boundary_.rows = " << global.macro_boundary_.rows() << std::endl);
  //  PLOG_(WHEREAMI << std::endl);
  //  global.log_dump();
  
  // Locate the macro vertex closest to the centre point.
  std::pair<double, double> latlong_centre(50.839496642, 4.376331828);
  int64_t macro_centre_vertex = -1;
  // Locate nearest macro vertex
  {
    fmesh::Point point;
    Converter::latlong_to_euclidean(latlong_centre, point);
    double min_distance = std::numeric_limits<double>::infinity();
    for (size_t k = 0; k < global.macro_mesh().nV(); ++k) {
      fmesh::Point diff = point;
      diff.accum(global.macro_mesh().S(k), -1.0);
      double distance = diff.length();
      if (distance < min_distance) {
        macro_centre_vertex = k;
        min_distance = distance;
      }
    }

    std::pair<double, double> latlong_centre_found;
    Converter::euclidean_to_latlong(global.macro_mesh().S(macro_centre_vertex), latlong_centre_found);
    std::streamsize sz = std::cout.precision(10);
    PLOG_("Target: "
          << latlong_centre.first << "N "
          << latlong_centre.second << "E"
          << std::endl <<
          "Found:  "
          << latlong_centre_found.first << "N "
          << latlong_centre_found.second << "E"
          << std::endl);
    std::cout.precision(sz);
  }
  
  std::set<int64_t> macro_vertices;
  macro_vertices.insert(1);
  macro_vertices.insert(macro_centre_vertex);
  MeshBilevel* domain = new MeshBilevel(timer, false);
  domain->define_subset_from_vertices(global, macro_vertices, true);
  LOG_("domain->macro_boundary_.rows = " << domain->macro_boundary_.rows() << std::endl);
  
  //  PLOG_(WHEREAMI << std::endl);
  //  domain->log_dump();

  LOG_("Make mesh domain: done\n");
  return domain;
}


void do_map_obs(ModelObservationSource* obs,
                const std::vector<double>& comp_weights,
                SpaceTimeMapper* mapper,
                Eigen::SparseMatrix<double>& Amat,
                TimerHierarchy* timer) {
  unused_arg(timer);
  //  if (timer) timer->tic("Map observations for day " + to_string(obs->get_day()));
  //  LOG_("Do map obs: start\n");
  // Size: numobs x mapper->num_basis()*2 ( mean&range )

  assert(comp_weights.size() == 2);
  std::vector<double> weights = comp_weights;
  switch (obs->analysis_input().Observable()) {
    case AnalysisInput::observable_tmean: weights[0] *= 1.0; weights[1] *= 0.0; break;
    case AnalysisInput::observable_tmin: weights[0] *= 1.0; weights[1] *= -0.5; break;
    case AnalysisInput::observable_tmax: weights[0] *= 1.0; weights[1] *= +0.5; break;
    default: break;
  }
  //  PLOG_("weights = " << weights[0] << ", " << weights[1] << std::endl);

  bool is_cell = obs->is_cell();
      
  std::vector<Eigen::Triplet<double> > triplets(0);
  triplets.reserve(obs->size() * 2 * 3 * 3);

  for (int64_t obs_num = 0;
       obs_num < obs->size();
       ++obs_num) {
    SpaceTimeMapper::offset_type offset(obs_num, 0);
    int64_t count = mapper->mapping(
        SpaceTimeMapper::Point(obs->latlong()[obs_num].latitude,
                               obs->latlong()[obs_num].longitude,
                               obs->get_day(),
                               is_cell),
        true, // solar_time
        offset, // target row/col index
        !obs->is_mobile(), // Should the point be remembered or not?
        &weights,
        triplets);
    unused_arg(count);
  }

  //  PLOG_("Reserved " << obs->size() * 2 * 3 * 3 <<
  //        "\tRequired " << triplets.size() <<
  //        std::endl);

  Amat.resize(obs->size(), mapper->num_basis() * 2);
  Amat.setFromTriplets(triplets.begin(), triplets.end());
  Amat.makeCompressed();
  //  LOG_("Do map obs: done\n");
  //  if (timer) timer->toc();
}


int mpi_world_rank() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}
  

#define RNK "#" << mpi_world_rank() << ": "

void do_initialise_weights(SpaceTimeMapper* sub_mapper,
                           int64_t time_halo_depth,
                           int64_t space_halo_depth,
                           int64_t num_components,
                           VectorVariable& halo_hat,
                           VectorVariable& weights) {
  VectorVariable time_hat;
  VectorVariable space_hat;
  VectorVariable time_weights;
  VectorVariable space_weights;
  VectorVariable season_fun = VectorVariable::Ones(sub_mapper->num_basis_season());
  VectorVariable comp_fun = VectorVariable::Ones(num_components);
  if (sub_mapper->maps().time_->mesh() == NULL) {
    time_hat = VectorVariable::Ones(1);
    time_weights = VectorVariable::Ones(1);
  } else {
    sub_mapper->maps().time_->mesh()->get_weights(time_halo_depth, time_hat, time_weights);
  }
  LOG_(RNK << "Time halo hat:\n" << time_hat << std::endl);
  LOG_(RNK << "Time weights:\n" << time_weights << std::endl);
  
  if (sub_mapper->maps().space_point_->mesh() == NULL) {
    LOG_(RNK << std::endl);
    space_hat = VectorVariable::Ones(1);
    space_weights = VectorVariable::Ones(1);
  } else {
    LOG_(RNK << std::endl);
    // Relative halo depth:
    double space_halo = double(space_halo_depth) /
        double(sub_mapper->maps().space_point_->mesh()->macro_edge_size() + 1);
    sub_mapper->maps().space_point_->mesh()->initialise_weights(-1, space_halo);
    space_hat = sub_mapper->maps().space_point_->mesh()->macro_halo_hat();
    space_weights = sub_mapper->maps().space_point_->mesh()->macro_weights();
    sub_mapper->maps().space_point_->mesh()->clear_weights();
  }
    LOG_(RNK << std::endl);

  //  halo_hat.resize(sub_mapper->num_basis() * num_components, 1);
  weights = Eigen::kroneckerProduct(Eigen::kroneckerProduct(comp_fun, time_weights).eval(),
                                    Eigen::kroneckerProduct(space_weights, season_fun).eval());
    LOG_(RNK << std::endl);
  halo_hat = Eigen::kroneckerProduct(Eigen::kroneckerProduct(comp_fun, time_hat).eval(),
                                     Eigen::kroneckerProduct(space_hat, season_fun).eval());
    LOG_(RNK << std::endl);
}









NCHelper& accumulate_to_disk(NCHelper& nch,
                             const std::string& name,
                             MPIIOPrioritiesVarm& mpiio_prio,
                             const VectorVariable& values,
                             const VectorVariable& weights)
{
  VectorVariable result_buf(values.rows(), values.cols());
  auto prio_begin = mpiio_prio.begin();
  for (auto prio = prio_begin;
       prio != mpiio_prio.end();
       ++prio) {
    // Read data
    iget_varm(nch, name, *prio, result_buf.data()).err();
    nch.wait_all().stat();
    // Add values * weights
    auto result = prio.begin<double>(result_buf.data());
    auto weight = prio.cbegin<double>(weights.data());
    auto value = prio.cbegin<double>(values.data());
    auto value_end = prio.cend<double>(values.data());
    while (value != value_end) {
      *result += (*value) * (*weight);
      ++result;
      ++weight;
      ++value;
    }
    // Write result.
    iput_varm(nch, name, *prio, result_buf.data()).err();
    nch.wait_all().stat();
    MPI_Barrier(nch.mpi_comm());
  }
  return nch;
}








class ComponentDaily : public ComponentBase {
 protected:
  int64_t epoch_ = 1850;
  TimeMapper* time_map_;
  SpaceMapper* space_point_map_;
  SpaceMapper* space_cell_map_;
  SpaceTimeMapper mapper_;
  TimeMapper* next_block_map_;
  bool initialised_weights_;
 public:
  ComponentDaily(TimeMapper* time_map,
                 SpaceMapper* space_point_map,
                 SpaceMapper* space_cell_map,
                 TimerHierarchy* timer)
      : ComponentBase("Daily", NULL, timer),
        time_map_(time_map),
        space_point_map_(space_point_map),
        space_cell_map_(space_cell_map),
        mapper_(time_map_, space_point_map_, space_cell_map_, NULL, epoch_),
        next_block_map_(NULL),
        initialised_weights_(false)
  {
    update_next_block_map();
  }
  ~ComponentDaily() {
    if (next_block_map_) {
      TimeMesh* time_mesh = next_block_map_->mesh();
      delete next_block_map_;
      next_block_map_ = NULL;
      delete time_mesh;
    }
  }

  class Block : public BlockBase {
    TimeMapper* time_map_;
    SpaceTimeMapper* mapper_;
    ComponentDaily* mycomp_;
    MPIIOPrioritiesVarm mpiio_varm_;
   public:
    Block(ComponentDaily* comp, TimeMapper* time_map)
        : BlockBase(comp),
          time_map_(time_map),
          mycomp_(ComponentDaily::as_comp(comp))
    {
      LOG_("Construct new Daily block" << std::endl);
      set_obs_span(time_map_->mesh()->day_span_full());
      mapper_ = new SpaceTimeMapper(time_map_,
                                    comp->space_point_map_,
                                    comp->space_cell_map_,
                                    NULL,
                                    comp->epoch_);
      size_ = mapper_->num_basis() * 2;
      init_var();

      if (!mycomp_->initialised_weights_) {
        mycomp_->initialise_weights(mapper_);
      }

      make_mpiio_varm();
      make_Qprior();
    }

    virtual ~Block() {
      delete mapper_;
      if (time_map_) {
        TimeMesh* time_mesh = time_map_->mesh();
        delete time_map_;
        delete time_mesh;
      }
    }

    /** Convert spatial MPIIOPriorities to info for get/put_varm */
    void make_mpiio_varm() {
    LOG_(RNK << std::endl);
      mpiio_varm_.construct(2,
                            std::vector<bool>{true, true, true, false},
                            mycomp_->space_point_map_->mesh()->mpiio_priorities(),
                            *mapper_);
    LOG_(RNK << std::endl);
      mpiio_varm_.log_dump(std::cout);
    LOG_(RNK << std::endl);
    }
    const MPIIOPrioritiesVarm& mpiio_varm() const {
    LOG_(RNK << std::endl);
      return mpiio_varm_;
    }
    MPIIOPrioritiesVarm& mpiio_varm() {
    LOG_(RNK << std::endl);
      return mpiio_varm_;
    }

    void make_Qprior() {
      LOG("Construct Qprior for Daily block" << std::endl);
      Qprior_ = new SparseVariable(size_, size_);
    }

    void map_obs(ModelObservationSource* obs_source,
                 Eigen::SparseMatrix<double>& Amat) {
      //      mycomp_->timer_tic("Map observation source");
      std::vector<double> weights{1.0, 1.0}; // Mean and range attaches to observations.
      do_map_obs(obs_source, weights, mapper_, Amat, mycomp_->timer());
      //      mycomp_->timer_toc();
    }
  };

  static ComponentDaily* as_comp(ComponentBase* comp) {
    return dynamic_cast<ComponentDaily*>(comp);
  }
  static Block* as_block(BlockBase* block) {
    return dynamic_cast<Block*>(block);
  }
  
  BlockBase* new_block() {
    Block* block = new Block(this, next_block_map_);
    update_next_block_map();
    return block;
  }

  /** Update the next_block_map information
   *
   * To be called when the next block has just been constructed.
   * Does not deallocate the old one; that's done when the block is destructed.
   */
  void update_next_block_map() {
    int64_t num_intervals = 20;
    TimeMeshUnitless::SuperMeshInfo info(num_intervals, num_intervals / 2, 0);
    int64_t extra_shift = 0;
    if (next_block_map_ != NULL) {
      TimeMesh* next_mesh = next_block_map_->mesh();
      if (next_mesh->super_offset() + next_mesh->num_basis()
          >=
          time_map_->mesh()->num_basis()) {
        next_block_map_ = NULL;
        next_obs_span_ = DaySpan(false);
        return;
      }
      extra_shift = next_mesh->super_info().shift() + 1;
    }
    TimeMesh* time_mesh = new TimeMesh(time_map_->mesh(),
                                       info,
                                       extra_shift);
    next_block_map_ = new TimeMapper(time_mesh, epoch_);
    next_obs_span_ = time_mesh->day_span_full();
  }

  void initialise_weights(SpaceTimeMapper* sub_mapper) {
    int64_t time_halo = 2;
    int64_t space_halo = 2;
    do_initialise_weights(sub_mapper,
                          time_halo,
                          space_halo,
                          2,
                          halo_hat_,
                          weights_);
    initialised_weights_ = true;
  }

};





class ComponentSeason : public ComponentBase {
  int64_t epoch_ = 1850;
  SpaceTimeMapper* mapper_;
  bool initialised_weights_;
 public:
  ComponentSeason(SpaceMapper* space_point_map,
                  SpaceMapper* space_cell_map,
                  TimeMapper* season_map,
                  TimerHierarchy* timer)
      : ComponentBase("Season", NULL, timer),
        mapper_(NULL),
        initialised_weights_(false)
  {
    mapper_ = new SpaceTimeMapper(NULL, space_point_map, space_cell_map, season_map, epoch_);
    next_obs_span_ = DaySpan(true);
  }
  ~ComponentSeason() {
    delete mapper_;
  }


  class Block : public BlockBase {
    SpaceTimeMapper* mapper_;
    ComponentSeason* mycomp_;
    MPIIOPrioritiesVarm mpiio_varm_;
   public:
    Block(ComponentSeason* comp)
        : BlockBase(comp),
          mapper_(comp->mapper_),
          mycomp_(ComponentSeason::as_comp(comp))
    {
      LOG_("Construct new Season block" << std::endl);
      set_obs_span(comp->next_obs_span_);
      size_ = mapper_->num_basis() * 2;
      init_var();

      if (!mycomp_->initialised_weights_) {
        mycomp_->initialise_weights(mapper_);
      }

      make_mpiio_varm();
      make_Qprior();
      summary_log();
    }
    ~Block() {
    }

    /** Convert spatial MPIIOPriorities to info for get/put_varm */
    void make_mpiio_varm() {
      mpiio_varm_.construct(2,
                            std::vector<bool>{true, false, true, true},
                            mapper_->maps().space_point_->mesh()->mpiio_priorities(),
                            *mapper_);
      mpiio_varm_.log_dump(std::cout);
    }
    const MPIIOPrioritiesVarm& mpiio_varm() const {
      return mpiio_varm_;
    }
    MPIIOPrioritiesVarm& mpiio_varm() {
      return mpiio_varm_;
    }

    void make_Qprior() {
      LOG("Construct Qprior for Season block" << std::endl);
      Qprior_ = new SparseVariable(size_, size_);
    }

    void map_obs(ModelObservationSource* obs_source,
                 Eigen::SparseMatrix<double>& Amat) {
      std::vector<double> weights{1.0, 0.0}; // Only the Mean attaches to observations.
      do_map_obs(obs_source, weights, mapper_, Amat, mycomp_->timer());
    }
  };

  static ComponentSeason* as_comp(ComponentBase* comp) {
    return dynamic_cast<ComponentSeason*>(comp);
  }
  static Block* as_block(BlockBase* block) {
    return dynamic_cast<Block*>(block);
  }
  
  BlockBase* new_block() {
    Block* block = new Block(this);
    return block;
  }

  void initialise_weights(SpaceTimeMapper* sub_mapper) {
    int64_t space_halo = 2;
    do_initialise_weights(sub_mapper,
                          0,
                          space_halo,
                          2,
                          halo_hat_,
                          weights_);
    initialised_weights_ = true;
  }

};


std::string get_block_name(BlockBase* block) {
  if (ComponentDaily::as_block(block)) {
    return std::string("Daily");
  } else if (ComponentSeason::as_block(block)) {
    return std::string("Season");
  } else {
    return std::string("Unknown");
  }
}


class HandlerApplyQprior : public HandlerBase {
  std::string vec_;
  std::string result_;
 public:
  HandlerApplyQprior(const std::string& vec,
                     const std::string& result)
      : HandlerBase(),
        vec_(vec),
        result_(result)
  {
    vec_names(vec);
    vec_names(result);
  }
  void pre_obs() {
    PLOG_("\tApplyQprior to a " << get_block_name(block_) << " block.\n");
    block_->vec(result_) = block_->Qprior() * block_->vec(vec_);
    LOG_("TODO: apply cross-block precisions." << std::endl);
    block_->summary_log();
  }
};
class HandlerApplyAQeA : public HandlerBase {
  std::string vec_;
  std::string result_;
 public:
  HandlerApplyAQeA(const std::string& vec,
                   const std::string& result)
      : HandlerBase(),
        vec_(vec),
        result_(result)
  {
    vec_names(vec);
    vec_names(result);
  }
  typedef std::vector<std::string> downstream_type;

  void obs_single_block() {
    for (size_t k = 0; k < Amat().size(); ++k) {
      if (get_obs()->sources()[k]->size() == 0) {
        continue;
      }
      block_->vec(result_) +=
          Amat()[k].transpose() *
          get_obs()->sources()[k]->prec_uncorrelated().cwiseProduct(
              Amat()[k] * block_->vec(vec_));
    }
  }
  void obs_collect_from_downstream(downstream_type& downstream) {
    for (size_t k = 0; k < Amat().size(); ++k) {
      if (get_obs()->sources()[k]->size() == 0) {
        continue;
      }
      for (downstream_type::iterator it = downstream.begin();
           it != downstream.end();
           ++it) {
        ComponentBase* downstream_comp = block_->get_comp(*it);
        BlockBase* downstream_block = downstream_comp->find_covering_block(block_->obs_span());
        if (downstream_block) {
          block_->vec(result_) +=
              Amat()[k].transpose() *
              get_obs()->sources()[k]->prec_uncorrelated().cwiseProduct(
                  downstream_block->Amat()[k] * downstream_block->vec(vec_));
        } else {
          LOG_("Downstream block not found\n");
        }
      }
    }
  }    
  void obs_distribute_to_downstream(downstream_type& downstream) {
    VectorVariable QeAWv;
    for (size_t k = 0; k < Amat().size(); ++k) {
      if (get_obs()->sources()[k]->size() == 0) {
        continue;
      }
      QeAWv = get_obs()->sources()[k]->prec_uncorrelated().cwiseProduct(
          Amat()[k] *
          block_->get_comp_base()->halo_hat().cwiseProduct(block_->vec(vec_)));
      for (downstream_type::iterator it = downstream.begin();
           it != downstream.end();
           ++it) {
        ComponentBase* downstream_comp = block_->get_comp(*it);
        for (blocks_type::iterator bit = downstream_comp->blocks().begin();
             bit != downstream_comp->blocks().end();
             ++bit) {
          (*bit)->vec(result_) += (*bit)->Amat()[k].transpose() * QeAWv;
        }
      }
    }
  }
  void obs() {
    downstream_type downstream;
    PLOG_("\tApplyAQeA to a " << get_block_name(block_) << " block.\n");
    obs_single_block();
    if (ComponentDaily::as_block(block_)) {
      downstream.push_back("Season");
    }
    obs_collect_from_downstream(downstream);
    obs_distribute_to_downstream(downstream);

    block_->summary_log();
  }
};
class HandlerAccumulateQ : public HandlerBase {
 public:
  HandlerAccumulateQ() : HandlerBase() {
    sparse_names("Qpost");
    vec_names("AQey");
    sparse_names("AQeA");
  }
  void pre_obs() {
    PLOG_("\tAccumulateQ: Set Qpost to Qprior for a " << get_block_name(block_) << " block.\n");
    sparse("Qpost") = block_->Qprior();
    LOG_("TODO: Accumulate Qprior cross-blocks to Qpost.\n");
    block_->summary_log();
  }
  void obs() {
    PLOG_("\tAccumulateQ: Accumulate AQA from a day to a " << get_block_name(block_) << " block.\n");
    SparseVariable& AQeA = sparse("AQeA");
    for (size_t k = 0; k < Amat().size(); ++k) {
      if (get_obs()->sources()[k]->size() == 0) {
        continue;
      }
      AQeA += Amat()[k].transpose() * Amat()[k];
    }
    block_->summary_log();
  }
};
class HandlerApplyAQy : public HandlerBase {
 public:
  HandlerApplyAQy() : HandlerBase() {
    vec_names("AQey");
  }
  void obs() {
    PLOG_("\tApplyAQy: Apply AQy for a day to a " << get_block_name(block_) << " block.\n");
    
    for (size_t k = 0; k < Amat().size(); ++k) {
      if (get_obs()->sources()[k]->size() == 0) {
        continue;
      }
      vec("AQey") +=
          Amat()[k].transpose() *
          get_obs()->sources()[k]->prec_uncorrelated().cwiseProduct(
              get_obs()->sources()[k]->measurement());
    }
    block_->summary_log();
  }
};





ModelObservationSources* make_obs_sources(const std::string& path,
                                         SpaceTimeMapper* mask,
                                         int64_t epoch_calyear_obs,
                                         const DaySpan& obs_span,
                                         TimerHierarchy* timer) {
  // Make an input manager with specified base path for input data
  ModelObservationSources* obs_sources =
      new ModelObservationSources(path,
                                  mask,
                                  epoch_calyear_obs,
                                  timer);
  
  (*obs_sources)
      .add(AnalysisInput(
          AnalysisInput::source_satellite_ocean, 
          AnalysisInput::observable_tmean, 
          obs_span.first(), obs_span.second()))
      .add(AnalysisInput(
          AnalysisInput::source_satellite_land, 
          AnalysisInput::observable_tmin,
          obs_span.first(), obs_span.second()))
      .add(AnalysisInput(
          AnalysisInput::source_satellite_land, 
          AnalysisInput::observable_tmax,
          obs_span.first(), obs_span.second()))
      .add(AnalysisInput(
          AnalysisInput::source_satellite_ice, 
          AnalysisInput::observable_tmin,
          obs_span.first(), obs_span.second()))
      .add(AnalysisInput(
          AnalysisInput::source_satellite_ice, 
          AnalysisInput::observable_tmax,
          obs_span.first(), obs_span.second()))
      .add(AnalysisInput(
          AnalysisInput::source_insitu_land, 
          AnalysisInput::observable_tmin, 
          obs_span.first(), obs_span.second()))
      .add(AnalysisInput(
          AnalysisInput::source_insitu_land, 
          AnalysisInput::observable_tmax, 
          obs_span.first(), obs_span.second()))
      .add(AnalysisInput(
          AnalysisInput::source_insitu_ocean, 
          AnalysisInput::observable_tmean, 
          obs_span.first(), obs_span.second()));
  return obs_sources;
}


int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  MPI_Comm world_comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &world_comm);

  VersionRegistry::log();
  //  int64_t epoch_calyear_obs;
  int64_t start_day;
  int64_t end_day;

  std::string config_filename_default("config.cfg");
  LOG_(std::endl);
  Commandline cmdline(argc, argv, config_filename_default,
                      std::string("obs"), std::string("$Rev: 1299 $"));
  LOG_(std::endl);
  if (cmdline.handle_options()) {
    VERSIONS_PURGE();
    return cmdline.ret();
  }
  LOG_(std::endl);
    
  start_day = DateTimeHelper(cmdline.epoch_calyear_obs()).get_day(cmdline.start_calyear(), 1, 1) +
      cmdline.start_dayoffset();
  LOG_(std::endl);
  end_day = DateTimeHelper(cmdline.epoch_calyear_obs()).get_day(cmdline.end_calyear(), 1, 1) +
      cmdline.end_dayoffset();
  LOG_(std::endl);

  //    if (cmdline.vm().count("raw_datafile_path")) {
  std::cout << "Raw datafile path is: " << cmdline.raw_datafile_path() << std::endl;
  std::cout << "Python mesh filename prefix is: " << cmdline.python_mesh_prefix() << std::endl;
  //    }


  TimerHierarchy timer;
  timer.tic("Total");

  // Start and end day numbers since 01/01/1850
  DaySpan obs_span(start_day, end_day);
  //  obs_span = DaySpan(57374, 57401);

  DaySpan model_span(start_day, end_day);
  //  model_span = DaySpan(57374, 57401);

  timer.tic("Setup meshes and mesh mappers");
  std::vector<int64_t> levels{cmdline.macro_level(), cmdline.micro_level()};
  LOG_(std::endl);
  MeshBilevel* domain = make_mesh_domain(cmdline.python_mesh_prefix(),
                                         levels,
                                         NULL);
  
  LOG_(std::endl);
  ModelMeshes* meshes = new ModelMeshes(domain,
                                        model_span,
                                        cmdline.epoch_calyear_obs(),
                                        world_comm,
                                        NULL);
  if (!meshes->is_active()) {
    int myrank;
    MPI_Comm_rank(world_comm, &myrank);
    LOG_("Process with world rank " << myrank << " is not needed. Exiting." << std::endl);
    delete meshes;
    delete domain;
    VERSIONS_PURGE();
    LOG_(std::endl);
    MPI_Comm_free(&world_comm);
    LOG_(std::endl);
    MPI_Finalize();
    LOG_(std::endl);
    return 0;
  }
  LOG_(std::endl);
  
  LOG_(std::endl);
  timer.toc();

  domain->log_dump();
  meshes->space->log_dump();  

  timer.tic("Setup observation input manager");

  ModelObservationSources* obs_sources =
      make_obs_sources(cmdline.raw_datafile_path(),
                       meshes->subdomain_map,
                       cmdline.epoch_calyear_obs(),
                       obs_span,
                       &timer);
  timer.toc();

  timer.tic("Setup components");
  BlockScanner block_scanner(obs_sources, true);
  block_scanner
      .add(new ComponentDaily(meshes->daily_map,
                              meshes->space_point_map,
                              meshes->space_cell_map,
                              &timer))
      .add(new ComponentSeason(meshes->space_point_map,
                               meshes->space_cell_map,
                               meshes->season_map,
                               &timer));

  timer.toc().tic("Setup handlers");
  Handlers handlers(true);
  handlers
      .add(new HandlerApplyQprior("RHS", "result"))
      .add(new HandlerAccumulateQ())
      .add(new HandlerApplyAQeA("RHS2", "result2"))
      .add(new HandlerApplyAQy());

  timer.toc().tic("Run block scanner");
  block_scanner.scanner(handlers, obs_span, true);
  timer.toc();

  timer.toc(-1);
  std::cout << "\nTimings\n";
  timer.print();

  delete obs_sources;
  delete domain;
  delete meshes;
  VERSIONS_PURGE();
  
  MPI_Comm_free(&world_comm);
  MPI_Finalize();
  return 0;
}


