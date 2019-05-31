#ifndef BLOCKS_HPP
#define BLOCKS_HPP

#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(blocks_hpp, "$Revision: 1290 $")

#include "ambitious/debuglog/debuglog.hpp"
#include <eustace/analysis/inputmanager.h>
#include <iostream>
#include <set>
#include <map>
#include <deque>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "ambitious/meshes/meshes.hpp"
#include "ambitious/timer/timer.hpp"
#include "ambitious/to_string/to_string.hpp"
#include "ambitious/observe/observe.hpp"

using namespace EUSTACE;

template <typename T> void unused_arg(const T&) {} // Helper function.




typedef std::set<std::string> var_names_type;

class ComponentBase;
class BlockBase;
class HandlerBase;
typedef std::map<std::string, ComponentBase*> components_type;
typedef std::deque<BlockBase*> blocks_type;

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorVariable;
typedef Eigen::SparseMatrix<double> SparseVariable;
typedef std::vector<Eigen::SparseMatrix<double> > Amats_type;

class HandlerBase {
 protected:
  var_names_type vec_names_;
  var_names_type sparse_names_;
  BlockBase* block_;
  components_type* comps_;
 public:
  HandlerBase()
      : vec_names_(),
        sparse_names_(),
        block_(NULL),
        comps_(NULL)
  {
  }
  virtual ~HandlerBase() {
  }
  HandlerBase& set_block(BlockBase* block) {
    block_ = block;
    return *this;
  }
  HandlerBase& set_comps(components_type* comps) { 
    comps_ = comps;
    return *this;
  }
  HandlerBase& set(BlockBase* block,
                   components_type* comps) {
    block_ = block;
    comps_ = comps;
    return *this;
  }
  HandlerBase& vec_names(const std::string& name) {
    vec_names_.insert(name);
    return *this;
  }
  HandlerBase& sparse_names(const std::string& name) {
    sparse_names_.insert(name);
    return *this;
  }
  var_names_type& vec_names() {
    return vec_names_;
  }
  var_names_type& sparse_names() {
    return sparse_names_;
  }

  VectorVariable& vec(const std::string& name);
  SparseVariable& sparse(const std::string& name);
  Amats_type& Amat();
  ModelObservationSources* get_obs();
  
  virtual void pre_obs() {
  }
  virtual void obs() {
  }
  virtual void post_obs() {
  }
};

class Handlers {
 public:
  typedef std::vector<HandlerBase*> handlers_type;
  typedef handlers_type::iterator iterator;
 protected:
  handlers_type handlers_;
  const bool deallocate_;
  components_type* comps_;
 public:
  Handlers(bool deallocate)
      : deallocate_(deallocate),
        comps_(NULL)
  {
  }
  ~Handlers() {
    clear();
  }

  /** Called by BlockScanner::scanner to add names to components
   *  before performing any operations
   */
  void var_names(ComponentBase* comp);

  Handlers& add(HandlerBase* handler) {
    handlers_.push_back(handler);
    return *this;
  }
  void clear();
  handlers_type::iterator begin() {
    return handlers_.begin();
  }
  handlers_type::iterator end() {
    return handlers_.end();
  }
  Handlers& set_comps(components_type* comps) { 
    comps_ = comps;
    for (iterator it = handlers_.begin();
         it != handlers_.end();
         ++it) {
      (*it)->set_comps(comps_);
    }
    return *this;
  }
  components_type* get_comps() {
    return comps_;
  }
  ComponentBase* get_comp(const std::string& name) {
    return (*comps_)[name];
  }
  void pre_obs(BlockBase* block);
  void obs(BlockBase* block);
  void post_obs(BlockBase* block);
};

/** Hold information about a model component
 *
 *
 *
 */
class ComponentBase {
 protected:
  components_type* comps_;
  std::string name_;
  blocks_type blocks_;
  var_names_type vec_names_;
  var_names_type sparse_names_;
  DaySpan next_obs_span_;
  ModelObservationSources* obs_;
  VectorVariable halo_hat_; /// Deriving class is responsible for filling
  VectorVariable weights_; /// Deriving class is responsible for filling
  TimerHierarchy* timer_;
  TIMER_HELPERS

 public:
  
  ComponentBase(const std::string& name,
                ModelObservationSources* obs = NULL,
                TimerHierarchy* the_timer = NULL)
      : comps_(NULL),
        name_(name),
        blocks_(0),
        vec_names_(),
        sparse_names_(),
        next_obs_span_(true), // Infinite span by default
        obs_(obs),
        halo_hat_(),
        weights_(),
        timer_(the_timer)
  {
  }
  virtual ~ComponentBase();

  const std::string& name() const {
    return name_;
  }
  blocks_type& blocks() {
    return blocks_;
  }
  TimerHierarchy* timer() {
    return timer_;
  }
  void add_vec(const std::string& name) {
    vec_names_.insert(name);
  }
  void add_sparse(const std::string& name) {
    sparse_names_.insert(name);
  }
  var_names_type& vec_names() {
    return vec_names_;
  }
  var_names_type& sparse_names() {
    return sparse_names_;
  }

  DaySpan obs_span();

  ComponentBase& set_obs(ModelObservationSources* obs) {
    obs_ = obs;
    return *this;
  }
  ModelObservationSources* get_obs() {
    return obs_;
  }
  ComponentBase& set_comps(components_type* comps) {
    comps_ = comps;
    return *this;
  }
  components_type* get_comps() {
    return comps_;
  }
  ComponentBase* get_comp(const std::string& name) {
    return (*comps_)[name];
  }


  /** Initialise new blocks that start on the current day
   *
   * If this is the first time blocks are setup, all blocks that
   * contain the specified day are initialised.
   *
   * After adding a new block, the Handlers::handle_pre_obs method is called.
   */
  void new_blocks(Handlers& handlers, int64_t current_day);
  void map_obs();
  void handle_obs(Handlers& handlers);
  void handle_post_obs(Handlers& handlers, int64_t current_day);

  virtual BlockBase* new_block();

  BlockBase* find_covering_block(const DaySpan& obs_span);

  VectorVariable& halo_hat() {
    return halo_hat_;
  }
  VectorVariable& weights() {
    return weights_;
  }
    
};


class BlockBase {
 public:
  typedef std::map<std::string, VectorVariable> VectorVariables;
  typedef std::map<std::string, SparseVariable> SparseVariables;
  typedef std::vector<Eigen::Triplet<double> > mapping_triplet_type;

 protected:
  ComponentBase* comp_; // What component does the block belong to?
  size_t size_;
  DaySpan obs_span_;
  VectorVariables vec_var_;
  SparseVariables sparse_var_;
  //  std::vector<mapping_triplet_type> mapping_triplets_;
  Amats_type Amat_;
  SparseVariable* Qprior_;

  void set_obs_span(const DaySpan& obs_span) {
    obs_span_ = obs_span;
  }
  
 public:
  BlockBase(ComponentBase* comp)
      : comp_(comp),
        size_(0),
        obs_span_(true),
        vec_var_(),
        sparse_var_(),
        Qprior_(NULL)
  {
  }
  
  virtual ~BlockBase() {
      if (Qprior_) {
        delete Qprior_;
      }
  }

  VectorVariable& vec(const std::string& name) {
    return vec_var_[name];
  };

  SparseVariable& sparse(const std::string& name) {
    return sparse_var_[name];
  };
  Amats_type& Amat() {
    return Amat_;
  }
  SparseVariable& Qprior() {
    assert(Qprior_ != NULL);
    return *Qprior_;
  }

  ComponentBase* component() {
    return comp_;
  }
  const std::string& component_name() {
    return comp_->name();
  }
  const DaySpan& obs_span() {
    return obs_span_;
  }
  ModelObservationSources* get_obs() {
    return comp_->get_obs();
  }
  
  /** Must be called by the constructor of a deriving class, after setting size_. */
  void init_var() {
    for (var_names_type::iterator it = comp_->vec_names().begin();
         it != comp_->vec_names().end();
         ++it) {
      //      input_.emplace(*it, VectorVariable(size_));
      vec_var_.insert(VectorVariables::value_type(*it, VectorVariable(size_)));
    }
    for (var_names_type::iterator it = comp_->sparse_names().begin();
         it != comp_->sparse_names().end();
         ++it) {
      sparse_var_.insert(SparseVariables::value_type(*it, SparseVariable(size_, size_)));
    }
  }

  /** Summarise block contents. */
  void summary_log() {
    PLOG_("Block summary:" << std::endl);
    PLOG_("\tSpan: " << obs_span() << "\n");
    for (Amats_type::iterator it = Amat_.begin();
         it != Amat_.end();
         ++it) {
      PLOG_("\t" << "Amat: "
            << (*it).rows()
            << " x "
            << (*it).cols()
            << ", nnz = "
            << (*it).nonZeros()
            << std::endl);
    }
    for (SparseVariables::iterator it = sparse_var_.begin();
         it != sparse_var_.end();
         ++it) {
      PLOG_("\t" << (*it).first << ": "
            << (*it).second.rows()
            << " x "
            << (*it).second.cols()
            << ", nnz = "
            << (*it).second.nonZeros()
            << std::endl);
    }
    for (VectorVariables::iterator it = vec_var_.begin();
         it != vec_var_.end();
         ++it) {
      PLOG_("\t" << (*it).first << ": "
            << (*it).second.rows()
            << " x "
            << (*it).second.cols()
            << ", nnz = "
            << (*it).second.nonZeros()
            << std::endl);
    }
  }

  virtual void map_obs(ModelObservationSource* obs_source,
                       Eigen::SparseMatrix<double>& Amat) {
    PLOG_("virtual map_obs not overridden" << std::endl);
    unused_arg(obs_source);
    unused_arg(Amat);
  }
  void map_obs() {
    ModelObservationSources::sources_type& src = get_obs()->sources();
    PLOG_("obs_->sources().size() = " << src.size() << "\n");
    Amat_.resize(src.size());
    std::vector<Eigen::SparseMatrix<double> >::iterator Aiter = Amat_.begin();
    for (ModelObservationSources::sources_type::iterator iter = src.begin();
         iter != src.end();
         ++iter) {
      map_obs(*iter, (*Aiter));
      ++Aiter;
    }
  }

  components_type* get_comps() {
    return comp_->get_comps();
  }
  ComponentBase* get_comp_base() {
    return comp_;
  }
  ComponentBase* get_comp(const std::string& name) {
    return (*get_comps())[name];
  }
};





/** Class to hold and iterate over model subblocks */
class BlockScanner {
 protected:
  components_type components_;
  ModelObservationSources* obs_sources_;
  bool deallocate_components_;

 public:
  BlockScanner(ModelObservationSources* obs_sources,
               bool deallocate_components)
      : components_(),
        obs_sources_(obs_sources),
        deallocate_components_(deallocate_components)
  {
  }
  ~BlockScanner() {
    clear();
  }
  void clear() {
    if (deallocate_components_) {
      for (typename components_type::iterator it = components_.begin();
           it != components_.end();
           ++it) {
        delete (*it).second;
      }
    }
    components_.clear();
  }

  BlockScanner& add(ComponentBase* component) {
    if (component) {
      components_.insert(
          components_type::value_type(component->name(),
                                      component));
      component->set_comps(&components_);
    }
    return *this;
  }

  /** Setup variable names from handler information
   */
  void set_var_names(Handlers& handlers) {
    for (typename components_type::iterator it = components_.begin();
         it != components_.end();
         ++it) {
      handlers.var_names((*it).second);
    }
  }


  /** Setup handlers
   */
  void init_handlers(Handlers& handlers) {
    set_var_names(handlers);
    handlers.set_comps(&components_);
  }

  /** Initialise new blocks that start on the current day
   *
   * If this is the first time blocks are setup, all blocks
   * that start on or before the current day initialised.
   */
  void new_blocks(Handlers& handlers, int64_t current_day) {
    for (typename components_type::iterator it = components_.begin();
         it != components_.end();
         ++it) {
      (*it).second->new_blocks(handlers, current_day);
    }
  }
  /** Construct observation A-matrices for every block
   */
  void map_obs() {
    for (typename components_type::iterator it = components_.begin();
         it != components_.end();
         ++it) {
      (*it).second->set_obs(obs_sources_).map_obs();
    }
  }
  /** Handle observation info for all blocks
   *
   * Perform observation driven operations and accumulate A matrix info
   */
  void handle_obs(Handlers& handlers) {
    for (typename components_type::iterator it = components_.begin();
         it != components_.end();
         ++it) {
      (*it).second->handle_obs(handlers);
    }
  }
  /** Handle blocks that end on or before the current day
   *
   * If current_day is -1, all blocks are handled.
   *
   * Each handled block is removed from its queue
   */
  void handle_post_obs(Handlers& handlers, int64_t current_day) {
    for (typename components_type::iterator it = components_.begin();
         it != components_.end();
         ++it) {
      (*it).second->handle_post_obs(handlers, current_day);
    }
  }

  DaySpan obs_span() {
    DaySpan span(false);
    for (typename components_type::iterator it = components_.begin();
         it != components_.end();
         ++it) {
      span.join((*it).second->obs_span());
    }
    return span;
  }

  /** Scan through data and perform operations
   *
   * @param handlers The operations to be performed
   * @param days The observation time span
   * @param finalise If true, compute post observation
   *   operations on all remaining blocks after the last day.
   *   If set to false, the scanner can continue in a
   *   subsequent call.
   */
  void scanner(Handlers& handlers,
               const DaySpan& days,
               bool finalise) {
    init_handlers(handlers); // Setup variable names, etc
    assert(!days.empty() && !days.infinite());
    for (int64_t current_day = days.first();
         current_day <= days.second();
         ++current_day) {
      PLOG_("Scanning day " << current_day << std::endl);
      new_blocks(handlers, current_day); // Initialise blocks and perform pre-data operations
      obs_sources_->retrieve_day(current_day); // Get observations
      map_obs(); // Construct A matrices
      handle_obs(handlers); // Perform observation based accumulation
      handle_post_obs(handlers, current_day); // Perform post-data operations on blocks
    }
    if (finalise) {
      handle_post_obs(handlers, -1);
    }
  }
};


#endif
