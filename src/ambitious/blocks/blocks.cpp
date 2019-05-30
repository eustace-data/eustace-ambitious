#include "ambitious/blocks/blocks.hpp"

#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(blocks_cpp, "$Revision: 1290 $")



using namespace EUSTACE;


ComponentBase::~ComponentBase() {
  while (!blocks_.empty()) {
    if (blocks_.front()) {
      delete blocks_.front();
    }
    blocks_.pop_front();
  }
}



DaySpan ComponentBase::obs_span() {
  DaySpan span(false);
  for (blocks_type::iterator bit = blocks_.begin();
       bit != blocks_.end();
       ++bit) {
    span.join((*bit)->obs_span());
  }
  return span;
}


BlockBase* ComponentBase::new_block() {
  return NULL;
}

void ComponentBase::new_blocks(Handlers& handlers,
                               int64_t current_day) {
  BlockBase* the_new_block;
  if (next_obs_span_.infinite()) { // Only a single block needed
    if (blocks_.empty()) {
      the_new_block = new_block();
      assert(the_new_block != NULL);
      blocks_.push_back(the_new_block);
      handlers.pre_obs(the_new_block);
    }
  } else {
    //  Add blocks while blocks_.back() is the last block the uses current_day
    while (!next_obs_span_.empty() &&
           (next_obs_span_.first() <= current_day)) {
      the_new_block = new_block();
      assert(the_new_block != NULL);
      blocks_.push_back(the_new_block);
      handlers.pre_obs(the_new_block);
    }
  }
}


BlockBase* ComponentBase::find_covering_block(const DaySpan& obs_span) {
  LOG_("Small span: " << obs_span << std::endl);
  for (blocks_type::reverse_iterator it = blocks_.rbegin();
       it != blocks_.rend();
       ++it) {
    LOG_("Block span: " << (*it)->obs_span() << std::endl);
    if ((*it)->obs_span().covers(obs_span)) {
      LOG_("Found span: " << (*it)->obs_span()<< std::endl);
      return (*it);
    }
  }
  LOG_("Span not found." << std::endl);
  return NULL;
}


void ComponentBase::map_obs() {
  timer_tic("Map locations");
  for (blocks_type::iterator bit = blocks_.begin();
       bit != blocks_.end();
       ++bit) {
    assert(*bit != NULL);
    (*bit)->map_obs();
  }
  timer_toc();
}

void ComponentBase::handle_obs(Handlers& handlers) {
  for (blocks_type::iterator bit = blocks_.begin();
       bit != blocks_.end();
       ++bit) {
    assert(*bit != NULL);
    handlers.obs(*bit);
  }
}

void ComponentBase::handle_post_obs(Handlers& handlers,
                                    int64_t current_day) {
  while (!blocks_.empty() &&
         ((!blocks_.front()->obs_span().infinite() &&
           (blocks_.front()->obs_span().empty() ||
            (blocks_.front()->obs_span().second() <= current_day))) ||
          (current_day == -1))) {
    assert(blocks_.front() != NULL);
    handlers.post_obs(blocks_.front());
    delete blocks_.front();
    blocks_.pop_front();
  }
}













VectorVariable& HandlerBase::vec(const std::string& name) {
  return block_->vec(name);
}
SparseVariable& HandlerBase::sparse(const std::string& name) {
  return block_->sparse(name);
}
Amats_type& HandlerBase::Amat() {
  return block_->Amat();
}
ModelObservationSources* HandlerBase::get_obs() {
  return block_->get_obs();
}




void Handlers::var_names(ComponentBase* comp) {
  for (handlers_type::iterator hit = handlers_.begin();
       hit != handlers_.end();
       ++hit) {
    for (var_names_type::iterator it = (*hit)->vec_names().begin();
         it != (*hit)->vec_names().end();
         ++it) {
      comp->add_vec(*it);
    }
    for (var_names_type::iterator it = (*hit)->sparse_names().begin();
         it != (*hit)->sparse_names().end();
         ++it) {
      comp->add_sparse(*it);
    }
  }
}


void Handlers::clear() {
  if (deallocate_) {
    for (Handlers::iterator it = handlers_.begin();
         it != handlers_.end();
         ++it) {
      delete(*it);
    }
  }
  handlers_.clear();
}


void Handlers::pre_obs(BlockBase* block) {
  PLOG_("Handlers::pre_obs, size = " << handlers_.size() << std::endl);
  for (Handlers::iterator it = handlers_.begin();
       it != handlers_.end();
       ++it) {
    (*it)->set_block(block).pre_obs();
  }
}
void Handlers::obs(BlockBase* block) {
  PLOG_("Handlers::obs, size = " << handlers_.size() << std::endl);
  for (Handlers::iterator it = handlers_.begin();
       it != handlers_.end();
       ++it) {
    (*it)->set_block(block).obs();
  }
}

void Handlers::post_obs(BlockBase* block) {
  PLOG_("Handlers::post_obs, size = " << handlers_.size() << std::endl);
  for (Handlers::iterator it = handlers_.begin();
       it != handlers_.end();
       ++it) {
    (*it)->set_block(block).post_obs();
  }
}
