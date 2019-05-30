/** Some file documentation
 * @file
 */
#ifndef RARANDOM_HPP
#define RARANDOM_HPP

#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(rarandom_hpp, "$Revision: 1292 $")


/* Test running times for random number generators
 */
#include <iostream>
#include <string>
#include <vector>
#include <boost/date_time.hpp>
#include <boost/random.hpp>
#include <random>
#include <Eigen/Dense>
#include <cmath>

// Read about old random number generator in std/tr1/boost here:
//   https://gcc.gnu.org/bugzilla/show_bug.cgi?id=44653

/** Random access random number generator
 *
 * Subsequences \f$S_i = [i N_s, (i+1)N_s-1]\f$
 * 
 * Reverse subsequence lookup \f$S(k)=i\f$ such that \f$k\in S_i\f$
 * 
 * Time to generate one number \f$T_g\f$
 *
 * Time to reset the generator \f$T_r\f$
 *
 * Expected time to find \f$k\f$ given an initial random state
 * \f{align*}{
 *   E[T] &=
 *   P[k\in S(k_0),\, k > k_0] (k-k_0) T_g + P[k\in S(k_0),\, k < k_0] (T_r+(k-S(k_0)N_s)T_g)
 *   \\&\phantom{= }
 *   + P[k\not\in S(k_0)] (T_r + (k-S(k)N_s)T_g)
 * \f}
 *
 * Needs to be tied to a specific distribution to allow random access seeks
 */
template <typename Scalar, typename Generator, typename Distribution>
class random_access_variates {
  typedef std::vector<Generator> GeneratorVector;
  struct FileHeader {
    uint64_t magic;		/**< magic file type ID */
    uint64_t sublength;		/**< the spacing between generators */
    uint64_t size;              /**< the number of generators */
    uint64_t generator_size;	/**< size of each generator entry */
  };

  Generator generator_;		/**<  */
  Distribution distribution_;	/**<  */
  uint64_t  next_; // Index of the next variate
  uint64_t sublength_; // Distance between saved states
  GeneratorVector states_; // Saved states at sublength_ * 0, 1, 2, 3, etc

public:
  /** Magic number for file type ID.
   *  sum(as.integer(charToRaw("EUSTACE")) * 256^(0:6)) calculated in full precision
   */
  static const uint64_t magic_filetype_ID_ = 19495721259717957L;
  static const uint64_t sizeof_FileHeader = sizeof(FileHeader);
  static const uint64_t sizeof_Generator = sizeof(Generator);
  
  random_access_variates(Distribution distribution,
			 uint64_t sublength) : generator_(),
					       distribution_(distribution),
					       next_(0),
					       sublength_(sublength),
					       states_(0) {
    generator_.seed((uint32_t)magic_filetype_ID_);
    states_.push_back(generator_);
  }
  random_access_variates(Generator generator,
			 Distribution distribution,
			 uint64_t sublength) : generator_(generator),
					       distribution_(distribution),
					       next_(0),
					       sublength_(sublength),
					       states_(0) {
    states_.push_back(generator_);
  }
  random_access_variates(std::string filename,
			 Distribution distribution) : generator_(),
						      distribution_(distribution),
						      next_(0),
						      sublength_(1),
						      states_(0) {
    load(filename);
  }

  void save(std::string filename) {
    std::ofstream file(filename, std::ofstream::binary | std::ifstream::out | std::ofstream::trunc);
    FileHeader header;
    header.magic = magic_filetype_ID_;
    header.sublength = sublength_;
    header.size = states_.size();
    header.generator_size = sizeof(Generator);
    file.write((char*)&header, sizeof(header));
    for (typename GeneratorVector::iterator i = states_.begin(); i != states_.end(); ++i) {
      file.write((char*)&(*i), sizeof(*i));
    }
  }

  void load(std::string filename) {
    std::ifstream file(filename, std::ifstream::binary | std::ifstream::in);
    FileHeader header;
    file.read((char*)&header, sizeof(header));
    if (header.magic != magic_filetype_ID_) {
      std::cerr << "Expected magic number " << magic_filetype_ID_
		<< ", found " << header.magic << std::endl;
    }
    if (header.generator_size != sizeof(Generator)) {
      std::cerr << "Expected generator size " << sizeof(Generator)
		<< ", found " << header.generator_size << std::endl;
    }
    sublength_ = header.sublength;
    std::cout << "sizeof(Generator) = " << sizeof(Generator) << std::endl;
    states_.resize(header.size);
    for (uint64_t i=0; i < header.size; ++i) {
      file.read((char*)&states_[i], sizeof(states_[i]));
    }
    generator_ = states_[0];
    next_ = 0;
  }

  /** Index of the next random number 
   * 
   * 
   * 
   * @return 
   */
  uint64_t next_index() const { return next_; }
  /** Number of stored waypoints 
   * 
   * 
   * 
   * @return 
   */
  uint64_t size() const { return states_.size(); };
  /** Length of the subsegments 
   * 
   * @return 
   */
  uint64_t sublength() const { return sublength_; };

  void seek(uint64_t index) {
    if (next_ == index) {
      return;
    }
    // Are we in the right subsequence, before the sought index?
    uint64_t subsequence_next = next_ / sublength_;
    uint64_t subsequence_index = index / sublength_;
    if (subsequence_index >= states_.size()) {
      if (subsequence_next+1 != states_.size()) {
	// Go to last available state
	generator_ = states_[states_.size()-1];
	next_ = (states_.size() - 1) * sublength_;
      }
    } else if ((subsequence_next != subsequence_index) | (next_ > index)) {
      // Go to the required subsequence state
	generator_ = states_[subsequence_index];
	next_ = subsequence_index * sublength_;
    }
    // Generate numbers until ready
    while (next_ < index) {
      (*this)();
    }
  }

  Scalar operator()() {
    if ((next_ % sublength_ == 0) &
	(next_ >= states_.size() * sublength_)) {
      states_.push_back(generator_);
    }
    ++next_;
    distribution_.reset();
    return distribution_(generator_);
  }

  Scalar operator()(uint64_t index) {
    seek(index);
    return (*this)();
  }
};


#endif
