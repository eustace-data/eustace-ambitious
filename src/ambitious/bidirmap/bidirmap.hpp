#ifndef BIDIRMAP_HPP
#define BIDIRMAP_HPP

#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(bidirmap_hpp, "$Revision: 1272 $")

#include "ambitious/debuglog/debuglog.hpp"

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <map>
#include "ambitious/timer/timer.hpp"





/* 
 * C+98 used non-const iterators for map insertion hints.
 * C+11 uses const iterators for map insertion hints.
 */
#ifdef OLDCXX
#define BIDIRMAP_CONST
#else
#define BIDIRMAP_CONST const
#endif
 




template<typename KeyA, typename KeyB>
class BidirectionalMap;


#ifdef OLDCXX
template<class KeyA,class KeyB>
std::ostream& operator<<(std::ostream& output,
			 BidirectionalMap<KeyA,KeyB>& M) {
#else
template<class KeyA,class KeyB>
std::ostream& operator<<(std::ostream& output,
			 BIDIRMAP_CONST BidirectionalMap<KeyA,KeyB>& M) {
#endif
  for (typename BidirectionalMap<KeyA,KeyB>::IteratorA iterA=M.Abegin();
       iterA != M.Aend();
       ++iterA) {
    output << "A -> B = " << iterA->first << " -> " << iterA->second << std::endl;
  }
  for (typename BidirectionalMap<KeyA,KeyB>::IteratorB iterB=M.Bbegin();
       iterB != M.Bend();
       ++iterB) {
    output << "B -> A = " << iterB->first << " -> " << iterB->second << std::endl;
  }
       
  return output;
}



template<typename KeyA, typename KeyB>
class BidirectionalMap
{
public:
  typedef typename std::map<KeyA, KeyB > StorageAB;
  typedef typename std::map<KeyB, KeyA > StorageBA;
  typedef typename std::pair<KeyA, KeyB > KeyPairAB;
  typedef typename std::pair<KeyB, KeyA > KeyPairBA;
#ifdef OLDCXX
  typedef typename StorageAB::iterator IteratorA;
  typedef typename StorageBA::iterator IteratorB;
#else
  typedef typename StorageAB::const_iterator IteratorA;
  typedef typename StorageBA::const_iterator IteratorB;
#endif
  typedef typename std::pair<IteratorA, IteratorB > IteratorPair;
protected:
  StorageAB AtoB_;
  StorageBA BtoA_;

public:
  BidirectionalMap() : AtoB_(), BtoA_() {}

  friend
  std::ostream& operator<< <> (std::ostream& output,
			       const BidirectionalMap<KeyA,KeyB>& M);

  typename StorageAB::size_type size() const {
    return AtoB_.size();
  }

  void clear() {
    AtoB_.clear();
    BtoA_.clear();
  }
  /** Insert a key pair and return a pair of iterators to the next
   *  elements, so they are useful as hints for a following call to
   *  a assign with hint, in an ordered assignment loop.
   */
  IteratorPair insert(const KeyA& A, const KeyB& B) {
    IteratorPair position = ABend();
    return insert(position, A, B);
  }
  /** Insert a key pair and return a pair of iterators to the next
   *  elements, so they are useful as hints for a following call to
   *  assign with hint, in an ordered assignment loop.
   */
  IteratorPair insert(IteratorPair& position, const KeyA& A, const KeyB& B) {
    IteratorA iterA = AtoB_.insert(position.first, KeyPairAB(A, B));
    IteratorB iterB = BtoA_.insert(position.second, KeyPairBA(B, A));
#ifdef OLDCXX
    return IteratorPair(iterA, iterB);
#else
    return IteratorPair(++iterA, ++iterB);
#endif
  }
  StorageAB& AtoB() {
    return AtoB_;
  }
  StorageBA& BtoA() {
    return BtoA_;
  }
  IteratorA AtoB_iter(const KeyA& A) {
    return AtoB_.find(A);
  }
  IteratorB BtoA_iter(const KeyB& B) {
    return BtoA_.find(B);
  }
  const KeyB& AtoB(const KeyA& A) {
    IteratorA iter = AtoB_.find(A);
    assert(iter != AtoB_.end());
    return iter->second;
  }
  const KeyA& BtoA(const KeyB& B) {
    IteratorB iter = BtoA_.find(B);
    assert(iter != BtoA_.end());
    return iter->second;
  }
  IteratorPair ABbegin() BIDIRMAP_CONST {
    return IteratorPair(AtoB_.begin(), BtoA_.begin());
  }
  IteratorPair ABend() BIDIRMAP_CONST {
    return IteratorPair(AtoB_.end(), BtoA_.end());
  }
  IteratorA Abegin() BIDIRMAP_CONST {
    return AtoB_.begin();
  }
  IteratorB Bbegin() BIDIRMAP_CONST {
    return BtoA_.begin();
  }
  IteratorA Aend() BIDIRMAP_CONST {
    return AtoB_.end();
  }
  IteratorB Bend() BIDIRMAP_CONST {
    return BtoA_.end();
  }
};


#endif
