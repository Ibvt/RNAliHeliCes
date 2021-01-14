/*
 * mfeanswer.hh
 *
 *  Created on: June 01, 2012
 *      Author: jhuang
 */

#ifndef MFE_ANSWER_HH_
#define MFE_ANSWER_HH_

#include "rtlib/adp.hh"

typedef Basic_Subsequence<char, unsigned> TUSubsequence;

struct mfeanswer {
  int energy;
  TUSubsequence firstStem;
  TUSubsequence lastStem;
  bool empty_;
  mfeanswer() : empty_(false) {}
bool operator>(const mfeanswer& other) const { return energy > other.energy; }
bool operator<(const mfeanswer& other) const { return energy < other.energy; }
bool operator==(const mfeanswer& other) const { return energy == other.energy; }
template <typename T> bool operator>(const T &other) const {return energy > other; }
template <typename T> bool operator<(const T &other) const {return energy < other; }
template <typename T> bool operator==(const T &other) const {return energy == other; }


mfeanswer(int i) : energy(i), empty_(false) {}
mfeanswer operator+(const mfeanswer &other) const
{
assert(!empty_); assert(!other.empty_);
return mfeanswer(energy + other.energy);
}
mfeanswer operator-(const mfeanswer &other) const
{
assert(!empty_);
if (other.empty_) return mfeanswer(energy);
return mfeanswer(energy - other.energy);
}
bool operator<=(const mfeanswer& other) const {
assert(!empty_); assert(!other.empty_);
return energy <= other.energy;
}

};

inline std::ostream &operator<<(std::ostream &o, const mfeanswer &tuple) {
  o << tuple.energy;
  return o;
}

inline void empty(mfeanswer &e) {e.empty_ = true; }
inline bool is_empty(const mfeanswer &e) { return e.empty_; }

#endif
