/*
 * mfe_answer_v4.hh
 *
 *  Created on: June 26, 2012
 *      Author: jhuang
 */

#ifndef MFE_ANSWER_V4_HH_
#define MFE_ANSWER_V4_HH_

//#include <map>
#include <vector>
#include "RNAStack.hh"

typedef Basic_Subsequence<char, unsigned> TUSubsequence;

class mfeanswer_v4 {
 public:
  // Constructor
  mfeanswer_v4(): /*before_is(true), _multiloop(false),*/ empty_(false) {
    //?? hi_i_j_vector = (std::vector<Rope>*)malloc(sizeof(std::vector<Rope>*));
  }

  // Destructor
  ~mfeanswer_v4() {
    //?? free (hi_i_j_vector);
  }
  
  
  int energy;
//  Rope hi_rope;  // hishape is [36.5m,(,27,41.5,)] vs. hi_rope (,27,41.5,),36.5m
//##  Rope triplet;  // in ml() and hl() is the first item of 24,20%28 | in is() the complete of 24,20%28
//  Rope triplet_bi;
  // std::vector<Rope*> hi_i_j_rope;
  //?? std::vector<Rope>* hi_i_j_vector;
  //?? in hishapes.gap: merge(res.hi_i_j_vector, le.hi_i_j_vector, re.hi_i_j_vector); vs. res.hi_i_j_vector = le.hi_i_j_vector + re.hi_i_j_vector;
//##  std::vector<Rope> hi_i_j_vector;
  
  
//##  unsigned int /*i0,j0,*/i1,j1;  // i1,j1 for outmost basepair of a hairpin-unit  vs.  the positions of bulge- and hairpin loop is not completely correct
//  bool before_is;
  TUSubsequence firstStem;
  TUSubsequence secondStem;
  TUSubsequence beforeLastStem;
  TUSubsequence lastStem;
  TUSubsequence subword;
  int next_ss_i, next_ss_j, outermost_initstem_i, outermost_initstem_j;
  int before_last_ss_i, before_last_ss_j, after_next_ss_i, after_next_ss_j;
  TUSubsequence prevSubword;
  bool prevIsMl;
  int prevE;
  //String rep;
  String pp;
//  std::map<stem, int> stemEnergyMap;
  std::vector<stem> stems;
  std::vector<int> energies;
//##  bool _multiloop;
  bool empty_;

bool operator>(const mfeanswer_v4& other) const { return energy > other.energy; }
bool operator<(const mfeanswer_v4& other) const { return energy < other.energy; }
bool operator==(const mfeanswer_v4& other) const { return energy == other.energy; }
template <typename T> bool operator>(const T &other) const {return energy > other; }
template <typename T> bool operator<(const T &other) const {return energy < other; }
template <typename T> bool operator==(const T &other) const {return energy == other; }


mfeanswer_v4(int i) : energy(i), empty_(false) {}
mfeanswer_v4 operator+(const mfeanswer_v4 &other) const
{
assert(!empty_); assert(!other.empty_);
return mfeanswer_v4(energy + other.energy);
}
mfeanswer_v4 operator-(const mfeanswer_v4 &other) const
{
assert(!empty_);
if (other.empty_) return mfeanswer_v4(energy);
return mfeanswer_v4(energy - other.energy);
}
bool operator<=(const mfeanswer_v4& other) const {
assert(!empty_); assert(!other.empty_);
return energy <= other.energy;
}

};




inline std::ostream &operator<<(std::ostream &o, const mfeanswer_v4 &tuple) {
  o << /*'(' <<*/ tuple.energy
//  unsigned int i;
//  for (i=0; i< tuple.hi_i_j_vector.size(); i++) {
//    o << "|" << tuple.hi_i_j_vector[i]; 
//  }
   //o << ", hi_rope=" << tuple.hi_rope;
   //<< ", " << tuple.hi_i_j_rope
//    << ", (" << tuple.subword.i
//    << ", " << tuple.outermost_initstem_i
//    << "), (" << tuple.next_ss_i
//    << "," << tuple.next_ss_j
//    << "), (" << tuple.after_next_ss_i
//    << "," << tuple.after_next_ss_j
//    << "), (" << tuple.before_last_ss_i
//    << "," << tuple.before_last_ss_j
//    << "), (" << tuple.outermost_initstem_j
//    << "," << tuple.subword.j
//    << "), " 
   //<< ", " << tuple.before_is
   //<< ", " << tuple.firstStem
   //<< ", " << tuple.lastStem
   //<< ", " << tuple.subword
   //<< ", " << tuple.rep
   //<< ", " << tuple.pp
   /*<< ')'*/ ;
  return o;
}

inline void empty(mfeanswer_v4 &e) {e.empty_ = true; }
inline bool is_empty(const mfeanswer_v4 &e) { return e.empty_; }






/*
inline
std::vector<Rope> operator+(const std::vector<Rope> ropeVec1, const std::vector<Rope> ropeVec2)
{
 // Rope r;
 // append(r, a);
 // append(r, b);
 // return r;
      std::vector<Rope> ropeRes;
      std::vector<Rope>::iterator ropeIt;
      for ( ropeIt=ropeVec1.begin() ; ropeIt < ropeVec1.end(); ropeIt++ )
        ropeRes.push_back(*ropeIt);
      for ( ropeIt=ropeVec2.begin() ; ropeIt < ropeVec2.end(); ropeIt++ )
        ropeRes.push_back(*ropeIt);
      return ropeRes;
}

inline
void append(std::vector<Rope> &a, const Rope &b)
{
      a.push_back(b);
}
*/

/*
inline
void mergeVectors(const std::vector<String*>& a, const std::vector<String*>& b, std::vector<String*> * res)
{
 // Rope r;
 // append(r, a);
 // append(r, b);
 // return r;
    //  std::vector<String> res;
      std::vector<String*>::iterator stringIt;
      // this doesn't work, because no 
      for ( stringIt = a.begin(); stringIt < a.end(); stringIt++ )
        (*res).push_back(*stringIt);
      for ( stringIt=b.begin(); stringIt < b.end(); stringIt++ )
        (*res).push_back(*stringIt);
      // return res;
}
*/

// another option for this is using: void mergeRopeVectors(const std::vector<Rope*>& a, const std::vector<Rope*>& b, std::vector<Rope*> *res)

// with pointer
/*
inline
std::vector<Rope*> operator+(const std::vector<Rope*>& a, const std::vector<Rope*>& b)
{
  unsigned int i;
  std::vector<Rope *> res;
  for (i = 0; i < a.size(); ++i)
    res.push_back(a[i]);  
  for (i = 0; i < b.size(); ++i)
    res.push_back(b[i]);

  return res;  // copy the rector, it is not very efficient
}
*/

// without pointer
/*
 * initially, it should be void operator+(std::vector<Rope> *a, const std::vector<Rope> &b)
 * but because of the limitation of gap file,
 * it requires of (const mfeanswer_v4&, const mfeanswer_v4&)
 * due to the const have to create a new mfeanswer_v4
 */
inline
std::vector<Rope> operator+(const std::vector<Rope> &a, const std::vector<Rope> &b)
{
  unsigned int i;
  std::vector<Rope> res;
  for (i = 0; i < a.size(); ++i)
    res.push_back(a[i]);  
  for (i = 0; i < b.size(); ++i)
    res.push_back(b[i]);

  return res;  
}

  /* 
   * copy the vector, it is not very efficient 
   * ==> TODO: efficiently implementation (merge(std::vector<Rope> *res, const std::vector<Rope> &a, const std::vector<Rope> &b))
   * problem is in gap file, it is impossible to use &res, ask Georg if it is possible to use pointer in gap file
   *
   */
//?? inline
//?? void merge(std::vector<Rope> *res, const std::vector<Rope> *a, const std::vector<Rope> *b)
//?? {
//??   unsigned int i;
//??   for (i = 0; i < (*a).size(); ++i)
//??     (*res).push_back((*a)[i]);  
//??   for (i = 0; i < (*b).size(); ++i)
//??     (*res).push_back((*b)[i]);
//?? }

/* 
 * compile correctly, but segmentation default, since it needs dynamic deletion of such pointers
 */


// another option for this is using: void mergeRopeVectors(const std::vector<Rope*>& a, const std::vector<Rope*>& b, std::vector<Rope*> *res)
/*
inline
std::vector<Rope*> operator+(const std::vector<Rope>& a, const std::vector<Rope>& b)
{
  unsigned int i;
  std::vector<Rope *> res;
  for (i = 0; i < a.size(); ++i)
    res.push_back(&(a[i]));  
  for (i = 0; i < b.size(); ++i)
    res.push_back(&(b[i]));

  return res;  // copy the rector, it is not very efficient
}*/

// with pointer
/*
inline
void append(std::vector<Rope*> &a, Rope r)
{
  //unsigned int i;
  //std::vector<Rope *> res;
  //for (i = 0; i < a.size(); ++i)
  //  a.push_back(a[i]);  
  return a.push_back(&r);
}*/


// without pointer
/* new create mfeanswer_v4 is not constant, last level of mfeanswer_v4 is constant
 * following is fine
 * 		mfeanswer_v4 res = e;
 * 		    append(res.hi_i_j_vector,new_triplet);
 */
inline
void append(std::vector<Rope> &a, const Rope r)
{
  a.push_back(r);
}

inline
void clear(std::vector<Rope> &a)
{
  a.clear();
}

inline
void append(std::vector<Rope> &a, const std::vector<Rope> &b)
{
  unsigned int i;
  for (i = 0; i < b.size(); ++i)
    a.push_back(b[i]);
}



//?? inline
//?? void append(std::vector<Rope> *a, Rope r)
//?? { 
//??   (*a).push_back(r);
//?? }



/*
int compare(const std::vector<Rope>& left, const std::vector<Rope>& right) {
  auto leftIt = left.begin();
  auto rightIt = right.begin();
  auto diff = 0;
  while (leftIt != left.end() && rightIt != right.end()) {
    if (*leftIt != *rightIt) {
      diff++;
    }
    leftIt++;
    rightIt++;
  }

  // Account for different length vector instances
  if (0 == diff && (leftIt != left.end() || rightIt != right.end())) {
    diff = 1;
  }

  return diff;
}
*/


/*
inline
bool isEmpty(const std::vector<Rope> &a)
{
  if (a.size()==0)
    return true;
  else
    return false;
}
*/

//////////////////////////////////////////////////////
///////////// add stem operators for stem ////////////
// inline
// std::vector<Rope> operator+(const std::vector<Rope> &a, const std::vector<Rope> &b)
// {
//   unsigned int i;
//   std::vector<Rope> res;
//   for (i = 0; i < a.size(); ++i)
//     res.push_back(a[i]);  
//   for (i = 0; i < b.size(); ++i)
//     res.push_back(b[i]);
// 
//   return res;  
// }
// 
// inline
// void append(std::vector<Rope> &a, const Rope r)
// {
//   a.push_back(r);
// }
// 
// inline
// void clear(std::vector<Rope> &a)
// {
//   a.clear();
// }
// 
// inline
// void append(std::vector<Rope> &a, const std::vector<Rope> &b)
// {
//   unsigned int i;
//   for (i = 0; i < b.size(); ++i)
//     a.push_back(b[i]);
// }

#endif