/* {{{

    This file is part of gapc (GAPC - Grammars, Algebras, Products - Compiler;
      a system to compile algebraic dynamic programming programs)

    Copyright (C) 2008-2011  Georg Sauthoff
         email: gsauthof@techfak.uni-bielefeld.de or gsauthof@sdf.lonestar.org

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

}}} */

#ifndef HASH_HH
#define HASH_HH

#include "bitops.hh"
#include "move.hh"

#include "ref.hh"

#include "pool.hh"

#include "vector_sparse.hh"

#include "hash_stats.hh"

#include <ostream>

#include <boost/cstdint.hpp>
#include <cassert>
#include <vector>

//#include "hishape_answer_v5.hh"

using std::swap;

// XXX remove
#include <iostream>


#ifndef HASH_INITIAL
#define HASH_INITIAL 16
#endif

#ifndef HASH_LOAD_FACTOR
#define HASH_LOAD_FACTOR 75
#endif

#ifndef HASH_SHRINK
#define HASH_SHRINK EnableShrink
#endif

#ifndef HASH_STAT
#define HASH_STAT NoStats
#endif


#include <typeinfo>

inline uint32_t hashable_value(int t) { return t; }

namespace Hash {


  // FIXME
  inline uint32_t hashable_value(size_t t) { return t; }


  struct Multhash {
    uint32_t operator()(uint32_t kk, uint32_t m) const
    {
      uint32_t i = 32 - count_leading_zeroes(m) - 1;
      assert(m == size_to_next_power(m));
      uint32_t s = 2654435769;  
      uint64_t k = kk;
      uint64_t t = s * k & 0x00000000FFFFFFFF;
      return t >> (32 - i);
    } 
  };

  template <typename T>
  struct Size2pow {
    T initial() const
    {
      //return 512;
      //return 4;
      //return 8;
      //return 32;
      //return 16;
      return HASH_INITIAL;
    }
    T expand(T size) const
    {
      return size * 2;
    }

    T filter(T size) const
    {
      return size_to_next_power(size);
    }
  };

  struct EnableShrink {
    enum { value = true };
  };
  struct DisableShrink {
    enum { value = false };
  };

  
  //////////////////////////////////////////////////////////////////
  //////////////////// class Default_Inspector /////////////////////
  // Inspector is a class 
  template <typename T, typename U = uint32_t>
  struct Default_Inspector {
    T init(const T &x) const { return x; }
    U hash(const T &x) const
    {
      return hashable_value(x);
    }
    void update(T &dst, const T &src) const
    {
    }
    bool equal(const T &a, const T &b) const
    {
      return a == b;
    }
    bool filter() const { return false; }
    bool filter(const T &x) const { assert(0); return false; }
    void finalize(T &src) const
    {
       //std::cerr << "dd" << std::endl;
    }

    uint32_t k() const { assert(0); return 0; }
    bool cutoff() const { return false; }
    bool equal_score(const T &a, const T &b) const
    {
      assert(0);;
      return false;
    }
    struct compare {
      bool operator()(const T &a, const T &b) const
      {
	assert(0);;
	return false;
      }
    };
  };

  
  ////////////////////////////////////////////////////
  //////////////////// class Set /////////////////////
  
  //## MOST IMPORTANT METHODS in this class ##
  /* 
   * 1.0 push_back(answers, ans); ==> 1.1 add(const T &t) ==> 1.2 void add(const T &t, bool update) ==> 1.3 bool insert(U index, const T &t, bool update)
   * 2.0 hash_filter( answers);   ==> 2.1 fileter() ==> filter sth. for example kbest, unique(...) does not work with kbest due to following code
   * 3.0 finalize( answers);      ==> 3.1 finalize() ==> inspector.finalize(*i); ==> DO NOTHING
   *
   */
  template <typename T,
            class Inspector = Default_Inspector<T>,
            typename U = uint32_t,
            typename Hash_Policy = Multhash,
            template<typename> class Resize_Policy = Size2pow,
            typename Shrink_Policy = HASH_SHRINK,
            typename Stat_Policy = HASH_STAT,
            unsigned int load_factor = HASH_LOAD_FACTOR
            >
  class Set {


    private:
      Vector_Sparse<T, U> array;  // array, U should be a hashcode
      std::vector<bool> init;
      U used_;
#ifndef NDEBUG
      bool finalized;
#endif

      Inspector inspector;
      Resize_Policy<U> resize_policy;
      Hash_Policy hash_policy;

      void rehash(U i)
      {
        assert(i > array.size());
        U size = resize_policy.filter(i);
        Vector_Sparse<T, U> a(size);
        swap(a, array);
        std::vector<bool> b(size);
        swap(init, b);
        used_ = 0;
        for (typename Vector_Sparse<T, U>::iterator i = a.begin();
             i != a.end(); ++i)
          add(*i, false);
      }
      bool loaded() const
      {
        unsigned int load = double(used_)/double(array.size()) * 100.0;
        return load >= load_factor;
      }
      void expand()
      {
        U i = resize_policy.expand(array.size()+1);
        rehash(i);
      }
      U hash(const U &index) const
      {
        return hash_policy(index, array.size());
      }
      
      
      ///////////////////////////////////////////////////////////
      //// 1.3 bool insert(U index, const T &t, bool update) ////
      
      //########################################################
      //################# FOR KEY FUNKTION #####################
      bool insert(U index, const T &t, bool update)
      {
#ifndef NDEBUG
        U check = 0;
#endif
        for (U i = index; ; i = (i+1)%array.size()) {
          assert(check++ < array.size());
          // 已存的元素，不在 array的末尾， 更新它就可以了
	  // How to ensure the last element have the lowest energy, the array before insert, it has already been sorted??
          if (init[i]) {
            if (inspector.equal(array(i), t)) {
	      //DEBUG std::cout << t.first.hi_rope<<"["<<t.first.pp<<"]"  << "(" << t.second << ")" << std::endl;
	     
              assert (update);
	      /*
	       * void update(type &dst, const type &src) 
		{
		  if ((src.second < dst.second))
		    {
		      dst.second = src.second;
		    }

		}
	       */
	      
              inspector.update(array(i), t);
	      //append(array(i).first.neighborHelixI, t.first.neighborHelixI);
              return false;
            }
          } else {
            init[i] = true;
	    // 存进去了一个value, used_ 要加 1 的，就是说新存的元素在 array的末尾
	    // if this record is alreay in category, update the energy value
            if (update)
              array.init(i, inspector.init(t));
            else
              array.init(i, t);
            return true;
          }
        }
        assert(0);
      }
 
 
      /////////////////////////////
      //// 1.2 add(const T &t) ////
      
      // #################################################################
      // ## at add-phase, alway check if it is used ==> unique function ##
      void add(const T &t, bool update)
      {
        assert(!finalized);
	// if all allocated space have been used up, expand the space 
        if (!array.size() || loaded())
          expand();
	//DEBUG std::cout << t.first.hi_rope << "(";
        U index = hash(inspector.hash(t));
	//DEBUG std::cout << "inspector.hash(t)=" << inspector.hash(t) << ")" << "(index=" << index << ")" << std::endl;
        bool r = insert(index, t, update);
        if (r)
          ++used_;
      }


      Set(const Set&);
      Set &operator=(const Set&);

    public:
      U ref_count;
      Set()
        : used_(0),
#ifndef NDEBUG
          finalized(false),
#endif
          ref_count(1)
      {
      }
      void resize(U i)
      {
        if (i < array.size())
          return;
        rehash(i);
      }
      
      /////////////////////////////
      //// 1.1 add(const T &t) ////
      void add(const T &t)
      {
        add(t, true);
      }
      
      bool is_empty() const { return !used_; }

      typedef typename Vector_Sparse<T, U>::iterator iterator;

      iterator begin() { assert(finalized); return array.begin(); }
      iterator end() { return array.end(); }


      //////////////////////
      //// 2.1 filter() ////
      void filter()
      {
        if (is_empty())
          return;

        if (inspector.filter()) {
	  //DEBUG std::cout << "inspector.filter()" << std::endl;
          Vector_Sparse<T, U> a(used_);
          swap(array, a);
          std::vector<bool> b(used_);
          swap(init, b);
          used_ = 0;
          for (typename Vector_Sparse<T, U>::iterator i = a.begin();
               i != a.end(); ++i)
            if (!inspector.filter(*i)) {
              array.init(used_, *i);
              init[used_] = true;
              ++used_;
            }
        }

        

        if (!inspector.cutoff() && Shrink_Policy::value) {
          //DEBUG std::cout << "===========================================" << std::endl;
          Vector_Sparse<T, U> a(used_);
          swap(array, a);
          std::vector<bool> b(used_);
          swap(init, b);
          U j = 0;
          for (typename Vector_Sparse<T, U>::iterator i = a.begin();
              i != a.end(); ++i, ++j) {
            array.init(j, *i);
            init[j] = true;
          }
        }

        // ## kbest ##
        if (inspector.cutoff()) {
          //DEBUG std::cout << "inspector.cutoff()" << std::endl;

          assert(Shrink_Policy::value);
          assert(array.end()-array.begin() == used_);
          std::sort(array.begin(), array.end(), typename Inspector::compare());
          if (array.begin() != array.end()) {
            typename Vector_Sparse<T, U>::iterator last = array.begin();
            typename Vector_Sparse<T, U>::iterator i = array.begin();
            ++i;
            U uniques = 1;
            U newend = 1;
            U k = inspector.k();
	    
            if (uniques<k)
            for (; i != array.end(); ++last, ++i) {
// without comments, sometimes more than 100, since it treat the hishapes with the same energy a record.
// with comments, it will precisely output 100 records
#if 0 
              if (!inspector.equal_score(*last, *i))
#endif
                ++uniques;
              if (uniques > k)
                break;
              ++newend;
            }
            used_ = newend;
            Vector_Sparse<T,U> a(used_);
            swap(array, a);
            std::vector<bool> b(used_);
            swap(init, b);
            U j = 0;
            for (typename Vector_Sparse<T, U>::iterator  x = a.begin(); j<newend; ++x, ++j) {
              array.init(j, *x);
              init[j] = true;
            }
          }
        }
      }

      ////////////////////////
      //// 3.1 finalize() ////
      void finalize()
      {
#ifndef NDEBUG
        assert(!finalized);
        finalized = true;
#endif

        if (is_empty())
          return;

        for (typename Vector_Sparse<T, U>::iterator i = array.begin();
            i != array.end(); ++i) {
	 // std::cout << typeid(*i).name() << std::endl;
	 // std::cout << (*i).hi_rope << std::endl;
	 // std::cout << typeid(inspector).name() << std::endl;
	  // ## DO NOTHING in hishape_mfe example ##
          inspector.finalize(*i);
	}
      }

      void *operator new(size_t t) throw (std::bad_alloc);
      void operator delete(void *b) throw ();
  };

#define SET_TEMPLATE_DECL \
  class T, \
  class I, \
  typename U, \
  typename Hash_Policy, \
  template<typename> class Resize_Policy, \
  typename Shrink_Policy, \
  typename Stat_Policy, \
  unsigned int load_factor
#define SET_TEMPLATE_ARGS \
T, I, U, Hash_Policy, Resize_Policy, Shrink_Policy, Stat_Policy, load_factor

  template<SET_TEMPLATE_DECL>
    struct Set_Dummy {
      static Pool<Set<SET_TEMPLATE_ARGS> > pool;
    };

  template<SET_TEMPLATE_DECL>
     Pool<Set<SET_TEMPLATE_ARGS> >
       Set_Dummy<SET_TEMPLATE_ARGS>::pool;

  template<SET_TEMPLATE_DECL>
  void *Set<SET_TEMPLATE_ARGS>::operator new(size_t t) throw (std::bad_alloc)
  {
    assert(sizeof(Set<SET_TEMPLATE_ARGS>) == t);
    Set<SET_TEMPLATE_ARGS> *r = Set_Dummy<SET_TEMPLATE_ARGS>::pool.malloc();
    return r;
  }

  template<SET_TEMPLATE_DECL>
  void Set<SET_TEMPLATE_ARGS>::operator delete(void *b) throw ()
  {
    if (!b)
      return;
    Set_Dummy<SET_TEMPLATE_ARGS>::pool.free(
        static_cast<Set<SET_TEMPLATE_ARGS>*>(b));
  }

  template<class T, class I>
  class Ref : public ::Ref::Lazy<Set<T, I> >
  {
    private:
    public:
  };

#undef SET_TEMPLATE_DECL
#undef SET_TEMPLATE_ARGS

}

#include "empty.hh"

///////////////////////////////////
//// 2.0 hash_filter( answers) ////
template<class T, class I>
inline void hash_filter(Hash::Ref<T, I> &x)
{
  // the function k best of every group is implemented here
  x->filter();
}


template<class T>
inline void finalize(T &x)
{
  //std::cerr << "aa" << std::endl;
}

////////////////////////////////
//// 3.0 finalize( answers) ////
// I is inspector is a class
template<class T, class I>
inline void finalize(Hash::Ref<T, I> &x)
{
  //  std::cerr << "bb" << std::endl;
  // 确定没有做任何事这里，空的，可以被删除。
  x->finalize();
        //for (typename Vector_Sparse<T, U>::iterator i = array.begin();
        //    i != array.end(); ++i)
        //  inspector.finalize(*i);
}

/////////////////////////////////////
//// 1.0 push_back(answers, ans) ////
template<class T, class I>
inline void push_back(Hash::Ref<T, I> &x, const T &e)
{
  assert(is_not_empty(e));
  x->add(e);
}

template<class T, class I>
inline void append(Hash::Ref<T, I> &x, Hash::Ref<T, I> &e)
{
  if (is_empty(e))
    return;
  assert(&x.ref() != &e.ref());
  for (typename Hash::Ref<T, I>::iterator i = e->begin(); i != e->end(); ++i)
    x->add(*i);
}

template<class T, class I, typename Iterator>
inline void append_filter(Hash::Ref<T, I> &x, std::pair<Iterator, Iterator> i)
{
  for (Iterator a = i.first; a != i.second; ++a)
    push_back(x, *a);
  hash_filter(x);
}

template<class T, class I>
inline void empty(Hash::Ref<T, I> &x)
{
}

template<class T, class I>
inline bool is_empty(const Hash::Ref<T, I> &x)
{
  return !x.l || x.const_ref().is_empty();
}

template<class T, class I>
inline void erase(Hash::Ref<T, I> &x)
{
}

template<class T, class U>
inline void update_filter(T &x, const U &a)
{
  x.update(a);
}

template<class T, class I>
inline
std::ostream &operator<<(std::ostream &out, Hash::Ref<T, I> &x)
{
  if (is_empty(x))
    return out;
  typename Hash::Ref<T, I>::Type &h = x.ref();
  for (typename Hash::Ref<T, I>::iterator i = h.begin(); i != h.end(); ++i)
    out << *i << '\n';
  return out;
}


#undef HASH_INITIAL
#undef HASH_LOAD_FACTOR
#undef HASH_SHRINK
#undef HASH_STAT

#endif
