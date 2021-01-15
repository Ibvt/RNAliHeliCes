
// A dynamic programming evaluator generated by GAP-C.
// 
//   GAP-C version:
//     bellmansgapc-2013.05.02
// 
//   GAP-C call:
//     /home/jhuang/local/gapc/bin/gapc -I ../../../ -p alg_ali_mfe * alg_ali_dotBracket_id ../../../eval_ali_microstate.gap -o eval_microstate.cc 
// 
// 


#ifndef eval_microstate_hh
#define eval_microstate_hh

#include "rtlib/adp.hh"

typedef Basic_Subsequence<M_Char, unsigned> TUSubsequence;

typedef Shape shape_t;
#include <rtlib/subopt.hh>
#include "rna.hh"
#include "Extensions/evalfold.hh"
#include "Extensions/alignment.hh"
#include "Extensions/typesRNAfolding.hh"
#include "Extensions/shapes.hh"

#include "Extensions/rnaoptions.hh"

class eval_microstate {

  public:
Basic_Sequence<M_Char> t_0_seq;
unsigned int t_0_left_most;
unsigned int t_0_right_most;

List_Ref<std::pair<mfecovar, String> > LBmfecovar_firstG_string_secondG_EM_zero;
std::pair<mfecovar, String>  Bmfecovar_firstG_string_secondG_E_zero;

class dangle_table_t {

private:

unsigned int t_0_left_most;
unsigned int t_0_right_most;
std::vector<List_Ref<std::pair<mfecovar, String> > > array;
std::vector<bool> tabulated;
unsigned int t_0_n;
List_Ref<std::pair<mfecovar, String> > zero;
unsigned int size()
{
  return (1 * ((((t_0_n * (t_0_n + 1)) / 2) + t_0_n) + 1));
}


public:

dangle_table_t()
{
  empty(zero);
}

void init(unsigned int t_0_n_, const std::string &tname)
{
t_0_n = t_0_n_;
t_0_left_most = 0;
t_0_right_most = t_0_n;
unsigned int newsize = size();
array.resize(newsize);
tabulated.clear();
tabulated.resize(newsize);
}
bool is_tabulated(unsigned int t_0_i, unsigned int t_0_j)
{
  if (((t_0_j - t_0_i) < 5))
    {
      return true;
    }

  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  return tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))];
}


void clear() { tabulated.clear(); }
List_Ref<std::pair<mfecovar, String> > &  get(unsigned int t_0_i, unsigned int t_0_j)
{
  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  if (((t_0_j - t_0_i) < 5))
    {
      return zero;
    }

  assert( tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))]);
  assert( ((0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i))) < size()));
  return array[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))];
}


void set(unsigned int t_0_i, unsigned int t_0_j, List_Ref<std::pair<mfecovar, String> > e)
{
  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  if (((t_0_j - t_0_i) < 5))
    {
      assert( 0);
    }

  assert( !is_tabulated(t_0_i, t_0_j));
  assert( ((0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i))) < size()));
  array[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))] = e;
  tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))] = true;
}


};
dangle_table_t dangle_table;


class iloop_table_t {

private:

unsigned int t_0_left_most;
unsigned int t_0_right_most;
std::vector<List_Ref<std::pair<mfecovar, String> > > array;
std::vector<bool> tabulated;
unsigned int t_0_n;
List_Ref<std::pair<mfecovar, String> > zero;
unsigned int size()
{
  return (1 * ((((t_0_n * (t_0_n + 1)) / 2) + t_0_n) + 1));
}


public:

iloop_table_t()
{
  empty(zero);
}

void init(unsigned int t_0_n_, const std::string &tname)
{
t_0_n = t_0_n_;
t_0_left_most = 0;
t_0_right_most = t_0_n;
unsigned int newsize = size();
array.resize(newsize);
tabulated.clear();
tabulated.resize(newsize);
}
bool is_tabulated(unsigned int t_0_i, unsigned int t_0_j)
{
  if (((t_0_j - t_0_i) < 9))
    {
      return true;
    }

  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  return tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))];
}


void clear() { tabulated.clear(); }
List_Ref<std::pair<mfecovar, String> > &  get(unsigned int t_0_i, unsigned int t_0_j)
{
  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  if (((t_0_j - t_0_i) < 9))
    {
      return zero;
    }

  assert( tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))]);
  assert( ((0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i))) < size()));
  return array[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))];
}


void set(unsigned int t_0_i, unsigned int t_0_j, List_Ref<std::pair<mfecovar, String> > e)
{
  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  if (((t_0_j - t_0_i) < 9))
    {
      assert( 0);
    }

  assert( !is_tabulated(t_0_i, t_0_j));
  assert( ((0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i))) < size()));
  array[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))] = e;
  tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))] = true;
}


};
iloop_table_t iloop_table;


class ml_comps_table_t {

private:

unsigned int t_0_left_most;
unsigned int t_0_right_most;
std::vector<List_Ref<std::pair<mfecovar, String> > > array;
std::vector<bool> tabulated;
unsigned int t_0_n;
List_Ref<std::pair<mfecovar, String> > zero;
unsigned int size()
{
  return (1 * ((((t_0_n * (t_0_n + 1)) / 2) + t_0_n) + 1));
}


public:

ml_comps_table_t()
{
  empty(zero);
}

void init(unsigned int t_0_n_, const std::string &tname)
{
t_0_n = t_0_n_;
t_0_left_most = 0;
t_0_right_most = t_0_n;
unsigned int newsize = size();
array.resize(newsize);
tabulated.clear();
tabulated.resize(newsize);
}
bool is_tabulated(unsigned int t_0_i, unsigned int t_0_j)
{
  if (((t_0_j - t_0_i) < 10))
    {
      return true;
    }

  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  return tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))];
}


void clear() { tabulated.clear(); }
List_Ref<std::pair<mfecovar, String> > &  get(unsigned int t_0_i, unsigned int t_0_j)
{
  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  if (((t_0_j - t_0_i) < 10))
    {
      return zero;
    }

  assert( tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))]);
  assert( ((0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i))) < size()));
  return array[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))];
}


void set(unsigned int t_0_i, unsigned int t_0_j, List_Ref<std::pair<mfecovar, String> > e)
{
  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  if (((t_0_j - t_0_i) < 10))
    {
      assert( 0);
    }

  assert( !is_tabulated(t_0_i, t_0_j));
  assert( ((0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i))) < size()));
  array[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))] = e;
  tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))] = true;
}


};
ml_comps_table_t ml_comps_table;


class ml_comps1_table_t {

private:

unsigned int t_0_left_most;
unsigned int t_0_right_most;
std::vector<List_Ref<std::pair<mfecovar, String> > > array;
std::vector<bool> tabulated;
unsigned int t_0_n;
List_Ref<std::pair<mfecovar, String> > zero;
unsigned int size()
{
  return (1 * ((((t_0_n * (t_0_n + 1)) / 2) + t_0_n) + 1));
}


public:

ml_comps1_table_t()
{
  empty(zero);
}

void init(unsigned int t_0_n_, const std::string &tname)
{
t_0_n = t_0_n_;
t_0_left_most = 0;
t_0_right_most = t_0_n;
unsigned int newsize = size();
array.resize(newsize);
tabulated.clear();
tabulated.resize(newsize);
}
bool is_tabulated(unsigned int t_0_i, unsigned int t_0_j)
{
  if (((t_0_j - t_0_i) < 5))
    {
      return true;
    }

  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  return tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))];
}


void clear() { tabulated.clear(); }
List_Ref<std::pair<mfecovar, String> > &  get(unsigned int t_0_i, unsigned int t_0_j)
{
  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  if (((t_0_j - t_0_i) < 5))
    {
      return zero;
    }

  assert( tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))]);
  assert( ((0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i))) < size()));
  return array[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))];
}


void set(unsigned int t_0_i, unsigned int t_0_j, List_Ref<std::pair<mfecovar, String> > e)
{
  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  if (((t_0_j - t_0_i) < 5))
    {
      assert( 0);
    }

  assert( !is_tabulated(t_0_i, t_0_j));
  assert( ((0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i))) < size()));
  array[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))] = e;
  tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))] = true;
}


};
ml_comps1_table_t ml_comps1_table;


class strong_table_t {

private:

unsigned int t_0_left_most;
unsigned int t_0_right_most;
std::vector<List_Ref<std::pair<mfecovar, String> > > array;
std::vector<bool> tabulated;
unsigned int t_0_n;
List_Ref<std::pair<mfecovar, String> > zero;
unsigned int size()
{
  return (1 * ((((t_0_n * (t_0_n + 1)) / 2) + t_0_n) + 1));
}


public:

strong_table_t()
{
  empty(zero);
}

void init(unsigned int t_0_n_, const std::string &tname)
{
t_0_n = t_0_n_;
t_0_left_most = 0;
t_0_right_most = t_0_n;
unsigned int newsize = size();
array.resize(newsize);
tabulated.clear();
tabulated.resize(newsize);
}
bool is_tabulated(unsigned int t_0_i, unsigned int t_0_j)
{
  if (((t_0_j - t_0_i) < 5))
    {
      return true;
    }

  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  return tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))];
}


void clear() { tabulated.clear(); }
List_Ref<std::pair<mfecovar, String> > &  get(unsigned int t_0_i, unsigned int t_0_j)
{
  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  if (((t_0_j - t_0_i) < 5))
    {
      return zero;
    }

  assert( tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))]);
  assert( ((0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i))) < size()));
  return array[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))];
}


void set(unsigned int t_0_i, unsigned int t_0_j, List_Ref<std::pair<mfecovar, String> > e)
{
  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  if (((t_0_j - t_0_i) < 5))
    {
      assert( 0);
    }

  assert( !is_tabulated(t_0_i, t_0_j));
  assert( ((0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i))) < size()));
  array[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))] = e;
  tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))] = true;
}


};
strong_table_t strong_table;


class struct_table_t {

private:

unsigned int t_0_left_most;
unsigned int t_0_right_most;
std::vector<List_Ref<std::pair<mfecovar, String> > > array;
std::vector<bool> tabulated;
unsigned int t_0_n;
List_Ref<std::pair<mfecovar, String> > zero;
unsigned int size()
{
  return (1 * ((t_0_n + 1) * 1));
}


public:

struct_table_t()
{
  empty(zero);
}

void init(unsigned int t_0_n_, const std::string &tname)
{
t_0_n = t_0_n_;
t_0_left_most = 0;
t_0_right_most = t_0_n;
unsigned int newsize = size();
array.resize(newsize);
tabulated.clear();
tabulated.resize(newsize);
}
bool is_tabulated(unsigned int t_0_i)
{
  unsigned int t_0_j = t_0_n;
  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  if ((t_0_j < (t_0_n - 0)))
    {
      return true;
    }

  unsigned int t_0_real_j = (t_0_n - t_0_j);
  return tabulated[(0 + (1 * (t_0_i + (t_0_real_j * (t_0_n + 1)))))];
}


void clear() { tabulated.clear(); }
List_Ref<std::pair<mfecovar, String> > &  get(unsigned int t_0_i)
{
  unsigned int t_0_j = t_0_n;
  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  if ((t_0_j < (t_0_n - 0)))
    {
      return zero;
    }

  unsigned int t_0_real_j = (t_0_n - t_0_j);
  assert( tabulated[(0 + (1 * (t_0_i + (t_0_real_j * (t_0_n + 1)))))]);
  assert( ((0 + (1 * (t_0_i + (t_0_real_j * (t_0_n + 1))))) < size()));
  return array[(0 + (1 * (t_0_i + (t_0_real_j * (t_0_n + 1)))))];
}


void set(unsigned int t_0_i, List_Ref<std::pair<mfecovar, String> > e)
{
  unsigned int t_0_j = t_0_n;
  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  if ((t_0_j < (t_0_n - 0)))
    {
      return;
    }

  unsigned int t_0_real_j = (t_0_n - t_0_j);
  assert( !is_tabulated(t_0_i));
  assert( ((0 + (1 * (t_0_i + (t_0_real_j * (t_0_n + 1))))) < size()));
  array[(0 + (1 * (t_0_i + (t_0_real_j * (t_0_n + 1)))))] = e;
  tabulated[(0 + (1 * (t_0_i + (t_0_real_j * (t_0_n + 1)))))] = true;
}


};
struct_table_t struct_table;


class weak_table_t {

private:

unsigned int t_0_left_most;
unsigned int t_0_right_most;
std::vector<List_Ref<std::pair<mfecovar, String> > > array;
std::vector<bool> tabulated;
unsigned int t_0_n;
List_Ref<std::pair<mfecovar, String> > zero;
unsigned int size()
{
  return (1 * ((((t_0_n * (t_0_n + 1)) / 2) + t_0_n) + 1));
}


public:

weak_table_t()
{
  empty(zero);
}

void init(unsigned int t_0_n_, const std::string &tname)
{
t_0_n = t_0_n_;
t_0_left_most = 0;
t_0_right_most = t_0_n;
unsigned int newsize = size();
array.resize(newsize);
tabulated.clear();
tabulated.resize(newsize);
}
bool is_tabulated(unsigned int t_0_i, unsigned int t_0_j)
{
  if (((t_0_j - t_0_i) < 5))
    {
      return true;
    }

  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  return tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))];
}


void clear() { tabulated.clear(); }
List_Ref<std::pair<mfecovar, String> > &  get(unsigned int t_0_i, unsigned int t_0_j)
{
  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  if (((t_0_j - t_0_i) < 5))
    {
      return zero;
    }

  assert( tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))]);
  assert( ((0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i))) < size()));
  return array[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))];
}


void set(unsigned int t_0_i, unsigned int t_0_j, List_Ref<std::pair<mfecovar, String> > e)
{
  assert( (t_0_i <= t_0_j));
  assert( (t_0_j <= t_0_n));
  if (((t_0_j - t_0_i) < 5))
    {
      assert( 0);
    }

  assert( !is_tabulated(t_0_i, t_0_j));
  assert( ((0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i))) < size()));
  array[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))] = e;
  tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))] = true;
}


};
weak_table_t weak_table;



void init(const gapc::Opts &opts)
{
const std::vector<std::pair<const char *, unsigned> > &inp = opts.inputs;
if(inp.size() != 1)
  throw gapc::OptException("Number of input sequences does not match.");

  t_0_seq.copy(inp[0].first, inp[0].second);
char_to_rna(t_0_seq);
  dangle_table.init( t_0_seq.size(), "dangle_table");
  iloop_table.init( t_0_seq.size(), "iloop_table");
  ml_comps_table.init( t_0_seq.size(), "ml_comps_table");
  ml_comps1_table.init( t_0_seq.size(), "ml_comps1_table");
  strong_table.init( t_0_seq.size(), "strong_table");
  struct_table.init( t_0_seq.size(), "struct_table");
  weak_table.init( t_0_seq.size(), "weak_table");
empty(LBmfecovar_firstG_string_secondG_EM_zero);
empty(Bmfecovar_firstG_string_secondG_E_zero);

t_0_left_most = 0;
t_0_right_most = t_0_seq.size();
}

  private:
    List_Ref<std::pair<mfecovar, String> > &  nt_dangle(unsigned int t_0_i, unsigned int t_0_j);
    std::pair<mfecovar, String>  nt_hairpin(unsigned int t_0_i, unsigned int t_0_j);
    List_Ref<std::pair<mfecovar, String> > &  nt_iloop(unsigned int t_0_i, unsigned int t_0_j);
    List_Ref<std::pair<mfecovar, String> > nt_leftB(unsigned int t_0_i, unsigned int t_0_j);
    List_Ref<std::pair<mfecovar, String> > &  nt_ml_comps(unsigned int t_0_i, unsigned int t_0_j);
    List_Ref<std::pair<mfecovar, String> > &  nt_ml_comps1(unsigned int t_0_i, unsigned int t_0_j);
    List_Ref<std::pair<mfecovar, String> > nt_multiloop(unsigned int t_0_i, unsigned int t_0_j);
    List_Ref<std::pair<mfecovar, String> > nt_rightB(unsigned int t_0_i, unsigned int t_0_j);
    List_Ref<std::pair<mfecovar, String> > nt_stack(unsigned int t_0_i, unsigned int t_0_j);
    List_Ref<std::pair<mfecovar, String> > &  nt_strong(unsigned int t_0_i, unsigned int t_0_j);
    List_Ref<std::pair<mfecovar, String> > &  nt_struct(unsigned int t_0_i);
    List_Ref<std::pair<mfecovar, String> > &  nt_weak(unsigned int t_0_i, unsigned int t_0_j);

    std::pair<mfecovar, String>  addss(const std::pair<mfecovar, String> & p_x, const TUSubsequence & p_r);
    std::pair<mfecovar, String>  bl(const TUSubsequence & p_lb, const TUSubsequence & p_lr, const std::pair<mfecovar, String> & p_x, const TUSubsequence & p_rb);
    std::pair<mfecovar, String>  br(const TUSubsequence & p_lb, const std::pair<mfecovar, String> & p_x, const TUSubsequence & p_rr, const TUSubsequence & p_rb);
    std::pair<mfecovar, String>  cadd(const std::pair<mfecovar, String> & p_x, const std::pair<mfecovar, String> & p_y);
    std::pair<mfecovar, String>  drem(const TUSubsequence & p_lb, const std::pair<mfecovar, String> & p_x, const TUSubsequence & p_rb);
    std::pair<mfecovar, String>  edl(const TUSubsequence & p_ldangle, const std::pair<mfecovar, String> & p_x, const TUSubsequence & p_rb);
    std::pair<mfecovar, String>  edlr(const TUSubsequence & p_ldangle, const std::pair<mfecovar, String> & p_x, const TUSubsequence & p_rdangle);
    std::pair<mfecovar, String>  edr(const TUSubsequence & p_lb, const std::pair<mfecovar, String> & p_x, const TUSubsequence & p_rdangle);
    List_Ref<std::pair<mfecovar, String> > h(List_Ref<std::pair<mfecovar, String> > i);
    std::pair<mfecovar, String>  hl(const TUSubsequence & p_lb, const TUSubsequence & p_r, const TUSubsequence & p_rb);
    std::pair<mfecovar, String>  il(const TUSubsequence & p_lb, const TUSubsequence & p_lr, const std::pair<mfecovar, String> & p_x, const TUSubsequence & p_rr, const TUSubsequence & p_rb);
    std::pair<mfecovar, String>  incl(const std::pair<mfecovar, String> & p_x);
    std::pair<mfecovar, String>  ml(const TUSubsequence & p_lb, const std::pair<mfecovar, String> & p_x, const TUSubsequence & p_rb);
    std::pair<mfecovar, String>  mldl(const TUSubsequence & p_lb, const TUSubsequence & p_dl, const std::pair<mfecovar, String> & p_x, const TUSubsequence & p_rb);
    std::pair<mfecovar, String>  mldlr(const TUSubsequence & p_lb, const TUSubsequence & p_dl, const std::pair<mfecovar, String> & p_x, const TUSubsequence & p_dr, const TUSubsequence & p_rb);
    std::pair<mfecovar, String>  mldr(const TUSubsequence & p_lb, const std::pair<mfecovar, String> & p_x, const TUSubsequence & p_dr, const TUSubsequence & p_rb);
    std::pair<mfecovar, String>  nil(const TUSubsequence & p_n);
    std::pair<mfecovar, String>  sadd(const TUSubsequence & p_lb, const std::pair<mfecovar, String> & p_x);
    std::pair<mfecovar, String>  sr(const TUSubsequence & p_lb, const std::pair<mfecovar, String> & p_x, const TUSubsequence & p_rb);


    mfecovar addss_l(const mfecovar & x, const TUSubsequence & r);
    mfecovar bl_l(const TUSubsequence & lb, const TUSubsequence & lr, const mfecovar & x, const TUSubsequence & rb);
    mfecovar br_l(const TUSubsequence & lb, const mfecovar & x, const TUSubsequence & rr, const TUSubsequence & rb);
    mfecovar cadd_l(const mfecovar & x, const mfecovar & y);
    mfecovar drem_l(const TUSubsequence & lb, const mfecovar & x, const TUSubsequence & rb);
    mfecovar edl_l(const TUSubsequence & ldangle, const mfecovar & x, const TUSubsequence & rb);
    mfecovar edlr_l(const TUSubsequence & ldangle, const mfecovar & x, const TUSubsequence & rdangle);
    mfecovar edr_l(const TUSubsequence & lb, const mfecovar & x, const TUSubsequence & rdangle);
    mfecovar h_l(List_Ref<mfecovar> i);
template <typename Iterator>
    mfecovar h_l(std::pair<Iterator, Iterator> i)
;
    mfecovar hl_l(const TUSubsequence & lb, const TUSubsequence & r, const TUSubsequence & rb);
    mfecovar il_l(const TUSubsequence & lb, const TUSubsequence & lr, const mfecovar & x, const TUSubsequence & rr, const TUSubsequence & rb);
    mfecovar incl_l(const mfecovar & x);
    mfecovar ml_l(const TUSubsequence & lb, const mfecovar & x, const TUSubsequence & rb);
    mfecovar mldl_l(const TUSubsequence & lb, const TUSubsequence & dl, const mfecovar & x, const TUSubsequence & rb);
    mfecovar mldlr_l(const TUSubsequence & lb, const TUSubsequence & dl, const mfecovar & x, const TUSubsequence & dr, const TUSubsequence & rb);
    mfecovar mldr_l(const TUSubsequence & lb, const mfecovar & x, const TUSubsequence & dr, const TUSubsequence & rb);
    mfecovar nil_l(const TUSubsequence & n);
    mfecovar sadd_l(const TUSubsequence & lb, const mfecovar & x);
    mfecovar sr_l(const TUSubsequence & lb, const mfecovar & x, const TUSubsequence & rb);


    String addss_r(const String & e, const TUSubsequence & rb);
    String bl_r(const TUSubsequence & lb, const TUSubsequence & lregion, const String & e, const TUSubsequence & rb);
    String br_r(const TUSubsequence & lb, const String & e, const TUSubsequence & rregion, const TUSubsequence & rb);
    String cadd_r(const String & le, const String & re);
    String drem_r(const TUSubsequence & lloc, const String & e, const TUSubsequence & rloc);
    String edl_r(const TUSubsequence & lb, const String & e, const TUSubsequence & loc);
    String edlr_r(const TUSubsequence & lb, const String & e, const TUSubsequence & rb);
    String edr_r(const TUSubsequence & loc, const String & e, const TUSubsequence & rb);
    List_Ref<String> h_r(List_Ref<String> i);
    String hl_r(const TUSubsequence & lb, const TUSubsequence & region, const TUSubsequence & rb);
    String il_r(const TUSubsequence & lb, const TUSubsequence & lregion, const String & e, const TUSubsequence & rregion, const TUSubsequence & rb);
    String incl_r(const String & e);
    String ml_r(const TUSubsequence & lb, const String & e, const TUSubsequence & rb);
    String mldl_r(const TUSubsequence & lb, const TUSubsequence & dl, const String & e, const TUSubsequence & rb);
    String mldlr_r(const TUSubsequence & lb, const TUSubsequence & dl, const String & e, const TUSubsequence & dr, const TUSubsequence & rb);
    String mldr_r(const TUSubsequence & lb, const String & e, const TUSubsequence & dr, const TUSubsequence & rb);
    String nil_r(const TUSubsequence & loc);
    String sadd_r(const TUSubsequence & lb, const String & e);
    String sr_r(const TUSubsequence & lb, const String & e, const TUSubsequence & rb);


 public:
   void cyk();

 public:
   List_Ref<std::pair<mfecovar, String> > &  run()
{
  return nt_struct(t_0_left_most);
}
void print_stats(std::ostream &o)
{
#ifdef STATS
      o << "\n\nN = " << seq.size() << '\n';
      dangle_table.print_stats(o, "dangle_table");
      iloop_table.print_stats(o, "iloop_table");
      ml_comps_table.print_stats(o, "ml_comps_table");
      ml_comps1_table.print_stats(o, "ml_comps1_table");
      strong_table.print_stats(o, "strong_table");
      struct_table.print_stats(o, "struct_table");
      weak_table.print_stats(o, "weak_table");
#endif
}

template <typename Value>   void  print_result(std::ostream &out, Value& res)

{
if (isEmpty(res))
  out << "[]\n";
else
  out << res << '\n';

}
template <typename Value>   void  print_backtrack(std::ostream &out, Value& value)

{
}
   void  print_subopt(std::ostream &out, int  delta = 0) {}

};

#ifndef NO_GAPC_TYPEDEFS
namespace gapc {
  typedef eval_microstate class_name;
  typedef List_Ref<std::pair<mfecovar, String> > &  return_type;
}
#endif

#endif
