#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <set>
#include <limits>
#include <cassert>
#include <math.h>

// #include <dlfcn.h>
// #include <execinfo.h>
// #include <typeinfo>
// #include <string>
// #include <memory>
// #include <cxxabi.h>
// #include <cstdlib>
// 
#include "hishapeh_mfe_pp.hh"
#include "hishapehplus_mfe_pp.hh"
#include "hishapem_mfe_pp.hh"
#include "hishapeb_mfe_pp.hh"

#include "rtlib/string.hh"
#include "rtlib/list.hh"
#include "rtlib/hash.hh"
#include "rtlib/asymptotics.hh"
#include "rtlib/generic_opts.hh"

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/replace.hpp>

#include "hited_util.hh"

namespace po = boost::program_options;

//using namespace std;




#define GRAPHSIZE 4096
#define MAX(a, b) ((a > b) ? (a) : (b))
#define HISHREPLENGTH 2048


//int e; /* The number of nonzero edges in the graph */

int n_node; /* The number of nodes in the graph */
float distances[GRAPHSIZE][GRAPHSIZE]; /* distances[i][j] is the distance between node i and j; or 0 if there is no direct connection */
float d[GRAPHSIZE]; /* d[i] is the length of the shortest path between the source (s) and node i */
int prev[GRAPHSIZE]; /* prev[i] is the node that comes right before i in the shortest path from the source to i*/
int anchors[100];



// ###########################################################################################
// ############ c functions from findpath.h, fold.h, fold_vars.h and utils.h #################
typedef struct path {
  double en;
  char *s;
} path_t;

extern "C" int find_saddle (char *seq, char *struc1, char *struc2, int max);
extern "C" path_t* get_path(char *seq, char *s1, char* s2, int maxkeep);


/* function from fold.c */
extern "C" float  fold(const char *sequence, char *structure);
/* calculate mfe-structure of sequence */
extern "C" float  energy_of_struct(const char *string, const char *structure);
/* calculate energy of string on structure */
extern "C" void   free_arrays(void);           /* free arrays for mfe folding */
extern "C" void   initialize_fold(int length); /* allocate arrays for folding */
extern "C" void   update_fold_params(void);    /* recalculate parameters */
extern "C" char  *backtrack_fold_from_pair(char *sequence, int i, int j);
extern "C" int loop_energy(short * ptable, short *s, short *s1, int i);
extern "C" void		export_fold_arrays(int **f5_p, int **c_p, int **fML_p, int **fM1_p, int **indx_p, char **ptype_p);

/* some circfold related functions...	*/
extern "C"	float	circfold(const char *string, char *structure);
extern "C"	float	energy_of_circ_struct(const char *string, const char *structure);
extern "C"	void	export_circfold_arrays(int *Fc_p, int *FcH_p, int *FcI_p, int *FcM_p, int **fM2_p, int **f5_p, int **c_p, int **fML_p, int **fM1_p, int **indx_p, char **ptype_p);




/* to use floats instead of doubles in pf_fold() comment next line */
#define LARGE_PF
#ifdef  LARGE_PF
#define FLT_OR_DBL double
#else
#define FLT_OR_DBL float
#endif

extern "C" int  noGU;           /* GU not allowed at all */
extern "C" int  no_closingGU;   /* GU allowed only inside stacks */
extern "C" int  tetra_loop;     /* Fold with specially stable 4-loops */
extern "C" int  energy_set;     /* 0 = BP; 1=any mit GC; 2=any mit AU-parameter */
extern "C" int  dangles;	    /* use dangling end energies (not in part_func!) */
/*@null@*/

extern "C" int oldAliEn;        /* use old alifold energies (with gaps) */
extern "C" int ribo;            /* use ribosum matrices */
extern "C" char *RibosumFile;   /* warning this variable will vanish in the future
			       ribosums will be compiled in instead */
extern "C" char *nonstandards;  /* contains allowed non standard bases */
extern "C" double temperature;   /* rescale parameters to this temperature */
extern "C" int  james_rule;     /* interior loops of size 2 get energy 0.8Kcal and
			       no mismatches, default 1 */
extern "C" int  logML;          /* use logarithmic multiloop energy function */
extern "C" int  cut_point;      /* first position of 2nd strand for co-folding */

typedef struct bond {               /* base pair */
   int i;
   int j;
} bondT;
extern "C" bondT  *base_pair; /* list of base pairs */

extern "C" FLT_OR_DBL *pr;          /* base pairing prob. matrix */
extern "C" int   *iindx;            /* pr[i,j] -> pr[iindx[i]-j] */
extern "C" double pf_scale;         /* scaling factor to avoid float overflows*/
extern "C" int    fold_constrained; /* fold with constraints */
extern "C" int    do_backtrack;     /* calculate pair prob matrix in part_func() */
extern "C" int    noLonelyPairs;    /* avoid helices of length 1 */
extern "C" char backtrack_type;     /* usually 'F'; 'C' require (1,N) to be bonded;
				   'M' seq is part of a multi loop */
char * option_string(void);



/* Header file for utils.c */
#ifdef HAVE_CONFIG_H
#include <config.h>
#ifndef HAVE_STRDUP
extern "C" char *strdup(const char *s);
#endif
#endif
#ifdef WITH_DMALLOC
/* use dmalloc library to check for memory management bugs */
#include "dmalloc.h"
#define space(S) calloc(1,(S))
#else
extern "C" /*@only@*/ /*@notnull@*/
void  *space(unsigned size) /*@ensures MaxSet(result) == (size-1);@*/;
			    /* allocate space safely */
extern "C" /*@only@*/ /*@notnull@*/
void  *xrealloc(/*@null@*/ /*@only@*/ /*@out@*/ /*@returned@*/ void *p, unsigned size) /*@modifies *p @*/ /*@ensures MaxSet(result) == (size-1) @*/;
#endif

extern "C" /*@exits@*/ void nrerror(const char message[]);  /* die with error message */
extern "C" void   init_rand(void);                /* make random number seeds */
extern "C" unsigned short xsubi[3];               /* current 48bit random number */
extern "C" double urn(void);                      /* random number from [0..1] */
extern "C" int    int_urn(int from, int to);      /* random integer */
extern "C" void   filecopy(FILE *from, FILE *to); /* inefficient `cp' */
extern "C" /*@observer@*/ char  *time_stamp(void);               /* current date in a string */
extern "C" /*@only@*/ /*@notnull@*/ char  *random_string(int l, const char symbols[]);
/* random string of length l using characters from symbols[] */
extern "C" int    hamming(const char *s1, const char *s2);
/* calculate hamming distance */
extern "C" /*@only@*/ /*@null@*/ char  *get_line(const FILE *fp); /* read one (arbitrary length) line from fp */


extern "C" char *pack_structure(const char *struc);
/* pack secondary secondary structure, 5:1 compression using base 3 encoding */
extern "C" char *unpack_structure(const char *packed);
/* unpack sec structure packed with pack_structure() */
extern "C" short *make_pair_table(const char *structure);
/* returns a newly allocated table, such that:  table[i]=j if (i.j) pair or
   0 if i is unpaired, table[0] contains the length of the structure. */

extern "C" int bp_distance(const char *str1, const char *str2);
/* dist = {number of base pairs in one structure but not in the other}
   same as edit distance with open-pair close-pair as move-set */

// ################################## end ####################################################
// namespace {
//   void * last_frames[20];
//   size_t last_size;
//   std::string exception_name;
// 
//   std::string demangle(const char *name) {
//     int status;
//     std::unique_ptr<char,void(*)(void*)> realname(abi::__cxa_demangle(name, 0, 0, &status), &std::free);
//     return status ? "failed" : &*realname;
//   }
// }
// 
// extern "C" {
//   void __cxa_throw(void *ex, void *info, void (*dest)(void *)) {
//     exception_name = demangle(reinterpret_cast<const std::type_info*>(info)->name());
//     last_size = backtrace(last_frames, sizeof last_frames/sizeof(void*));
// 
//     static void (*const rethrow)(void*,void*,void(*)(void*)) __attribute__ ((noreturn)) = (void (*)(void*,void*,void(*)(void*)))dlsym(RTLD_NEXT, "__cxa_throw");
//     rethrow(ex,info,dest);
//   }
// }







typedef struct move {
    int moveI;  
    int moveJ;
    int moveE;
    int moveWhen;  
} moveT;

int getBasePairDistance(char *ss1, char *ss2) {
    short *pairTable1, *pairTable2;
    moveT *moveTList;
    int i, length, distance=0;
    pairTable1 = make_pair_table(ss1);
    pairTable2 = make_pair_table(ss2);
    length = (int) strlen(ss1);
    moveTList = (moveT *) space(sizeof(moveT)*length); 

    for ( i=1; i<=length; i++ ) {
        if (pairTable1[i] != pairTable2[i]) {
            if ( i<pairTable1[i] ) {
	        moveTList[distance].moveI = -i;
	        moveTList[distance].moveJ = -pairTable1[i];
	        moveTList[distance++].moveWhen = 0;
            }
            if ( i<pairTable2[i] ) {
	        moveTList[distance].moveI = i;
	        moveTList[distance].moveJ = pairTable2[i];
	        moveTList[distance++].moveWhen = 0;
            }
       }
  }
  free(pairTable1);
  free(pairTable2);
  free(moveTList);
  return distance;
}

void printD() {
	int i;

	printf("Distances:\n");
	for (i = 0; i < n_node; i++)
		printf("%10d", i);
	printf("\n");
	for (i = 0; i < n_node; i++) {
		printf("%10f", d[i]);
	}
	printf("\n");
}

/*
 * Prints the shortest path from the source to dest.
 *
 * dijkstra_all_targets(int) MUST be run at least once BEFORE
 * this is called
 */
void printPathIntoFile(FILE *file, int dest) {
	if (prev[dest] != -1)
		printPathIntoFile(file, prev[dest]);
	fprintf(file, "%d ", dest);
}
void printPath(int dest) {
	if (prev[dest] != -1)
		printPath(prev[dest]);
	printf("%d ", dest);
}


void dijkstra_all_targets(int s) {
	int i, k, mini;
	int visited[GRAPHSIZE];

	// initialize d[], prev[], visited[]
	// empty the three arrays
	for (i = 0; i < n_node; ++i) {
		d[i] = -INFINITY;
		prev[i] = -1; /* no path has yet been found to i */
        	visited[i] = 0; /* the i-th element has not yet been visited */
	}
	//d[s] = 0;



	d[s] = -INFINITY;
	prev[s] = -1;
	for (i = 0; i < n_node; i++) {
	    if (i==s)
	    {
		continue;
	    }
            d[i] = distances[s][i];
	    prev[i] = s;
	    //printf("Path to %f[0][%d]: ", distances[0][i], i);
	}
	visited[s] = 1;
	
	
	for (k = 0; k < n_node; ++k) {
	        if (k==s)
		{
		  continue;
		}
		mini = -1;
		for (i = 0; i < n_node; ++i)
			if (!visited[i] && ((mini == -1) || (d[i] < d[mini])))
				mini = i;
  
		visited[mini] = 1;
		//printf("%d\t",mini);
  
		for (i = 0; i < n_node; ++i)
		{
			if (i==s)
			{
			    continue;
			}
			if (distances[mini][i])
			{
				if (MAX(d[mini],distances[mini][i]) < d[i])
				{
					d[i] = MAX(d[mini],distances[mini][i]);
					prev[i] = mini;
				}
/*
				if (d[mini] + dist[mini][i] < d[i]) {
					d[i] = d[mini] + dist[mini][i];
					prev[i] = mini;
				}
*/
				/*
				float barrier_energy_through_i = (d[i] > distances[k][i]) ? d[i] : distances[k][i];
				if (barrier_energy_through_i < d[k]) {  // if the higher is on the lowest path
					d[k] = barrier_energy_through_i;// set the lowest path on array
					prev[k] = i;
					printf("prev[%d]=%d", k, i);
				}*/
			}
		}
	}
	
}

float calculateDistance(const std::string& sequence, const std::vector<std::string>& structures, int i, int j, int maxkeep)
{

      int k_max=0;
      float en_max=-INFINITY;
      char struc_max[HISHREPLENGTH];
      int dist;
      path_t *route, *r;
      int k;

      route=get_path( const_cast<char*>(sequence.c_str()), const_cast<char*>(structures[i].c_str()), const_cast<char*>(structures[j].c_str()), maxkeep);
      dist=getBasePairDistance(const_cast<char*>(structures[i].c_str()), const_cast<char*>(structures[j].c_str()));

      for (k=0;k<dist+1;k++){
	  if (route[k].en > en_max)
	  {
	      k_max = k;
	      en_max = route[k].en;
	      strcpy(struc_max, route[k].s);
	  }
      }

      for (r=route; r->s; r++) {
	//printf("%s %6.2f - %6.2f\n", r->s, energy_of_struct(seq,r->s), r->en);
	free(r->s);
      }
      free(route);
      distances[i][j] = en_max;
      return distances[i][j];
}
int dijkstra_single_target(const std::string& sequence, const std::vector<std::string>& structures, int s, int t, int maxkeep) {
	int i, k, mini;
	int visited[GRAPHSIZE];
	int dist;

	// initialize d[], prev[], visited[]
	// empty the three arrays
	for (i = 0; i < n_node; ++i) {
		d[i] = -INFINITY;
		prev[i] = -1; /* no path has yet been found to i */
        	visited[i] = 0; /* the i-th element has not yet been visited */
	}
	//d[s] = 0;


        // ## first round calculation ##
	d[s] = -INFINITY;
	prev[s] = -1;
	for (i = 0; i < n_node; i++) {
	    if (i==s)
	    {
		continue;
	    }
	    // =======================================================================
	    d[i] = calculateDistance(sequence, structures, 0, i, maxkeep);    
	    
	    // =======================================================================

            //d[i] = distances[0][i];
	    prev[i] = s;
	    //$$ printf("Path to %f[0][%d]: ", distances[0][i], i);
	}
	visited[s] = 1;
	
	
	for (k = 0; k < n_node; ++k) {
	        if (k==s)
		{
		  continue;
		}
		mini = -1;
		for (i = 0; i < n_node; ++i)
			if (!visited[i] && ((mini == -1) || (d[i] < d[mini])))
				mini = i;
  
		visited[mini] = 1;
		//printf("%d\t",mini);
		
		
		for (i = 0; i < n_node; ++i)
		{
			if (i==s)
			{
			    continue;
			}
			
                        //float distance_mini_i = distances[mini][i];
			//if (distance_mini_i == -INFINITY)
			//{
			      float distance_mini_i = calculateDistance(sequence, structures, mini, i, maxkeep);
			//}
			if (MAX(d[mini],distance_mini_i) < d[i])
			{
				d[i] = MAX(d[mini],distance_mini_i);
				prev[i] = mini;
			}
		}
		
		if (mini==t)
		{
			int anchor_index = 0;
			
			int x=t;
			while (prev[x] != -1)
			{
			    //printf("%d ", x);
			    //printf("%s ", strings[x]);
			    anchors[anchor_index] = x;
			    anchor_index++;
			    x = prev[x];
			}
			
			// add the start item
			anchors[anchor_index] = s;
			anchor_index++;
			return anchor_index;
		}

	}
	
}


struct ToLower : public std::unary_function<char,char> {
        char operator()(char a)
        {
                return tolower(a);
        };
};
struct ToUpper : public std::unary_function<char,char> {
        char operator()(char a)
        {
                return toupper(a);
        };
};

template<typename A, typename B>
std::pair<B,A> flip_pair(const std::pair<A,B> &p)
{
    return std::pair<B,A>(p.second, p.first);
}

template<typename A, typename B>
std::map<B,A> flip_map(const std::map<A,B> &src)
{
    std::map<B,A> dst;
    std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()), 
                   flip_pair<A,B>);
    return dst;
}

/* Auxiliary functions for checking input for validity. */

/* Function used to check that 'opt1' and 'opt2' are not specified
   at the same time. */
// void conflicting_options(const po::variables_map& vm,
//                          const char* opt1, const char* opt2)
// {
//     if (vm.count(opt1) && !vm[opt1].defaulted()
//         && vm.count(opt2) && !vm[opt2].defaulted())
//         throw std::logic_error(std::string("Conflicting options '")
//                           + opt1 + "' and '" + opt2 + "'.");
// }

/* Function used to check that of 'for_what' is specified, then
   'required_option' is specified too. */
// void option_dependency(const po::variables_map& vm,
//                         const char* for_what, const char* required_option)
// {
//     if (vm.count(for_what) && !vm[for_what].defaulted())
//         if (vm.count(required_option) == 0 || vm[required_option].defaulted())
//             throw std::logic_error(std::string("Option '") + for_what
//                               + "' requires option '" + required_option + "'.");
// }

/*
void tokenize(const std::string& str,
                      std::vector<std::string>& tokens,
                      const std::string& delimiters = ",")
{
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}
*/

std::string ss2Hishapeh(const std::string ss)
{
  int i, last_upbracket_position=0;
  std::stack<int> st;
  std::ostringstream hishapehOss;

  for (i=0; i<ss.length(); i++) 
  { 
    if (ss[i]=='(')
    {
      last_upbracket_position = i;
      st.push(i);
    }
    else if (ss[i]==')')
    {
      if (last_upbracket_position==st.top() )
      {
	st.pop();
	int double_hishapeh = (last_upbracket_position + 1 + i+1);
        if (double_hishapeh == double_hishapeh/2*2) // even
        {
	  hishapehOss << double_hishapeh/2;
	}
	else
	{
	  hishapehOss << double_hishapeh/2 << ".5";
	}
	hishapehOss << ",";
      }
    }
  }
  std::string hishapehStr = hishapehOss.str();
  return hishapehStr.substr(0, hishapehStr.length() - 1);
  //butt.erase( butt.size() - 1 );
  //return butt;
  //return hishapehOss.str();
}


// calculating topology of landscape regarding hishape
void generateMatchStr(const std::string& hishape1,const std::string& hishape2,std::set<std::string>* unionSet)
{  
    // delete the pair of square brackets []
    std::string strHishape1 = hishape1;  // deep copy, create a new string
    std::string strHishape2 = hishape2;
    strHishape1.erase(strHishape1.begin()+strHishape1.size()-1);
    strHishape1.erase(strHishape1.begin());
    strHishape2.erase(strHishape2.begin()+strHishape2.size()-1);
    strHishape2.erase(strHishape2.begin());

    std::vector<std::string> helixIndices1, helixIndices2;
    tokenize(strHishape1, helixIndices1); 
    tokenize(strHishape2, helixIndices2);
    int size1 = (signed)helixIndices1.size();
    int size2 = (signed)helixIndices2.size();
    int i = 0;
    
    // delete '(' and ')'
    for (i=size1-1; i>=0; i--)
    {
      if ("(" == helixIndices1[i] || ")" == helixIndices1[i])
          helixIndices1.erase (helixIndices1.begin()+i);
    }
    for (i=size2-1; i>=0; i--)
    {
      if ("(" == helixIndices2[i] || ")" == helixIndices2[i])
          helixIndices2.erase (helixIndices2.begin()+i);
    }
     
    // http://stackoverflow.com/questions/6955578/subtraction-and-intersection-of-two-vectors-of-pointers-in-c
    std::sort(helixIndices1.begin(), helixIndices1.end());
    std::sort(helixIndices2.begin(), helixIndices2.end());

    std::vector<std::string> unionVec; //Union of V1 and V2
    set_union(helixIndices1.begin(),helixIndices1.end(), helixIndices2.begin(),helixIndices2.end(), std::back_inserter(unionVec)); //myvec3: 1 3 10
    
    std::vector<std::string>::iterator it;
    for ( it=unionVec.begin() ; it != unionVec.end(); it++ )
    {
      (*unionSet).insert(*it);
    }
}

void generateMatchStr2(const std::string& hishapeh1,const std::string& hishapeh2,std::set<std::string>* unionSet)
{  
    // delete the pair of square brackets []
    std::string strHishapeh1 = hishapeh1;  // deep copy, create a new string
    std::string strHishapeh2 = hishapeh2;

    std::vector<std::string> helixIndices1, helixIndices2;
    tokenize(strHishapeh1, helixIndices1); 
    tokenize(strHishapeh2, helixIndices2);
    
    // http://stackoverflow.com/questions/6955578/subtraction-and-intersection-of-two-vectors-of-pointers-in-c
    std::sort(helixIndices1.begin(), helixIndices1.end());
    std::sort(helixIndices2.begin(), helixIndices2.end());

    std::vector<std::string> unionVec; //Union of V1 and V2
    set_union(helixIndices1.begin(),helixIndices1.end(), helixIndices2.begin(),helixIndices2.end(), std::back_inserter(unionVec)); //myvec3: 1 3 10
    
    std::vector<std::string>::iterator it;
    for ( it=unionVec.begin() ; it != unionVec.end(); it++ )
    {
      (*unionSet).insert(*it);
    }
}

void calculateHishapes(int hishape_type, std::vector<std::pair<const char *, unsigned> > &inputs, uint32_t kbest, const std::string &match_str, std::vector<std::string>* hishreps/*, std::vector<std::string>* hishapes*/, bool printDetails = true)
{
	//std::vector<std::string> hishreps;
	//std::vector<float> hishrepenergies;
        std::map<std::string, std::string> pp_hishape;
	
	   
        // ########################
	// ## calculate hishapes ##
	// ########################
	try {
	    std::vector<std::string> tokens;
	    if (hishape_type==4) {
		gapc::hishapeh_pp_cls obj;
		
		obj.init(inputs, kbest, match_str, false);
		
		// placeholder
		obj.cyk();
		gapc::hishapeh_pp_ret res = obj.run();
		//## gapc::add_event("end_computation");
		
		//obj.print_result(std::cout, res);
		//##gapc::add_event("end_result_pp");

                intrusive_ptr<Backtrace<String, unsigned int> >  bt   = obj.backtrack(obj.t_0_left_most);
		intrusive_ptr<Backtrace_List<String, unsigned int> > l =
		  boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
		assert(!bt || (bt && l));
		if (l) {
		    /* save the data in a vector 'tokens' */
		    for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
			std::stringstream ss;
			(*i)->print(ss);

			std::copy(std::istream_iterator<std::string>(ss),
				  std::istream_iterator<std::string>(),
				  std::back_inserter<std::vector<std::string> >(tokens) );
		    }	 
		}	  
	    } else if (hishape_type==3) {
		gapc::hishapehplus_pp_cls obj;
		obj.init(inputs, kbest, match_str, false);
		obj.cyk();
		gapc::hishapehplus_pp_ret res = obj.run();
		// obj.print_result(std::cout, res);   
                intrusive_ptr<Backtrace<String, unsigned int> >  bt   = obj.backtrack(obj.t_0_left_most);
		intrusive_ptr<Backtrace_List<String, unsigned int> > l =
		  boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
		assert(!bt || (bt && l));
		if (l) {
		    /* save the data in a vector 'tokens' */
		    for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
			std::stringstream ss;
			(*i)->print(ss);

			std::copy(std::istream_iterator<std::string>(ss),
				  std::istream_iterator<std::string>(),
				  std::back_inserter<std::vector<std::string> >(tokens) );
		    }	 
		}		  
	    } else if (hishape_type==2) {
		gapc::hishapem_pp_cls obj;
		obj.init(inputs, kbest, match_str, false);
		obj.cyk();
		gapc::hishapem_pp_ret res = obj.run();
		//TODO: what is function of print_result(...)? for nothing?!
		//obj.print_result(std::cout, res);  
                intrusive_ptr<Backtrace<String, unsigned int> >  bt   = obj.backtrack(obj.t_0_left_most);
		intrusive_ptr<Backtrace_List<String, unsigned int> > l =
		  boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
		assert(!bt || (bt && l));
		if (l) {
		    /* save the data in a vector 'tokens' */
		    for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
			std::stringstream ss;
			(*i)->print(ss);

			std::copy(std::istream_iterator<std::string>(ss),
				  std::istream_iterator<std::string>(),
				  std::back_inserter<std::vector<std::string> >(tokens) );
		    }	 
		}		 
	    } else if (hishape_type==1) {
		gapc::hishapeb_pp_cls obj;
		obj.init(inputs, kbest, match_str, false);
		obj.cyk();
		gapc::hishapeb_pp_ret res = obj.run();
		//obj.print_result(std::cout, res); 
		intrusive_ptr<Backtrace<String, unsigned int> >  bt   = obj.backtrack(obj.t_0_left_most);
		intrusive_ptr<Backtrace_List<String, unsigned int> > l =
		  boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
		assert(!bt || (bt && l));
		if (l) {
		    /* save the data in a vector 'tokens' */
		    for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
			std::stringstream ss;
			(*i)->print(ss);

			std::copy(std::istream_iterator<std::string>(ss),
				  std::istream_iterator<std::string>(),
				  std::back_inserter<std::vector<std::string> >(tokens) );
		    }	 
		}	 
	    } else {
		std::cout << "It does not exist hishape abstract type " << hishape_type << "." << std::endl;
	    }

	    
	    
	    // #########################
	    // ## processing commonly ##
	    unsigned int i;
	    // print sequence as well as its length
	    if (printDetails) {
		std::vector<std::pair<const char*, unsigned> >::iterator it;
		for ( it=inputs.begin() ; it < inputs.end(); it++ ) {
		    std::cout << "length = " << (*it).second << std::endl;
		    for (i = 0; i < (*it).second; i++)
		    {
			putchar(toupper((*it).first[i]));
		    }
		    std::cout << std::endl;
		}
	    }
	        
	    // calculate the longest hishape and mfe
	    unsigned int length_longest_hishape = 0, length_hishape = 0, length_longest_mfe = 0, length_mfe = 0;
	    for (i=0; 9*i < tokens.size(); i++) {
	      	tokens[9*i+2].erase (tokens[9*i+2].end()-1); 
			
		length_hishape = tokens[9*i+2].length();
		length_mfe = tokens[9*i+4].length();
		if (length_longest_hishape < length_hishape) {
		    length_longest_hishape = length_hishape;
		    length_longest_mfe = length_mfe;
		}
	    }
	    
	    // prepare format_pattern and output the results
	    std::string format_pattern = ( boost::format("%%%d.2f %%%ds") % (length_longest_mfe+5) % (length_longest_hishape+5) ).str();
	    i = 0;
	    for (i=0; 9*i < tokens.size(); i++) {  
														    // #### wrap hishape with "[...]" ####
		if (printDetails)  std::cout << tokens[9*i+7] << (boost::format(format_pattern) % ((boost::lexical_cast<int>(tokens[9*i+4]))/100.0f) % ("["+tokens[9*i+2]+"]")).str() << std::endl;
		(*hishreps).push_back(tokens[9*i+7]);
	    }

	    //======for (unsigned int i = 0; i<opts.repeats; ++i) obj.print_backtrack(std::cout, res);======                      
	} catch (std::exception &e) {
	    std::cerr << "Exception: " << e.what() << '\n';
	    std::exit(1);
	}
	

      #ifdef STATS
	obj.print_stats(std::cerr);
	// FIXME delete
	//Singleton<Hash::Set<Shape> >::ref().print_stats(std::cerr);
      #endif
      gapc::print_events(std::cerr);

      // Since static String::pool is in another translation unit
      // order of destruction is not defined -> hash needs to be destroyed before
      // FIXME delete
      //Singleton<Hash::Set<String> >::ref().purge();

}    


float hipathPairwise(const std::string &header, const std::string &seq, const std::string &ss1, const std::string &ss2, 
		    int hishape_type, std::vector<std::pair<const char *, unsigned> > &inputs, uint32_t kbest, bool printDetails = true) {
					      // calculate related hishape
					      std::set<std::string> related_hishape;
					      //$$ std::cout << ss2Hishapeh(ss1) << "|" << ss2Hishapeh(ss2) << std::endl;
					      //std::cout << "ss1=" << ss1 << std::endl;
					      //std::cout << "sSToHishape(ss1, 2)=" << sSToHishape(header+"\n"+seq+"\n"+ss2+"\n", hishape_type) << std::endl;
                                              unsigned found1 = ss1.find('(');
                                              unsigned found2 = ss2.find('(');
                                              std::string hishape1="[]", hishape2="[]";
					      //std::cout << found1 << "," << found2 << "," << std::string::npos << std::endl;
                                              if (found1<ss1.length())
                                              {
						  //std::cout << "1true" << std::endl;
                                                  hishape1=sSToHishape(header+"\n"+seq+"\n"+ss1+"\n", hishape_type);
					      }
					      if (found2<ss2.length())
					      {
						  //std::cout << "2true" << std::endl;
                                                  hishape2=sSToHishape(header+"\n"+seq+"\n"+ss2+"\n", hishape_type); 
                                              }
                                              generateMatchStr(hishape1, hishape2, &related_hishape);
					      //generateMatchStr2(ss2Hishapeh(ss1), ss2Hishapeh(ss2), &related_hishape);
					      //std::cout << pp_hishape[hishreps[i]] << "+" << pp_hishape[hishreps[j]] << "=";
					      std::set<std::string>::iterator it;
					      //std::cout << "myset contains:";
					      std::ostringstream match_oss;
					      for ( it=related_hishape.begin() ; it != related_hishape.end(); it++ )
					      {
						// update the match_str
						match_oss << *it << ",";
						//std::cout << " " << *it;
					      }
					      std::string match_str = match_oss.str();
					      //DEBUG	std::cout << match_str << std::endl;
					      
					      std::vector<std::string> related_hishreps;
					      // calculate related hishapes
					      calculateHishapes(hishape_type, inputs, kbest, match_str, &related_hishreps, printDetails);
					            
					      
					      //s = 0;
					      //t = 1;
					      
					      int s = -1;
					      int t = -1;
					      int i=0;
					      for (i=0; i<related_hishreps.size(); i++ )
					      {
						  if (related_hishreps[i] == ss1)  s=i;
						  if (related_hishreps[i] == ss2)  t=i;
					      }
					      //assert(s!=-1);
					      //assert(t!=-1);
					      if (s==-1) {
					          related_hishreps.push_back(ss1);
						  s = related_hishreps.size()-1;
					      }
					      if (t==-1) {
					          related_hishreps.push_back(ss2);
						  t = related_hishreps.size()-1;
					      }
					      //$$ printf("s=%d, t=%d\n", s, t);
					      std::string tmp_s = related_hishreps[s];
					      std::string tmp_t = related_hishreps[t];
					      std::vector<std::string>::iterator related_hishreps_begin = related_hishreps.begin();
					      if (s>t) {
                                                  related_hishreps.erase(related_hishreps_begin+s);
						  related_hishreps.erase(related_hishreps_begin+t);
					      } else {
						  related_hishreps.erase(related_hishreps_begin+t);
						  related_hishreps.erase(related_hishreps_begin+s);
					      }
					      related_hishreps.insert(related_hishreps_begin, tmp_t);
					      related_hishreps.insert(related_hishreps_begin, tmp_s);
                                              // with the swap above now s=0 and t=1
					      s = 0;
					      t = 1;
					      
					      int maxkeep, dist;
					      path_t *route;
					      maxkeep=10;
					      int k;
					      // update n_node, because n_node used in dijkstra_single_target(...)
					      n_node = related_hishreps.size();
					      
					      // given: n_node  
					      // return: anchor_no and anchors[]
					      int anchor_no = dijkstra_single_target(seq,related_hishreps,s,t,maxkeep);
					      //$$ printf("anchor_no=%d\n", anchor_no);
					      //$$ std::cout << "================================" << std::endl;
					      //$$ std::cout << seq << std::endl;
					      //$$ for (i = 0; i < n_node; ++i) {
                                              //$$     std::cout << related_hishreps[i] << std::endl;
					      //$$ }
					      					     
					      int k_max_saddle=0;
					      float en_max_saddle=-INFINITY, en_max_0=-INFINITY;
					      char struc_max_saddle[HISHREPLENGTH];
					      for (i = anchor_no-1; i >= 1; i--) {
						  route=get_path(const_cast<char*>(seq.c_str()), const_cast<char*>(related_hishreps[anchors[i]].c_str()), const_cast<char*>(related_hishreps[anchors[i-1]].c_str()), maxkeep);
						  dist=getBasePairDistance(const_cast<char*>(related_hishreps[anchors[i]].c_str()), const_cast<char*>(related_hishreps[anchors[i-1]].c_str()));
						  if (i==anchor_no-1) {
						      for (k=0;k<dist+1;k++){
							  if (printDetails)  printf("%2d: %8g  %s\n", k, route[k].en, route[k].s);
							  if (route[k].en > en_max_saddle)
							  {
							      k_max_saddle = k;
							      en_max_saddle = route[k].en;
							      strcpy(struc_max_saddle, route[k].s);
							  }
							  if (i==(anchor_no-1) && k==0) 
							  {
							      en_max_0 = route[k].en;
							  }
						      }
						  } else {
						      for (k=1;k<dist+1;k++){
							  if (printDetails)  printf("%2d: %8g  %s\n", k, route[k].en, route[k].s);
							  if (route[k].en > en_max_saddle)
							  {
							      k_max_saddle = k;
							      en_max_saddle = route[k].en;
							      strcpy(struc_max_saddle, route[k].s);
							  }
							  if (i==(anchor_no-1) && k==0) 
							  {
							      en_max_0 = route[k].en;
							  }
						      }
						  }
						  for (k=0;k<dist+1;k++){
						      free(route[k].s);
						  }
						  free(route);
					      }      
					      if (printDetails)  printf("Saddle structure: %s\n", struc_max_saddle);
					      //printf("%g %g kcal/mol\n", en_max_saddle, en_max_0);
					      if (printDetails)  printf("Saddle energy regarding start structure: %g kcal/mol\n", (en_max_saddle-en_max_0));
					      return en_max_saddle;
}
					      
					      

int main(int argc, char *argv[]) {
  

    // ########################
    // ## options            ##
    // ########################
        /* options declaration */
    //bool window_mode = false;
    //#ifdef WINDOW_MODE
        //    bool window_mode = true;
    //#endif

    #ifdef WINDOW_MODE
            unsigned int window_size;
            unsigned int window_increment;
    #endif
    unsigned int repeats = 1;
    typedef std::vector<std::pair<const char*, unsigned> > inputs_t;
    inputs_t inputs;
    //char *input = 0;                             // whether a sequence is given by user
    unsigned int hishape_type;
    uint32_t kbest;//, number;//, floor;  // , uint32_t kbest
    double T = 37.0;
    std::string par_filename = "";
    bool all = false;
    //HBT bool tree = false;
    bool readss = false;
    std::string ss_filename = "";

  
    try {
      
        // idea is the option is not compulsory and it will be given a relexed code  
        //po::options_description desc("Allowed options", 120, 60);
        po::options_description desc("Allowed options");
	
	desc.add_options()
	#ifdef WINDOW_MODE
		      ("windowsize,w)", po::value<unsigned int>(&window_size)->default_value(0), "Specify window size")
		      ("windowincrement,p", po::value<unsigned int>(&window_increment)->default_value(0), "Specify window position increment (use with -w)")
	#endif
	//("delta,d", po::value<unsigned int>(&delta)->default_value(100), "Set energy range (kcal/100mols)")            
	//("repeats,r", po::value<unsigned int>(&repeats)->default_value(1), "Sampling with specified times")
	("file,f", po::value< std::vector<std::string> >()/*->required()*/, "Read sequence(s) from specified file")
	//("sequence,s", po::value< std::vector<std::string> >(), "Read sequence from specified sequence")
	("type,t", po::value<unsigned int>(&hishape_type)->default_value(1), "Specify hishape abstract type")
	("kbest,k", po::value<unsigned int>(&kbest), "Choose the k-best related hishapes as anchor points")  //[only make sense in all-agsinst-all model]    ->default_value(20)
        ("Temperature,T", po::value<double>(&T)->implicit_value(37.0), "Specifiy the temperature for energy parameters in degree Celsius.") 
	("Parameter,P", po::value<std::string>(&par_filename), "Specifiy the path of energy parameter file.")
	("all,a", po::value<bool>(&all)->zero_tokens(), "Calculate folding pathways all-to-all")
	("readss,r", po::value<bool>(&readss)->zero_tokens(), "Use anchor secondary structures from specified file instead of generating them from sequence(s)")
	("ssfile,s", po::value<std::string>(&ss_filename), "Read secondary structures from specified file")
	//HBT ("tree,b", po::value<bool>(&tree)->zero_tokens(), "Generate barrier tree for hienergies")
	//("number,n", po::value<unsigned int>(&number)->default_value(20), "Specify minimal number of anchor hishapes required to trigger short path calculation")	
	//("floor,l", po::value<unsigned int>(&floor)->default_value(2), "Specify a hishape abstraction type that this program furthest reaches, in other words, the option sets a stop point "
        //"during the iterative calculating candidate hishape anchors from the most abstraction type to the least abstraction type")
	("help,h", "Produce help message")
	("version,v", "Show version");
	
        /*
         * it can add at most additional 2 records, that means,
         * it is still probable that the sum of the records >= 2
         *
         * ./main --help -f aaa dddd ccc -s bbb
         * will get the result
         * input files are: aaa
         * input sequences are: dddd ccc bbb
         * positional options are: dddd ccc bbb
         */
        po::positional_options_description p;
        p.add("file", -1);


        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(desc).positional(p).run(), vm);
        //po::store(parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        // option exclusions and implications
	conflicting_options(vm, "help", "version");
	conflicting_options(vm, "file", "sequence");
	//HBT option_dependency(vm, "tree", "all");
        option_dependency(vm, "ssfile", "readss");
	
	
        //std::cout << "vm.size()=" << vm.size() << std::endl;
        if (vm.count("help") || vm.size()==1) {
            std::cout << "Usage: HiPath -f INPUT-file [options]\n";  // (-s INPUT|)
	    std::cout << "Examples:" << std::endl;
	    std::cout << "  ./HiPath -f ../examples/riboswitches.faa" << std::endl;
	    std::cout << "  ./HiPath -f ../examples/switches_4.faa -k 40" << std::endl;
	    std::cout << "  ./HiPath -f ../examples/collosoma_slrna.fa -a -k 20" << std::endl;
	    std::cout << "  ./HiPath -f ../examples/xbix.fa -a -k 40 -r -s ../examples/xbix.hi" << std::endl;
            std::cout << desc; //<< '\n';
            //std::cout << "input files are: " << vm["file"].as< std::vector<std::string> >() << '\n';
            //std::cout << "input sequences are: " << vm["sequence"].as< std::vector<std::string> >() << '\n';
            //std::cout << "positional options are: " << vm["sequence"].as< std::vector<std::string> >() << '\n';  // if it is empty, return bad_any_cast error
            return 0;
        }

        if (vm.count("version")) {
            std::cout << "HiPath 2.0.12 (Feb. 09, 2012)\n";  // TODO: use the C++ buildin date function 
            return 0;
        }

        if (vm.count("type"))
	{
	    if (hishape_type < 1 || hishape_type > 4) 
	    {
		std::cout << "The hishape abstract type can only be specified between 1 and 4." << std::endl;  // For multiple hishape alignment, 
		return 1;
	    }
	}
	
	if (vm.count("Temperature"))
	{
	    if( T < -273.15 ){
	      fprintf(stderr, "Value of --Temp must be > -273.15\n");
	      exit (EXIT_FAILURE);
	    }
#ifdef RNALIB_H
            temperature = T;
#endif
	}

	if (vm.count("Parameter"))
	{
#ifdef RNALIB_H
            librna_read_param_file(/*par_filename*/par_filename.c_str());
#endif
	}
	else
	{
#ifdef RNALIB_H
	    librna_read_param_file(0);
#endif
	}


        if (vm.count("file"))
        {        
	    std::vector<std::string> files = vm["file"].as< std::vector<std::string> >();    
            if (!files.empty()) {    
	        // in option --all the input file should be fasta-format
	        if (files.size()==1 && all==1) 
		{
			std::string line;
			std::fstream fiInputFile;
			if (false == fiInputFile.is_open())
			{       
			      // get first input file
			      fiInputFile.open(files[0].c_str(), std::ios::in);
			      // checks if the file was opened
			      if ((false == fiInputFile) || (false == fiInputFile.is_open()))
			      {
				      // not opened, throws an exception
				      throw "ERROR: loading input file (\"" + files[0] + "\")!";
			      }
			      
			      std::string header = "";
			      std::string seq = "";
			      int line_no = 0;
			      while ( getline(fiInputFile, line) )
			      {

				      if ( (line_no%2 == 0) ) {
					  if (line[0] == '>')
					  {
					      // get header only
					      //header = fasta.substr(1);
					      header = line;
					  }
					  else
					  {
					      throw "ERROR: wrong faa-format in input file (\"" + files[0] + "\")!";
					  }
				      }
				      else 
				      {
					  if (line_no%2 == 1)
					      seq = line;
				      }
					      

				      if (header != "" && seq != "") {
					      std::cout << header << std::endl;
					      std::cout << seq << std::endl;
					     
					      if (!vm.count("kbest")) {
						  kbest=(uint32_t)round(124000*pow(seq.length(),(-1.5)));
					      }
					      // ####################################################################################
					      // ##################### calculate hipath all-to-all #########################  
					      inputs.clear();
					      inputs.push_back(std::make_pair(seq.c_str(), seq.length()));
					      std::vector<std::string> structures;
					      if (!readss) {
					          calculateHishapes(hishape_type, inputs, kbest, "", &structures);
					      } else {
						  std::string ssLine;
						  std::fstream ssFileStream;
						  if (false == ssFileStream.is_open())
						  {       
							ssFileStream.open(ss_filename.c_str(), std::ios::in);
							if ((false == ssFileStream) || (false == ssFileStream.is_open()))
							{
								// not opened, throws an exception
								throw "ERROR: loading input file (\"" + ss_filename + "\")!";
							}
							//int t=0;
							while ( getline(ssFileStream, ssLine) )
							{
							  structures.push_back(ssLine.substr(10,seq.length()));
							  //printf("%6d    ",t+1);
							  std::cout << ssLine << std::endl;
							  //t++;
							}
						  }
					      }
					      
					      int ss_number = structures.size();
					      
					      float hi_energies[ss_number][ss_number];
					      // print the Hienergy
					      int i, j;
					      for (i = 0; i < ss_number; ++i) {
						      for (j = 0; j < ss_number; ++j) {
							      if (i!=j)
							      {
								      //TODO: consider to delete 2D array distances, it isn't used any more
								      // NOTE: tatally wrong
								      hi_energies[i][j] = hipathPairwise("", const_cast<char*>(seq.c_str()), const_cast<char*>(structures.at(i).c_str()), const_cast<char*>(structures.at(j).c_str()), hishape_type, inputs, kbest, false);
							      }
							      else
							      {
								      hi_energies[i][j] = 0.0f;
							      }
						      }
					      }
					      

					      // relax the table
					      for (i = 0; i < ss_number; ++i) {
						      for (j = i+1; j < ss_number; ++j) {
							      if (hi_energies[i][j] < hi_energies[j][i])
								  hi_energies[j][i] = hi_energies[i][j];
							      else  // hi_energies[i][j] >= hi_energies[j][i]
								  hi_energies[i][j] = hi_energies[j][i];
						      }
					      }
					      
					      // output the result
					      //puts("----------hi_energies---------------------------\n");
					      for (i = 0; i < ss_number; ++i) {
						      std::cout << "S";
						      printf("%10d",i);
						      //for (j = 0; j < i; ++j) {
							//      printf("%10.4g", hi_energies[i][j]);
						      //}
						      for (j = 0; j < ss_number; ++j) {
							  if (j<i) {
							      printf("%10.4g",hi_energies[i][j]);
							  } else if (j==i) {
							      printf("%10.4g",0.0);
							  } else {
							      printf("%10.4g",hi_energies[j][i]);
							  }
						      }
						      printf("\n");
					      }
					      
					      //HBT if(tree){
						//HBT puts("Generating HiEnergy Tree\n");
					      //HBT }
					      // ####################################################################################
					      header = "";
					      seq = "";
				      }
				      line_no ++;
			      }
			}
			else
			{
				std::cout << "ERROR: unable to open input file (\"" + files[0] + "\")!\n";
			}
			
			if (true == fiInputFile.is_open())
			{
				fiInputFile.close();
			}
		}
		else if (files.size()==1)
		{
			std::string line;
			std::fstream fiInputFile;
			if (false == fiInputFile.is_open())
			{       
			      // get first input file
			      fiInputFile.open(files[0].c_str(), std::ios::in);
			      // checks if the file was opened
			      if ((false == fiInputFile) || (false == fiInputFile.is_open()))
			      {
				      // not opened, throws an exception
				      throw "ERROR: loading input file (\"" + files[0] + "\")!";
			      }
			      
			      std::string header = "";
			      std::string seq = "";
			      std::string ss1 = "";
			      std::string ss2 = "";
			      int line_no = 0;
			      while ( getline(fiInputFile, line) )
			      {

				      if ( (line_no%4 == 0) ) {
					  if (line[0] == '>')
					  {
					      // get header only
					      //header = fasta.substr(1);
					      header = line;
					  }
					  else
					  {
					      throw "ERROR: wrong faa-format in input file (\"" + files[0] + "\")!";
					  }
				      }
				      else 
				      {
					  if (line_no%4 == 1)
					      seq = line;
					  if (line_no%4 == 2)
					      ss1 = line;
					  if (line_no%4 == 3)
					      ss2 = line;
				      }


				      if (header != "" && seq != "" && ss1 != "" && ss2 != "") {
					      std::cout << header << std::endl;
					      std::cout << seq << std::endl;
					      std::cout << ss1 << std::endl;
					      std::cout << ss2 << std::endl;
					      
					      if (!vm.count("kbest")) {
						  kbest=(uint32_t)round(124000*pow(seq.length(),(-1.5))); 
					      }
					      // ####################################################################################
					      // ##################### calculate hipath pairwise #########################  
					      inputs.clear();
					      inputs.push_back(std::make_pair(seq.c_str(), seq.length()));
					      
					      //std::cout << header <<","<<seq<<","<<ss1<<","<<ss2<<","<<hishape_type<<","<<kbest<<std::endl; 
					      hipathPairwise(header, seq, ss1, ss2, hishape_type, inputs, kbest);

					      // ####################################################################################
					      header = "";
					      seq = "";
					      ss1 = "";
					      ss2 = "";
				      }
				      line_no ++;
			      }
			}
			else
			{
				std::cout << "ERROR: unable to open input file (\"" + files[0] + "\")!\n";
			}
			
			if (true == fiInputFile.is_open())
			{
				fiInputFile.close();
			}  
		}
		else
		{
                    throw "ERROR: loading input file (\"" + files[0] + "\")!";
		}
	    }
        }
		


/*
        if (vm.count("sequence"))
        {
	    std::vector<std::string> sequences = vm["sequence"].as< std::vector<std::string> >();
	        
	    unsigned int optind = 0;
	    for (; optind < sequences.size(); ++optind) {
		    std::string strBase=sequences[optind];
		    // convert to small letters
		    transform(strBase.begin(),strBase.end(),strBase.begin(),ToLower());
		    // t -> u
		    replace(strBase.begin(),strBase.end(),'t','u');
		    // delete '.'  and '-' from alignment files
		    remove(strBase.begin(),strBase.end(),'.');
		    remove(strBase.begin(),strBase.end(),'-');

		    input = new char[std::strlen(strBase.c_str())+1];
		    std::strcpy(input, strBase.c_str());
		    unsigned n = std::strlen(input);
		    inputs.push_back(std::make_pair(input, n));

	    }
	}
	*/

	

	
    } catch (std::exception &e) {
        std::cerr << "Exception: " << e.what() << '\n';
	    // http://stackoverflow.com/questions/9766568/unique-ptr-compile-error
	    //-std=c++0x  -Wall -Wextra -pedantic //g++ -std=c++0x -DHAVE_CONFIG_H -I.   -I./rtlib -I./librna   -g -O2 -MT HiPath-hipath.o -MD -MP -MF .deps/HiPath-hipath.Tpo -c -o HiPath-hipath.o `test -f 'hipath.cc' || echo './'`hipath.cc
	    //LDFLAGS=-ldl            //c++ -g -O2 -Wl,-R/usr/local/lib -Wl,-R/usr/local/lib -o .libs/HiPath HiPath-hipath.o HiPath-hited_util.o HiPath-hishapeh_mfe_pp.o HiPath-hishapehplus_mfe_pp.o HiPath-hishapem_mfe_pp.o HiPath-hishapeb_mfe_pp.o  -L/usr/local/lib -lboost_date_time-mt -lboost_program_options-mt ./libs/libRNA.a librna/.libs/librna.so rtlib/.libs/rtlib.so comp/.libs/comp.so -lm -ldl
// 	    std::cerr << "Caught a: " << exception_name << std::endl;
// 	    // print to stderr
// 	    backtrace_symbols_fd(last_frames, last_size, 2);
        std::exit(1);
    }
  
  


  
  

							

	

      
      /* delete all allocated memory */
      // it is not necessary because only std::string is used in this case
      //for (inputs_t::iterator i = inputs.begin(); i != inputs.end(); ++i)
      //	delete[] (*i).first;
	
	
      return 0;
}
