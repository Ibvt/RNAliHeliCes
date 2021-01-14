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
void conflicting_options(const po::variables_map& vm,
                         const char* opt1, const char* opt2)
{
    if (vm.count(opt1) && !vm[opt1].defaulted()
        && vm.count(opt2) && !vm[opt2].defaulted())
        throw std::logic_error(std::string("Conflicting options '")
                          + opt1 + "' and '" + opt2 + "'.");
}

/* Function used to check that of 'for_what' is specified, then
   'required_option' is specified too. */
void option_dependency(const po::variables_map& vm,
                        const char* for_what, const char* required_option)
{
    if (vm.count(for_what) && !vm[for_what].defaulted())
        if (vm.count(required_option) == 0 || vm[required_option].defaulted())
            throw std::logic_error(std::string("Option '") + for_what
                              + "' requires option '" + required_option + "'.");
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
    uint32_t kbest, number;//, floor;  // , uint32_t kbest
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
	("type,t", po::value<unsigned int>(&hishape_type)->default_value(2), "In all-agsinst-all model, it specifies hishape abstract type. "
	"In pairwise model, it specifies a hishape abstraction type that this program furthest reaches, in other words, the option sets a stop point "
	"during the iterative calculating candidate hishape anchors from the most abstraction type to the least abstraction type")
	("kbest,k", po::value<unsigned int>(&kbest)->default_value(3), "Choose the k-best classes [only make sense in all-agsinst-all model]") 
        ("Temperature,T", po::value<double>(&T)->implicit_value(37.0), "Specifiy the temperature for energy parameters in degree Celsius.") 
	("Parameter,P", po::value<std::string>(&par_filename), "Specifiy the path of energy parameter file.")
	("all,a", po::value<bool>(&all)->zero_tokens(), "Calculate folding pathways all-to-all")
	("readss,r", po::value<bool>(&readss)->zero_tokens(), "Use anchor secondary structures from specified file instead of generating them from sequence(s)")
	("ssfile,s", po::value<std::string>(&ss_filename), "Read secondary structures from specified file")
	//HBT ("tree,b", po::value<bool>(&tree)->zero_tokens(), "Generate barrier tree for hienergies")
	("number,n", po::value<unsigned int>(&number)->default_value(20), "Specify minimal number of anchor hishapes required to trigger short path calculation")	
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
        if (vm.count("help") || vm.size()==3) {
            std::cout << "Usage: BFSPath -f INPUT-file [options]\n";  // (-s INPUT|)
	    std::cout << "Examples:" << std::endl;
	    std::cout << "  ./BFSPath -f ../examples/xbix.fa -a -r -s ../examples/xbix.ss" << std::endl;
            std::cout << desc; //<< '\n';
            //std::cout << "input files are: " << vm["file"].as< std::vector<std::string> >() << '\n';
            //std::cout << "input sequences are: " << vm["sequence"].as< std::vector<std::string> >() << '\n';
            //std::cout << "positional options are: " << vm["sequence"].as< std::vector<std::string> >() << '\n';  // if it is empty, return bad_any_cast error
            return 0;
        }

        if (vm.count("version")) {
            std::cout << "BFSPath 2.0.12 (Feb. 09, 2012)\n";  // TODO: use the C++ buildin date function 
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
					      //std::cout << header << std::endl;
					      //std::cout << seq << std::endl;
					      
					      // ####################################################################################
					      // ##################### calculate BFSPath all-to-all #########################  
					      inputs.clear();
					      inputs.push_back(std::make_pair(seq.c_str(), seq.length()));
					      std::vector<std::string> structures;
// 					      if (!readss) {
// 					          calculateHishapes(hishape_type, inputs, kbest, "", &structures);
// 					      } else {
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
							
							while ( getline(ssFileStream, ssLine) )
							{
							  structures.push_back(ssLine);
							}
						  }
//					      }
					      
					      int ss_number = structures.size();
					      
					      float hi_energies[ss_number][ss_number];
					      // print the Hienergy
					      int i, j;
					      for (i = 0; i < ss_number; ++i) {
						      for (j = 0; j < ss_number; ++j) {
							      if (i!=j)
							      {
								      int k_max=0;
								      float en_max=-INFINITY;
								      char struc_max[HISHREPLENGTH];
								      int dist;
								      path_t *route, *r;
								      int k;

								      route=get_path( const_cast<char*>(seq.c_str()), const_cast<char*>(structures.at(i).c_str()), const_cast<char*>(structures.at(j).c_str()), 10);
								      //dist=getBasePairDistance(const_cast<char*>(structures[i].c_str()), const_cast<char*>(structures[j].c_str()));

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
								      
								      hi_energies[i][j] = en_max;
								      //printf("%g\n", en_max);
								      
							      }
							      else
							      {
								      hi_energies[i][j] = 0.0f;
							      }
						      }
					      }
					      

// 					      // relax the table
// 					      for (i = 0; i < ss_number; ++i) {
// 						      for (j = i+1; j < ss_number; ++j) {
// 							      if (hi_energies[i][j] < hi_energies[j][i])
// 								  hi_energies[j][i] = hi_energies[i][j];
// 							      else  // hi_energies[i][j] >= hi_energies[j][i]
// 								  hi_energies[i][j] = hi_energies[j][i];
// 						      }
// 					      }
					      
					      // check the intermediate result
					      //puts("----------hi_energies---------------------------\n");
					      for (i = 0; i < ss_number; ++i) {
						      for (j = 0; j < ss_number; ++j) {
							      printf("%0.3f\t", hi_energies[i][j]);
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
					      //## std::cout << seq << std::endl;
					      //## std::cout << ss1 << std::endl;
					      //## std::cout << ss2 << std::endl;
					      
					      // ####################################################################################
					      // ##################### calculate BFSPath pairwise #########################  
					      //inputs.clear();
					      //inputs.push_back(std::make_pair(seq.c_str(), seq.length()));
					      //BFSPathPairwise(header, seq, ss1, ss2, hishape_type, inputs, number);
					      
					      
					      int k_max=0;
					      float en_max=-INFINITY;
					      char struc_max[HISHREPLENGTH];
					      int dist;
					      path_t *route, *r;
					      int k;

					      route=get_path( const_cast<char*>(seq.c_str()), const_cast<char*>(ss1.c_str()), const_cast<char*>(ss2.c_str()), 10);
					      dist=getBasePairDistance(const_cast<char*>(ss1.c_str()), const_cast<char*>(ss2.c_str()));

					      for (k=0;k<dist+1;k++){
						  if (route[k].en > en_max)
						  {
						      k_max = k;
						      en_max = route[k].en;
						      strcpy(struc_max, route[k].s);
						  }
					      }

					      for (r=route; r->s; r++) {
						printf("%s %6.2f - %6.2f\n", r->s, energy_of_struct(seq.c_str(),r->s), r->en);
						free(r->s);
					      }
					      free(route);
					      
					      

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
		
	
    } catch (std::exception &e) {
        std::cerr << "Exception: " << e.what() << '\n';
        std::exit(1);
    }  
  
  
  


  
  

							

	

      
      /* delete all allocated memory */
      // it is not necessary because only std::string is used in this case
      //for (inputs_t::iterator i = inputs.begin(); i != inputs.end(); ++i)
      //	delete[] (*i).first;
	
	
      return 0;
}
