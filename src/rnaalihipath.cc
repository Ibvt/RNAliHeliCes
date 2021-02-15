#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>

#include "H/data_structures.h"
extern "C" {
#include "H/fold.h"
#include "H/fold_vars.h"
#include "H/utils.h"
#include "H/pair_mat.h"
#include "H/aln_util.h"
#include "H/alifold.h"
}
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <set>
#include <limits>
#include <cassert>

#include "probs_h_microstate.hh"
#include "probs_hplus_microstate.hh"
#include "probs_m_microstate.hh"
#include "probs_a_microstate.hh"
#include "eval_microstate.hh"

#include "rtlib/string.hh"
#include "rtlib/list.hh"
#include "rtlib/hash.hh"
#include "rtlib/asymptotics.hh"
#include "rtlib/generic_opts.hh"

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/replace.hpp>

//#include "hited_util.hh"
//#include "findpath.hh"
#include "Extensions/rnaoptions.hh"


namespace po = boost::program_options;



#define MAX_NUM_NAMES    500


//======for short path calculation======
#define GRAPHSIZE 2048
//#define INFINITY GRAPHSIZE*GRAPHSIZE
#define MAX(a, b) ((a > b) ? (a) : (b))
#define HISHREPLENGTH 1024
#define HISHREPSIZE 2048

//int e; /* The number of nonzero edges in the graph */
int n_node; /* The number of nodes in the graph */
int n_struc; /* The number of nodes in the graph */
float distances[GRAPHSIZE][GRAPHSIZE]; /* distances[i][j] is the distance between node i and j; or 0 if there is no direct connection */
float d[GRAPHSIZE]; /* d[i] is the length of the shortest path between the source (s) and node i */
int prev[GRAPHSIZE]; /* prev[i] is the node that comes right before i in the shortest path from the source to i*/
//char struc1[HISHREPLENGTH], struc2[HISHREPLENGTH];  //consensus_seq[HISHREPLENGTH], 
int anchors[100];

//TODO: AS and n_seq are global, not in accord with good programming convention ==> dangerous
int n_seq;
char     *AS[MAX_NUM_NAMES];          /* aligned sequences, a array of sequences */

/*
#################################
# PRIVATE VARIABLES of findpath.c           #
#################################
*/
const char  *seq=NULL;
short       *S=NULL, *S1=NULL;
int         BP_dist;
move_t      *path=NULL;
int         path_fwd; /* 1: struc1->struc2, else struc2 -> struc1 */


struct ToLower : public std::unary_function<char,char> {
        char operator()(char a)
        {
                return tolower(a);
        };
};

typedef char** stringArray;
stringArray MallocStringArray(size_t SizeOfOneString, size_t StringCount)
{
  char** t=(char**)malloc(StringCount*sizeof(char*));
  size_t i;
  for(i=0;i<StringCount;++i)
    t[i]=(char*)malloc(SizeOfOneString);
  return t;
}
//stringArray structures;  //, hishapes;


void FreeStringArray(stringArray StringArray, size_t StringCount)
{
  size_t i;
  for(i=0;i<StringCount;++i)
    free(StringArray[i]);
  free(StringArray);
} 

void replace_char (char* s, char find, char replace) {
  while (*s != 0) {
    if (*s == find)
      *s = replace;
    s++;
  }
}


//std::string sSToHishapeh(const std::string&);
int getBasePairDistance(char *ss1, char *ss2);
void printD();
void printPathIntoFile(FILE *file, int dest);
void printPath(int dest);
void dijkstra_all_targets(int s);
float calculateDistance(const std::string& sequence, const std::vector<std::string>& structures, int i, int j, int maxkeep);
int dijkstra_single_target(const std::string& sequence, const std::vector<std::string>& structures, int s, int t, int maxkeep);
template<typename A, typename B>
std::pair<B,A> flip_pair(const std::pair<A,B> &p);
template<typename A, typename B>
std::map<B,A> flip_map(const std::map<A,B> &src);
void conflicting_options(const po::variables_map& vm,
                         const char* opt1, const char* opt2);
void option_dependency(const po::variables_map& vm,
                        const char* for_what, const char* required_option);
void tokenize(const std::string& str,
                      std::vector<std::string>& tokens,
                      const std::string& delimiters = ",");
double energy_of_struc_on_as(const char* struc);
int energy_of_struc_on_as(const short* pt);
std::string sS2Hishapeh(const std::string& sS);
void generateMatchStr(const std::string& hishape1,const std::string& hishape2,std::set<std::string>* unionSet);
void calculateHishapes(int hishape_type, std::vector<std::pair<const char *, unsigned> > &inputs, uint32_t kbest, const std::string &match_str, std::vector<std::string>* hishreps/*, std::vector<std::string>* hishapes*/, bool printDetails);
float hipathPairwise(const std::string &header, const std::string &seq, const std::string &ss1, const std::string &ss2, 
		    int hishape_type, std::vector<std::pair<const char *, unsigned> > &inputs, uint32_t kbest, bool printDetails);

/*
#################################
# PUBLIC FUNCTION DECLARATIONS of findpath.h #
#################################
*/
/**
 *  \brief Find energy of a saddle point between 2 structures
 *  (serch only direct path)
 *
 *  \param seq RNA sequence
 *  \param struc1 A pointer to the character array where the first
 *         secondary structure in dot-bracket notation will be written to
 *  \param struc2 A pointer to the character array where the second
 *         secondary structure in dot-bracket notation will be written to
 *  \param max integer how many strutures are being kept during the search
 *  \returns the saddle energy in 10cal/mol
 */
int     find_saddle(const char *seq,
                    const char *struc1,
                    const char *struc2,
                    int max);


/**
 *  \brief Find refolding path between 2 structures
 *  (serch only direct path)
 *
 *  \param seq RNA sequence
 *  \param s1 A pointer to the character array where the first
 *         secondary structure in dot-bracket notation will be written to
 *  \param s2 A pointer to the character array where the second
 *         secondary structure in dot-bracket notation will be written to
 *  \param maxkeep integer how many strutures are being kept during the search
 *  \returns direct refolding path between two structures
 */
path_t* get_path( const char *seq,
                  const char *s1,
                  const char* s2,
                  int maxkeep);

/**
 *  \brief Free memory allocated by get_path() function
 *
 *  \param path pointer to memory to be freed
 */
void    free_path(path_t *path);

/*
#################################
# PRIVATE FUNCTION DECLARATIONS of findpath.c #
#################################
*/
int     *pair_table_to_loop_index (short *pt);
move_t  *copy_moves(move_t *mvs);
int     compare_ptable(const void *A, const void *B);
int     compare_energy(const void *A, const void *B);
int     compare_moves_when(const void *A, const void *B);
void    free_intermediate(intermediate_t *i);
int     find_path_once(const char *struc1, const char *struc2, int maxE, int maxl);
int     try_moves(intermediate_t c, int maxE, intermediate_t *next, int dist);


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
    char *input = 0;            // whether a sequence is given by user
    char *consensus_seq = 0;    // consensus sequence
    char *consensus_mis = 0;
    
    int /*n_seq,*/ length, i;
    //char     *AS[MAX_NUM_NAMES];          /* aligned sequences, a array of sequences */
    char     *names[MAX_NUM_NAMES];       /* sequence names */

    unsigned int hishape_type;
    uint32_t kbest;//, number;//, floor;  // , uint32_t kbest
    double T = 37.0;
    std::string par_filename = "";
    bool all = false;
    //HBT bool tree = false;
    std::string ss_filename = "";
    std::string startStruc="", targetStruc="";
    //std::vector<std::string> structures;
  
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
	("file,f", po::value< std::vector<std::string> >()/*->required()*/, "Read sequence alignment from specified file")
	//SEQ ("sequence,s", po::value< std::vector<std::string> >(), "Read sequence from specified sequence")
	("type,t", po::value<unsigned int>(&hishape_type)->default_value(1), "Specify hishape abstract type")
	("kbest,k", po::value<unsigned int>(&kbest), "Choose the k-best related hishapes as anchor points")  //[only make sense in all-agsinst-all model]    ->default_value(20)
        ("Temperature,T", po::value<double>(&T)->implicit_value(37.0), "Specifiy the temperature for energy parameters in degree Celsius.") 
	("Parameter,P", po::value<std::string>(&par_filename), "Specifiy the path of energy parameter file.")
	("ssfile,F", po::value<std::string>(&ss_filename), "Read start and target secondary structures from specified file")
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
	//SEQ conflicting_options(vm, "file", "sequence");
	option_dependency(vm, "file", "ssfile");
	//SEQ option_dependency(vm, "sequence", "ssfile");
	
        //std::cout << "vm.size()=" << vm.size() << std::endl;
        if (vm.count("help") || vm.size()==1) {
            std::cout << "Usage: RNAliHiPath -f INPUT-file [options]\n";  // (-s INPUT|)
	    std::cout << "Examples:" << std::endl;
	    //std::cout << "  ./RNAliHiPath -f ../examples/riboswitches.faa" << std::endl;
	    //std::cout << "  ./RNAliHiPath -f ../examples/switches_4.faa -k 10" << std::endl;
	    //std::cout << "  ./RNAliHiPath -f ../examples/test4.aln -F ../examples/test4.ss -P ./librna/vienna/rna_turner1999.par -t 2" << std::endl;
	    std::cout << "  ./RNAliHiPath -f ../examples/test2.aln -k 10 -F ../examples/test2.ss" << std::endl;
	    //SEQ std::cout << "  ./RNAliHiPath -s GGGUAACCACUAAAAUCCCGAAAGGGUGGGCU_GUGGGGACCCUCC#GGGUAACCACUAAAAUCCCGAAAGGGUGGGCU_GUGGUGACCUUCC#GGGUAACCACUAAAAUCCCGAAAGGGUGGGCU_GUGGUGACCUUCC#GGGUAACCACUAAAAUCCCGAAAGGGUGGGCUAGUGGCGACCCUCC#GGGUAACCACUAAAAUCCCGAAAGGGUGGGCU_GUGGUGACCUUCC# -k 10 -F ../examples/test2.ss" << std::endl;
            std::cout << desc; //<< '\n';
            //std::cout << "input files are: " << vm["file"].as< std::vector<std::string> >() << '\n';
            //std::cout << "input sequences are: " << vm["sequence"].as< std::vector<std::string> >() << '\n';
            //std::cout << "positional options are: " << vm["sequence"].as< std::vector<std::string> >() << '\n';  // if it is empty, return bad_any_cast error
            return 0;
        }

        if (vm.count("version")) {
            std::cout << "RNAliHiPath 1.2.2 (Feb. 14, 2021)\n";  // TODO: use the C++ buildin date function 
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
	
	if (vm.count("ssfile"))
	{
            if (ss_filename!="") {   
		std::string ssLine;
		std::fstream ssFileStream;
		if (false == ssFileStream.is_open())
		{       
		      ssFileStream.open(ss_filename.c_str(), std::ios::in);
		      //if ((false == ssFileStream) || (false == ssFileStream.is_open()))
		      if ( (ssFileStream.rdstate() & std::ifstream::failbit ) != 0 )
		      {
			      throw "ERROR: loading input file (\"" + ss_filename + "\")!";
		      }
		      std::string line;
		      int line_no = 0;
		      while ( getline(ssFileStream, line) )
		      {
			  if (line_no%2 == 0) 
			  {
			      startStruc = line;
			  }
			  else if(line_no%2 == 1) 
			  {
			      targetStruc = line;
			  }
			  line_no ++;
		      }
		}
	    }
	}
	
        //std::cout << vm.count("file") << "," << vm.count("ssfile") << std::endl;
        if (vm.count("file"))
        {         
	    std::vector<std::string> files = vm["file"].as< std::vector<std::string> >();   
            if (!files.empty()) {    
		if (files.size()==1)
		{
		    FILE     *clust_file = fopen ( files[0].c_str(), "r");  //stdin;
		        
		    n_seq = read_clustal(clust_file, AS, names);

		    //if (clust_file != stdin) 
		    fclose(clust_file);
		    if (n_seq==0)
		      puts("no sequences found");
		    length = (int) strlen(AS[0]);  // read structures
		    
		    input = (char*) malloc ((length+1)*n_seq+1);
		    if (input==NULL) exit (1);
		      
		    strcpy(input, AS[0]); /* copy name into the new var */
		    strcat(input, "#"); /* add the extension */
		    printf("%s    #%s\n", AS[0], names[0]);
		    for (i=1; i<n_seq; i++) {
			strcat(input,AS[i]);
			strcat(input, "#");
			printf("%s    #%s\n", AS[i], names[i]);
		    }
		    strcat(input, "\0");
		    replace_char(input, '-', '_');
		    inputs.push_back(std::make_pair(input, strlen(input))); 
		    
// 		    //#### invoke the eval_microstate to calculate the secondary structure energy based on alignment of sequences #### 
// 		    // allow Lonely Basepairs
//                     gapc::Opts::getOpts()->allowLonelyBasepairs = 1;
//                     gapc::Opts::getOpts()->alifold_minscore_basepair = -2000;
// 		    //const char* benchmark_struc = "(((((((..((((.........)))).(((((.......))))).....(((((.......)))))))))))).";  //".........((((.........)))).(((((.......)))))................";
// 		    //const char* benchmark_struc = ".(((....))).....(((....)))....................";  //".(((....))).....(((....))).........((......)).";  -3.600000, 0.000000, (null), .(((....))).....(((....)))....................
//                     const char* benchmark_struc = "((((.(((((......(((....))).......))))).))))...";  
// 		    //printf(" (%6.2f = %6.2f + %6.2f) \n", min_en, real_en, min_en-real_en );
// 		    
// 		    //std::cout << energy_of_struc_on_as(example_struc.c_str()) << std::endl;
// 		    //std::cout << energy_of_struc_on_as(make_pair_table(example_struc.c_str())) << std::endl;
// 		    //std::cout << "gapc::Opts::getOpts()->allowLonelyBasepairs=" << gapc::Opts::getOpts()->allowLonelyBasepairs << std::endl;
// 
// 		    //#### benchmark eval_microstate vs alifold ####
// 		    Pairs::getGivenPairs()->setStructure(std::make_pair(benchmark_struc, strlen(benchmark_struc)));
// 		    gapc::eval_cls obj;
// 		    obj.init(inputs);
// 		    obj.cyk();
// 		    gapc::eval_ret res = obj.run();
// 		    std::string out = boost::lexical_cast<std::string>(res);
// 		    std::cout << out << std::endl;
// 
// 		    float *ens  = (float *)space(2*sizeof(float));
// 		    energy_of_alistruct((const char **)AS, benchmark_struc, n_seq, ens);
//                     printf("%f, %f, %s, %s\n", ens[0], ens[1], consensus_seq, benchmark_struc);
		    
		    // TODO: 1. improvement performance 
		    //       2, both leak checks 
		    //       3. ENERGY is total DIFFERENT! check the total energy calculated by substraction or addition
		    //       4, transform the .... to xxxx ==> it is not necessary!
		        
		}
		else
		{
                    throw "ERROR: loading input file (\"" + files[0] + "\")!";
		}
	    }
        }

//SEQ 
// 	if (vm.count("sequence"))
//         {
// 	    std::vector<std::string> sequences = vm["sequence"].as< std::vector<std::string> >();
// 	        
// 	    unsigned int optind = 0;
// 	    for (; optind < sequences.size(); ++optind) {
// 		std::string strBase=sequences[optind];
// 		// convert to small letters
// 		transform(strBase.begin(),strBase.end(),strBase.begin(),ToLower());
// 		// t -> u
// 		replace(strBase.begin(),strBase.end(),'t','u');
// 		// delete '.'  and '-' from alignment files
// 		//remove(strBase.begin(),strBase.end(),'.');
// 		//remove(strBase.begin(),strBase.end(),'-');
// 		replace(strBase.begin(),strBase.end(),'-','_');
// 
// 		//input = new char[std::strlen(strBase.c_str())+1];
// 		input = (char*) malloc (std::strlen(strBase.c_str())+1);
// 		std::strcpy(input, strBase.c_str());
// 		unsigned n = std::strlen(input);
// 		inputs.push_back(std::make_pair(input, n));
// 
// 	    }
// 	}
	

	//consensus_mis = consens_mis((const char **) AS);
	//printf("%s    #most informative sequence\n", consensus_mis);
	consensus_seq = consensus((const char **) AS);
	printf("%s    #consensus\n", consensus_seq);
	replace_char(consensus_seq, '-', '_');  // preparing the consensus sequence for further usage
	std::string consensus_seq_str(consensus_seq);
	
	//std::cout << consensus_seq_str << "," << startStruc << "," << targetStruc << "," << hishape_type << "," << kbest << std::endl;
	hipathPairwise("", consensus_seq_str, startStruc, targetStruc, hishape_type, inputs, kbest, true);
    } catch (std::exception &e) {
        std::cerr << "Exception1: " << e.what() << '\n';
        std::exit(1);
    }
  
  


    //#### free resources ####      
    /* delete all allocated memory */
    //for (inputs_t::iterator i = inputs.begin(); i != inputs.end(); ++i)
    //    delete[] (*i).first;
    if (input != NULL)
        free(input);
    free(consensus_seq);
    free(consensus_mis);
     
    //KEY
    //if (cstruc!=NULL) free(cstruc);
    //if (cstruc2!=NULL) free(cstruc2);
    //free(base_pair);
    (void) fflush(stdout);

    //FreeStringArray(structures,HISHREPSIZE);
    //FreeStringArray(hishapes,HISHREPSIZE);
    //free(structure);
    //free(structure2);
    for (i=0; AS[i]; i++) {
      free(AS[i]); free(names[i]);
      //printf ("i=%d\n", i);
    }
    
    return(EXIT_SUCCESS);
}




//############################ private methods #########################

int getBasePairDistance(char *ss1, char *ss2) {
    short *pairTable1, *pairTable2;
    move_t *move_tList;
    int i, length, distance=0;
    pairTable1 = make_pair_table(ss1);
    pairTable2 = make_pair_table(ss2);
    length = (int) strlen(ss1);
    move_tList = (move_t *) space(sizeof(move_t)*length); 

    for ( i=1; i<=length; i++ ) {
        if (pairTable1[i] != pairTable2[i]) {
            if ( i<pairTable1[i] ) {
	        move_tList[distance].i = -i;
	        move_tList[distance].j = -pairTable1[i];
	        move_tList[distance++].when = 0;
            }
            if ( i<pairTable2[i] ) {
	        move_tList[distance].i = i;
	        move_tList[distance].j = pairTable2[i];
	        move_tList[distance++].when = 0;
            }
       }
  }
  free(pairTable1);
  free(pairTable2);
  free(move_tList);
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
	free(r->s);
      }
      free(route);
      //free_path(route);
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


void tokenize(const std::string& str,
                      std::vector<std::string>& tokens,
                      const std::string& delimiters)
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
 
    
double energy_of_struc_on_as(const char* struc){
  double energy=0.0;
  std::vector<std::string> tokens;
  
  try {
    Pairs::getGivenPairs()->setStructure(std::make_pair(struc, strlen(struc)));
  } catch (std::exception &e) {
    std::cerr << "Exception2: " << e.what() << '\n';
    std::exit(1);
  }
  gapc::eval_cls obj;

  try {
    obj.init(gapc::Opts::getOpts()->inputs);  //*gapc::Opts::getOpts());
  } catch (std::exception &e) {
    std::cerr << "Exception3: " << e.what() << '\n';
    std::exit(1);
  }

  // actual performance gains like 20%
  // see also http://www.ddj.com/cpp/184401305
 
  // workaround stupid Sun CC std::cout to fd0 after sync_with_stdio 
  // with -m64 and stlport4 bug:
  // http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=6804239
#if defined(__SUNPRO_CC) && __SUNPRO_CC <= 0x5100
  #warning Enable sync_with_stdio because of Sun CC 12 Compiler Bug
#else
  std::ios_base::sync_with_stdio(false);
#endif
  std::cin.tie(0);

#ifdef WINDOW_MODE
  unsigned n = obj.t_0_seq.size();
  for (unsigned int i = 0; ; i+=gapc::Opts::getOpts()->window_increment) {
    unsigned int right = std::min(n, i+gapc::Opts::getOpts()->window_size);
    gapc::return_type res = obj.run();
    std::cout << "Answer ("
      << i << ", " << right << ") :\n";
    obj.print_result(std::cout, res);
    for (unsigned int j = 0; j<gapc::Opts::getOpts()->repeats; ++j)
      obj.print_backtrack(std::cout, res);
    if (i+gapc::Opts::getOpts()->window_size >= n)
      break;
    obj.window_increment();
  }
#else
  gapc::add_event("start");

  obj.cyk();
  gapc::eval_ret res = obj.run();

  gapc::add_event("end_computation");

  //std::cout << "Answer: \n";
  //obj.print_result(std::cout, res);
  std::string out = boost::lexical_cast<std::string>(res);
  std::cout << out << std::endl;
  tokenize(out, tokens, " ");
  try {
      energy = boost::lexical_cast<double>(tokens[2]);
  } catch (std::exception &e) {
      std::cout << "out=" << out << std::endl;
      std::cout << "tokens[2]=" << tokens[2] << std::endl;
      std::cerr << "Exception13: " << e.what() << '\n';
      std::exit(1);
  }   
  
  gapc::add_event("end_result_pp");

#ifdef TRACE
  std::cerr << "start backtrack\n";
#endif
  for (unsigned int i = 0; i<gapc::Opts::getOpts()->repeats; ++i)
    obj.print_backtrack(std::cout, res);
  obj.print_subopt(std::cout, gapc::Opts::getOpts()->delta);

  gapc::add_event("end");
#endif

#ifdef STATS
  obj.print_stats(std::cerr);
#endif

  gapc::print_events(std::cerr);

  //return 0;
  return energy/100.0;
}


int energy_of_struc_on_as(const short* pt){
  double energy=0;
  std::vector<std::string> tokens;
  
  try {
    Pairs::getGivenPairs()->setStructure(pt);
  } catch (std::exception &e) {
    std::cerr << "Exception4: " << e.what() << '\n';
    std::exit(1);
  }
  gapc::eval_cls obj;

  try {
    obj.init(gapc::Opts::getOpts()->inputs);  //*gapc::Opts::getOpts());
  } catch (std::exception &e) {
    std::cerr << "Exception5: " << e.what() << '\n';
    std::exit(1);
  }

  // actual performance gains like 20%
  // see also http://www.ddj.com/cpp/184401305
 
  // workaround stupid Sun CC std::cout to fd0 after sync_with_stdio 
  // with -m64 and stlport4 bug:
  // http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=6804239
#if defined(__SUNPRO_CC) && __SUNPRO_CC <= 0x5100
  #warning Enable sync_with_stdio because of Sun CC 12 Compiler Bug
#else
  std::ios_base::sync_with_stdio(false);
#endif
  std::cin.tie(0);

#ifdef WINDOW_MODE
  unsigned n = obj.t_0_seq.size();
  for (unsigned int i = 0; ; i+=gapc::Opts::getOpts()->window_increment) {
    unsigned int right = std::min(n, i+gapc::Opts::getOpts()->window_size);
    gapc::return_type res = obj.run();
    std::cout << "Answer ("
      << i << ", " << right << ") :\n";
    obj.print_result(std::cout, res);
    for (unsigned int j = 0; j<gapc::Opts::getOpts()->repeats; ++j)
      obj.print_backtrack(std::cout, res);
    if (i+gapc::Opts::getOpts()->window_size >= n)
      break;
    obj.window_increment();
  }
#else
  gapc::add_event("start");

  obj.cyk();
  gapc::eval_ret res = obj.run();

  gapc::add_event("end_computation");

  //std::cout << "Answer: \n";
  //obj.print_result(std::cout, res);
  std::string out = "";
  try {
      out = boost::lexical_cast<std::string>(res);
      std::cout << out << std::endl;
  } catch (std::exception &e) {
      std::cerr << "Exception11: " << e.what() << '\n';
      std::exit(1);
  }   
  tokenize(out, tokens, " ");
//   if (tokens.size()!=13) {
//     std::cout << "tokens.size()=" << tokens.size() << std::endl;
//     std::cout << out << std::endl;
//   }

  try {
      energy = boost::lexical_cast<double>(tokens[2]);
  } catch (std::exception &e) {
      std::cout << "out=" << out << std::endl;
      std::cout << "tokens[2]=" << tokens[2] << std::endl;
      std::cerr << "Exception12: " << e.what() << '\n';
      std::exit(1);
  }   

  gapc::add_event("end_result_pp");

#ifdef TRACE
  std::cerr << "start backtrack\n";
#endif
  for (unsigned int i = 0; i<gapc::Opts::getOpts()->repeats; ++i)
    obj.print_backtrack(std::cout, res);
  obj.print_subopt(std::cout, gapc::Opts::getOpts()->delta);

  gapc::add_event("end");
#endif

#ifdef STATS
  obj.print_stats(std::cerr);
#endif

  gapc::print_events(std::cerr);

  //return 0;
  return (int)energy;
}
    

std::string sS2Hishapeh(const std::string& sS)
{
  int i, last_upbracket_position=0;
  std::stack<int> st;
  std::ostringstream hishapehOss;

  for (i=0; i<sS.length(); i++) 
  { 
    if (sS[i]=='(')
    {
      last_upbracket_position = i;
      st.push(i);
    }
    else if (sS[i]==')')
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
    //strHishape1.erase(strHishape1.begin()+strHishape1.size()-1);
    //strHishape1.erase(strHishape1.begin());
    //strHishape2.erase(strHishape2.begin()+strHishape2.size()-1);
    //strHishape2.erase(strHishape2.begin());

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

void calculateHishapes(int hishape_type, std::vector<std::pair<const char *, unsigned> > &inputs, uint32_t kbest, const std::string &match_str, std::vector<std::string>* hishreps/*, std::vector<std::string>* hishapes*/, bool printDetails = true)
{
        bool exact = false;
	int thresh = 1000;
	//std::vector<std::string> hishreps;
	//std::vector<float> hishrepenergies;
        std::map<std::string, std::string> pp_hishape;
	
	// ########################
	// ## calculate hishapes ##
	// ########################
  
	try {
	    std::vector<std::string> tokens;

	    if (hishape_type==4) {
		
		gapc::hishapeh_pfx_cls obj;
		//TODO: implement pfx with obj.init(inputs, kbest, match_str, exact, thresh*100, locmin, true);
		obj.init(inputs, kbest, match_str, exact, thresh*100);
  // actual performance gains like 20%
  // see also http://www.ddj.com/cpp/184401305

  // workaround stupid Sun CC std::cout to fd0 after sync_with_stdio 
  // with -m64 and stlport4 bug:
  // http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=6804239
#if defined(__SUNPRO_CC) && __SUNPRO_CC <= 0x5100
  #warning Enable sync_with_stdio because of Sun CC 12 Compiler Bug
#else
  std::ios_base::sync_with_stdio(false);
#endif
  std::cin.tie(0);
		obj.cyk();
		gapc::hishapeh_pfx_ret res = obj.run();
		obj.print_result(std::cout, res);
		
		//======for (unsigned int i = 0; i<gapc::Opts::getOpts()->repeats; ++i) obj.print_backtrack(std::cout, res);======
		intrusive_ptr<Backtrace<std::pair<String, double> , unsigned int> >  bt    = obj.backtrack(obj.t_0_left_most);
		intrusive_ptr<Backtrace_List<std::pair<String, double> , unsigned int> > l =
		    boost::dynamic_pointer_cast<Backtrace_List<std::pair<String, double> , unsigned int> > (bt);
		assert(!bt || (bt && l));
		//if (l) {
		for (Backtrace_List<std::pair<String, double> , unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
		    std::stringstream ss;
		    (*i)->print(ss);
		    std::copy(std::istream_iterator<std::string>(ss),
			      std::istream_iterator<std::string>(),
			      std::back_inserter<std::vector<std::string> >(tokens) );
		}
	    } 
	    else if (hishape_type==3) {
		gapc::hishapehplus_pfx_cls obj;
		obj.init(inputs, kbest, match_str, exact, thresh*100);
  // actual performance gains like 20%
  // see also http://www.ddj.com/cpp/184401305

  // workaround stupid Sun CC std::cout to fd0 after sync_with_stdio 
  // with -m64 and stlport4 bug:
  // http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=6804239
#if defined(__SUNPRO_CC) && __SUNPRO_CC <= 0x5100
  #warning Enable sync_with_stdio because of Sun CC 12 Compiler Bug
#else
  std::ios_base::sync_with_stdio(false);
#endif
  std::cin.tie(0);
		obj.cyk();
		gapc::hishapehplus_pfx_ret res = obj.run();
		obj.print_result(std::cout, res);   
		//======for (unsigned int i = 0; i<gapc::Opts::getOpts()->repeats; ++i) obj.print_backtrack(std::cout, res);======
		intrusive_ptr<Backtrace<std::pair<String, double> , unsigned int> >  bt    = obj.backtrack(obj.t_0_left_most);
		intrusive_ptr<Backtrace_List<std::pair<String, double> , unsigned int> > l =
		    boost::dynamic_pointer_cast<Backtrace_List<std::pair<String, double> , unsigned int> > (bt);
		assert(!bt || (bt && l));
		//if (l) {
		for (Backtrace_List<std::pair<String, double> , unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
		    std::stringstream ss;
		    (*i)->print(ss);
		    std::copy(std::istream_iterator<std::string>(ss),
			      std::istream_iterator<std::string>(),
			      std::back_inserter<std::vector<std::string> >(tokens) );
		}
	    } else if (hishape_type==2) {
		gapc::hishapem_pfx_cls obj;
		obj.init(inputs, kbest, match_str, exact, thresh*100);
  // actual performance gains like 20%
  // see also http://www.ddj.com/cpp/184401305

  // workaround stupid Sun CC std::cout to fd0 after sync_with_stdio 
  // with -m64 and stlport4 bug:
  // http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=6804239
#if defined(__SUNPRO_CC) && __SUNPRO_CC <= 0x5100
  #warning Enable sync_with_stdio because of Sun CC 12 Compiler Bug
#else
  std::ios_base::sync_with_stdio(false);
#endif
  std::cin.tie(0);
		obj.cyk();
		gapc::hishapem_pfx_ret res = obj.run();
		obj.print_result(std::cout, res);  
		//======for (unsigned int i = 0; i<gapc::Opts::getOpts()->repeats; ++i) obj.print_backtrack(std::cout, res);======
		intrusive_ptr<Backtrace<std::pair<String, double> , unsigned int> >  bt    = obj.backtrack(obj.t_0_left_most);
		intrusive_ptr<Backtrace_List<std::pair<String, double> , unsigned int> > l =
		    boost::dynamic_pointer_cast<Backtrace_List<std::pair<String, double> , unsigned int> > (bt);
		assert(!bt || (bt && l));
		//if (l) {
		for (Backtrace_List<std::pair<String, double> , unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
		    std::stringstream ss;
		    (*i)->print(ss);
		    std::copy(std::istream_iterator<std::string>(ss),
			      std::istream_iterator<std::string>(),
			      std::back_inserter<std::vector<std::string> >(tokens) );
		}
	    } else if (hishape_type==1) {
		gapc::hishapeb_pfx_cls obj;
		obj.init(inputs, kbest, match_str, exact, thresh*100);
  // actual performance gains like 20%
  // see also http://www.ddj.com/cpp/184401305

  // workaround stupid Sun CC std::cout to fd0 after sync_with_stdio 
  // with -m64 and stlport4 bug:
  // http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=6804239
#if defined(__SUNPRO_CC) && __SUNPRO_CC <= 0x5100
  #warning Enable sync_with_stdio because of Sun CC 12 Compiler Bug
#else
  std::ios_base::sync_with_stdio(false);
#endif
  std::cin.tie(0);
		obj.cyk();
		gapc::hishapeb_pfx_ret res = obj.run();
		obj.print_result(std::cout, res); 
		//======for (unsigned int i = 0; i<gapc::Opts::getOpts()->repeats; ++i) obj.print_backtrack(std::cout, res);======
		intrusive_ptr<Backtrace<std::pair<String, double> , unsigned int> >  bt    = obj.backtrack(obj.t_0_left_most);
		intrusive_ptr<Backtrace_List<std::pair<String, double> , unsigned int> > l =
		    boost::dynamic_pointer_cast<Backtrace_List<std::pair<String, double> , unsigned int> > (bt);
		assert(!bt || (bt && l));
		//if (l) {
		for (Backtrace_List<std::pair<String, double> , unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
		    std::stringstream ss;
		    (*i)->print(ss);
		    std::copy(std::istream_iterator<std::string>(ss),
			      std::istream_iterator<std::string>(),
			      std::back_inserter<std::vector<std::string> >(tokens) );
		}
	    } else {
		std::cout << "It does not exist hishape abstract type " << hishape_type << "." << std::endl;
	    }

	    
	    
	    // #########################
	    // ## processing commonly ##
	    unsigned int i;
	    // print sequence as well as its length
// 	    if (printDetails) {
// 		std::vector<std::pair<const char*, unsigned> >::iterator it;
// 		for ( it=inputs.begin() ; it < inputs.end(); it++ ) {
// 		    //std::cout << "length = " << (*it).second << std::endl;
// 		    for (i = 0; i < (*it).second; i++)
// 		    {
// 			putchar(toupper((*it).first[i]));
// 		    }
// 		    std::cout << std::endl;
// 		}
// 	    }
	    

	    // calculate the longest hishape and mfe
	    unsigned int length_longest_hishape = 0, length_hishape = 0, length_longest_mfe = 0, length_mfe = 0;
	    for (i=0; 25*i < tokens.size(); i++) {
		if (tokens[25*i+2]!="_")
		    tokens[25*i+2].erase (tokens[25*i+2].end()-1); 
		
		length_hishape = tokens[25*i+2].length();
		length_mfe = tokens[25*i+6].length();
		if (length_longest_hishape < length_hishape) {
		    length_longest_hishape = length_hishape;
		}
		if (length_longest_mfe < length_mfe) {
		    length_longest_mfe = length_mfe;
		}
	    }

	    // prepare format_pattern and output the results
	    std::string format_pattern = ( boost::format("%%%d.2f %%%ds") % (length_longest_mfe+5) % (length_longest_hishape+5) ).str();
	    i = 0;
	    for (i=0; 25*i < tokens.size(); i++) {  
														    // #### wrap hishape with "[...]" ####
		if (printDetails) std::cout << tokens[25*i+20] << (boost::format(format_pattern) % ((boost::lexical_cast<double>(tokens[25*i+6]))/100.0f) % ("["+tokens[25*i+2]+"]")).str() << std::endl;
		(*hishreps).push_back(tokens[25*i+20]);
	    }     
	} catch (std::exception &e) {
	    std::cerr << "Exception6: " << e.what() << '\n';
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
		    int hishape_type, std::vector<std::pair<const char *, unsigned> > &inputs, uint32_t kbest, bool printDetails = true) {  // int n_seq, char *AS[MAX_NUM_NAMES]
    // calculate related hishape
    std::set<std::string> related_hishape;
    //$$ std::cout << sS2Hishapeh(ss1) << "|" << sS2Hishapeh(ss2) << std::endl;
    //std::cout << "ss1=" << ss1 << std::endl;
    //std::cout << "sS2Hishapeh(ss1)=" << sS2Hishapeh(header+"\n"+seq+"\n"+ss2+"\n") << std::endl;
    unsigned found1 = ss1.find('(');
    unsigned found2 = ss2.find('(');
    std::string hishape1="[]", hishape2="[]";
    //std::cout << found1 << "," << found2 << "," << std::string::npos << std::endl;
    if (found1<ss1.length())
    {
	//std::cout << "1true" << std::endl;
	hishape1=sS2Hishapeh(ss1);
    }
    //std::cout << "hishape1=" << hishape1 << std::endl;
    if (found2<ss2.length())
    {
	//std::cout << "2true" << std::endl;
	hishape2=sS2Hishapeh(ss2); 
    }
    //std::cout << "hishape2=" << hishape2 << std::endl;
    generateMatchStr(hishape1, hishape2, &related_hishape);
    //generateMatchStr2(sS2Hishapeh(ss1), sS2Hishapeh(ss2), &related_hishape);
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
    //std::cout << "match_str=" << match_str << std::endl;
    
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
// 					    
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
        //free_path(route);
    }      
    if (printDetails)  printf("Saddle structure: %s\n", struc_max_saddle);
    //printf("%g %g kcal/mol\n", en_max_saddle, en_max_0);
    if (printDetails)  printf("Saddle energy regarding start structure: %g kcal/mol\n", (en_max_saddle-en_max_0));
    return en_max_saddle;
}


/*
#################################
# BEGIN OF FUNCTION DEFINITIONS of findpath.c #
#################################
*/
void free_path(path_t *path){
  path_t *tmp = path;
  if(tmp){
    while(tmp->s){ free(tmp->s); tmp++;}
    free(path);
  }
}

int find_saddle(const char *sequence, const char *struc1, const char *struc2, int max) {
  int maxl, maxE, i;
  const char *tmp;
  move_t *bestpath=NULL;
  int dir;

  path_fwd = 0;
  maxE = INT_MAX - 1;
  seq = sequence;

  update_fold_params();
  make_pair_matrix();

  /* nummerically encode sequence */
  S   = encode_sequence(seq, 0);
  S1  = encode_sequence(seq, 1);

  maxl=1;
  do {
    int saddleE;
    path_fwd = !path_fwd;
    if (maxl>max) maxl=max;
    if(path) free(path);
    //printf("before finding saddle(1675):%s, %s, %d, %d\n", struc1, struc2, maxE, maxl);
    saddleE  = find_path_once(struc1, struc2, maxE, maxl);
    //printf("after finding saddle(1675):%s, %s\n", struc1, struc2);
    if (saddleE<maxE) {
      maxE = saddleE;
      if (bestpath) free(bestpath);
      bestpath = path;
      path = NULL;
      dir = path_fwd;
    } else{
      free(path);path=NULL;
    }
    tmp=struc1;
    struc1=struc2;
    struc2=tmp;
    maxl *=2;
  } while (maxl<2*max);

  free(S); free(S1);
  /* (re)set some globals */
  path=bestpath;
  path_fwd = dir;
  return maxE;
}

void print_path(const char *seq, const char *struc) {
  int d;
  char *s;
  s = strdup(struc);

//REPLACE      
// try {
  printf("%s\n%s %6.2f\n", seq, s, energy_of_alistruct2((const char **)AS,s, n_seq));  //energy_of_struc_on_as(s), energy_of_structure(seq,s, 0);

// } catch (std::exception &e) {
//   std::cerr << "Exception_ali_error_1: " << e.what() << '\n';
//   std::exit(1);
// } 
  qsort(path, BP_dist, sizeof(move_t), compare_moves_when);
  for (d=0; d<BP_dist; d++) {
    int i,j;
    i = path[d].i; j=path[d].j;
    if (i<0) { /* delete */
      s[(-i)-1] = s[(-j)-1] = '.';
    } else {
      s[i-1] = '('; s[j-1] = ')';
    }
//REPLACE
// try {
  printf("%s %6.2f - %6.2f\n", s, energy_of_alistruct2((const char **)AS,s, n_seq), path[d].E/100.0);  // energy_of_struc_on_as(s),energy_of_structure(seq,s, 0)
// } catch (std::exception &e) {
//   std::cerr << "Exception_ali_error_2: " << e.what() << '\n';
//   std::exit(1);
// } 
  }
  
  free(s);
}

path_t *get_path(const char *seq, const char *s1, const char* s2, int maxkeep) {
  int E, d;
  path_t *route=NULL;

  E = find_saddle(seq, s1, s2, maxkeep);

  route = (path_t *)space((BP_dist+2)*sizeof(path_t));

  qsort(path, BP_dist, sizeof(move_t), compare_moves_when);

  if (path_fwd) {
    /* memorize start of path */
    route[0].s  = strdup(s1);
//REPLACE
    //printf("seq=%s\n", seq); Note that seq is consensus sequence
//  try {
    route[0].en = energy_of_alistruct2((const char **)AS,s1, n_seq);  //energy_of_struc_on_as(s1);  //energy_of_structure(seq, s1, 0);
//   } catch (std::exception &e) {
//     std::cerr << "Exception_ali_error_3: " << e.what() << '\n';
//     std::exit(1);
//   } 
    for (d=0; d<BP_dist; d++) {
      int i,j;
      route[d+1].s = strdup(route[d].s);
      i = path[d].i; j=path[d].j;
      if (i<0) { /* delete */
        route[d+1].s[(-i)-1] = route[d+1].s[(-j)-1] = '.';
      } else {
        route[d+1].s[i-1] = '('; route[d+1].s[j-1] = ')';
      }
      route[d+1].en = path[d].E/100.0;
    }
  }
  else {
    /* memorize start of path */

    route[BP_dist].s  = strdup(s2);
//REPLACE
    route[BP_dist].en = energy_of_alistruct2((const char **)AS,s2, n_seq);  //energy_of_struc_on_as(s2);  //energy_of_structure(seq, s2, 0);
    for (d=0; d<BP_dist; d++) {
      int i,j;
      route[BP_dist-d-1].s = strdup(route[BP_dist-d].s);
      i = path[d].i;
      j = path[d].j;
      if (i<0) { /* delete */
        route[BP_dist-d-1].s[(-i)-1] = route[BP_dist-d-1].s[(-j)-1] = '.';
      } else {
        route[BP_dist-d-1].s[i-1] = '('; route[BP_dist-d-1].s[j-1] = ')';
      }
      route[BP_dist-d-1].en = path[d].E/100.0;
    }
  }

#if _DEBUG_FINDPATH_
  fprintf(stderr, "\n%s\n%s\n%s\n\n", seq, s1, s2);
  for (d=0; d<=BP_dist; d++)
    fprintf(stderr, "%s %6.2f\n", route[d].s, route[d].en);
  fprintf(stderr, "%d\n", *num_entry);
#endif

  free(path);path=NULL;
  return (route);
}

int try_moves(intermediate_t c, int maxE, intermediate_t *next, int dist) {
  int *loopidx, len, num_next=0, en, oldE;
  move_t *mv;
  short *pt;

  len = c.pt[0];

  loopidx = pair_table_to_loop_index(c.pt);
  oldE = c.Sen;
  for (mv=c.moves; mv->i!=0; mv++) {
    int i,j;
    if (mv->when>0) continue;
    i = mv->i; j = mv->j;
    pt = (short *) space(sizeof(short)*(len+1));
    //# pt_structure = (char *) space(sizeof(char)*(len+1));
    memcpy(pt, c.pt,(len+1)*sizeof(short));
    //# memcpy(pt_structure, c.
    if (j<0) { /*it's a delete move */
      pt[-i]=0;
      pt[-j]=0;
    } else { /* insert move */
      if ((loopidx[i] == loopidx[j]) && /* i and j belong to same loop */
          (pt[i] == 0) && (pt[j]==0)     /* ... and are unpaired */
          ) {
        pt[i]=j;
        pt[j]=i;
      } else {
        free(pt);
        continue; /* llegal move, try next; */
      }
    }
#ifdef LOOP_EN
    en = c.curr_en + energy_of_move_pt(c.pt, S, S1, i, j);
#else
    //REPLACE                            S is sequence encode with 0, how to transform pairTable to structure (const char *struc)
    en = energy_of_alistruct_pt2((const char **)AS,pt, n_seq);  //energy_of_structure_pt(seq, pt, S, S1, 0);  //energy_of_struc_on_as(pt);
#endif
    if (en<maxE) {
      next[num_next].Sen = (en>oldE)?en:oldE;
      next[num_next].curr_en = en;
      next[num_next].pt = pt;
      mv->when=dist;
      mv->E = en;
      next[num_next++].moves = copy_moves(c.moves);
      mv->when=0;
    }
    else free(pt);
  }
  
  free(loopidx);
  return num_next;
}

int find_path_once(const char *struc1, const char *struc2, int maxE, int maxl) {
  short *pt1, *pt2;
  move_t *mlist;
  int i, len, d, dist=0, result;
  intermediate_t *current, *next;

  pt1 = make_pair_table(struc1);
  pt2 = make_pair_table(struc2);
  len = (int) strlen(struc1);

  mlist = (move_t *) space(sizeof(move_t)*len); /* bp_dist < n */

  /* 
   * save the bp difference of start structure from target structure, 
   * this is driver of actions
   */
  for (i=1; i<=len; i++) {
    if (pt1[i] != pt2[i]) {
      if (i<pt1[i]) { /* need to delete this pair */
        mlist[dist].i = -i;
        mlist[dist].j = -pt1[i];
        mlist[dist++].when = 0;
      }
      if (i<pt2[i]) { /* need to insert this pair */
        mlist[dist].i = i;
        mlist[dist].j = pt2[i];
        mlist[dist++].when = 0;
      }
    }
  }
  free(pt2);
  BP_dist = dist;
  current = (intermediate_t *) space(sizeof(intermediate_t)*(maxl+1));
  current[0].pt = pt1;
  //REPLACE
  current[0].Sen = current[0].curr_en = energy_of_alistruct_pt2((const char **)AS,pt1, n_seq);  //energy_of_structure_pt(seq, pt1, S, S1, 0); //energy_of_struc_on_as(pt1);  //energy_of_structure_pt(seq, pt1, S, S1, 0); 
  current[0].moves = mlist;
  next = (intermediate_t *) space(sizeof(intermediate_t)*(dist*maxl+1));

  for (d=1; d<=dist; d++) { /* go through the distance classes */
    int c, u, num_next=0;
    intermediate_t *cc;

    for (c=0; current[c].pt != NULL; c++) {
      num_next += try_moves(current[c], maxE, next+num_next, d);
    }
    if (num_next==0) {
      for (cc=current; cc->pt != NULL; cc++) free_intermediate(cc);
      current[0].Sen=INT_MAX;
      break;
    }
    /* remove duplicates via sort|uniq
       if this becomes a bottleneck we can use a hash instead */
    qsort(next, num_next, sizeof(intermediate_t),compare_ptable);
    for (u=0,c=1; c<num_next; c++) {
      if (memcmp(next[u].pt,next[c].pt,sizeof(short)*len)!=0) {
        next[++u] = next[c];
      } else {
        free_intermediate(next+c);
      }
    }
    num_next = u+1;
    qsort(next, num_next, sizeof(intermediate_t),compare_energy);
    /* free the old stuff */
    for (cc=current; cc->pt != NULL; cc++) free_intermediate(cc);
    for (u=0; u<maxl && u<num_next; u++) {
      current[u] = next[u];
    }
    for (; u<num_next; u++)
      free_intermediate(next+u);
    num_next=0;
  }
  
  free(next);
  path = current[0].moves;
  result = current[0].Sen;
  free(current[0].pt); free(current);
  return(result);
}


int *pair_table_to_loop_index (short *pt){
  /* number each position by which loop it belongs to (positions start
     at 1) */
  int i,hx,l,nl;
  int length;
  int *stack = NULL;
  int *loop = NULL;

  length = pt[0];
  stack  = (int *) space(sizeof(int)*(length+1));
  loop   = (int *) space(sizeof(int)*(length+2));
  hx=l=nl=0;

  for (i=1; i<=length; i++) {
    if ((pt[i] != 0) && (i < pt[i])) { /* ( */
      nl++; l=nl;
      stack[hx++]=i;
    }
    loop[i]=l;

    if ((pt[i] != 0) && (i > pt[i])) { /* ) */
      --hx;
      if (hx>0)
        l = loop[stack[hx-1]];  /* index of enclosing loop   */
      else l=0;                 /* external loop has index 0 */
      if (hx<0) {
        nrerror("unbalanced brackets in make_pair_table");
      }
    }
  }
  loop[0] = nl;
  free(stack);

#ifdef _DEBUG_LOOPIDX
  fprintf(stderr,"begin loop index\n");
  fprintf(stderr,
          "....,....1....,....2....,....3....,....4"
          "....,....5....,....6....,....7....,....8\n");
  print_structure(pt, loop[0]);
  for (i=1; i<=length; i++)
    fprintf(stderr,"%2d ", loop[i]);
  fprintf(stderr,"\n");
  fprintf(stderr, "end loop index\n");
  fflush(stderr);
#endif

  return (loop);
}

void free_intermediate(intermediate_t *i) {
   free(i->pt);
   free(i->moves);
   i->pt = NULL;
   i->moves = NULL;
   i->Sen = INT_MAX;
 }

int compare_ptable(const void *A, const void *B) {
  intermediate_t *a, *b;
  int c;
  a = (intermediate_t *) A;
  b = (intermediate_t *) B;

  c = memcmp(a->pt, b->pt, a->pt[0]*sizeof(short));
  if (c!=0) return c;
  if ((a->Sen - b->Sen) != 0) return (a->Sen - b->Sen);
  return (a->curr_en - b->curr_en);
}

int compare_energy(const void *A, const void *B) {
  intermediate_t *a, *b;
  a = (intermediate_t *) A;
  b = (intermediate_t *) B;

  if ((a->Sen - b->Sen) != 0) return (a->Sen - b->Sen);
  return (a->curr_en - b->curr_en);
}

int compare_moves_when(const void *A, const void *B) {
  move_t *a, *b;
  a = (move_t *) A;
  b = (move_t *) B;

  return(a->when - b->when);
}


move_t* copy_moves(move_t *mvs) {
  move_t *new_;
  new_ = (move_t *) space(sizeof(move_t)*(BP_dist+1));
  memcpy(new_,mvs,sizeof(move_t)*(BP_dist+1));
  return new_;
}
