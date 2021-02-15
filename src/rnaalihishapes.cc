/*
 * main.cc
 *
 *  Created on: Feb 3, 2011
 *      Author: jhuang
 * 
 * Install the core libraries:
 * :~/RNABarrier-distribution/trunk$ sudo apt-get install libboost-all-dev
 * 
 * valgrind --leak-check=full --show-reachable=yes check,
 *     step 1: modify std::ios_base::sync_with_stdio(_false); ==> std::ios_base::sync_with_stdio(_true);
 * valgrind --leak-check=full --show-reachable=yes ./RNAalihishapes -f ../../examples/test4.aln -k 10
 */

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
//#include <functional>
#include <cassert>

extern "C" {
#include "H/aln_util.h"
}
//#include "pf_answer.hh"
//#include "p_func.hh"
//#include "kinetic_hishapem_mfe_pp.hh"
//#include "kinetic_hishapem_mfe_pp.hh"
#include "probs_h_microstate.hh"
#include "probs_hplus_microstate.hh"
#include "probs_m_microstate.hh"
#include "probs_a_microstate.hh"
#include "pfall_microstate.hh"
#include "rep_consensus.hh"


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


#define MAX_NUM_NAMES    500

// hishape abstraction type
const unsigned int HISHAPE_B                = 1;
const unsigned int HISHAPE_M                = 2;
const unsigned int HISHAPE_H_PLUS           = 3;
const unsigned int HISHAPE_H                = 4;


struct ToLower : public std::unary_function<char,char> {
        char operator()(char a)
        {
                return tolower(a);
        };
};

template <typename Value>   void  print_backtrack(std::ostream &out, Value& value)
{

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

// also known as select1st in SGI STL implementation
template<typename T_PAIR>
struct GetKey: public std::unary_function<T_PAIR, typename T_PAIR::first_type>
{
    const typename T_PAIR::first_type& operator()(const T_PAIR& item) const
    {
	return item.first;
    }
};


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

void replace_char (char* s, char find, char replace) {
  while (*s != 0) {
    if (*s == find)
      *s = replace;
    s++;
  }
}

// END FUNCTIONS
//******************************************************************************


int main(int argc, char **argv)
{
  
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
    char *input = 0;                             // whether a sequence is given by user
    unsigned int hishape_type;
    uint32_t kbest;  // , uint32_t kbest
    std::string match_str = "";
    unsigned int exact = 0;
    unsigned int proba = 0;
    //unsigned int rapid = 0;
    unsigned int partition = 0;
    //unsigned int Partition = 0;
    double thresh = 10000.0;
    double T = 37.0;
    std::string P = "";


    try {
      
        // idea is the option is not compulsory and it will be given a relexed code  
        //po::options_description desc("Allowed options", 120, 60);
        po::options_description desc("Allowed options");
	
	desc.add_options()
	#ifdef WINDOW_MODE
		      ("windowsize,w)", po::value<unsigned int>(&window_size)->default_value(0), "Specify window size")
		      ("windowincrement,p", po::value<unsigned int>(&window_increment)->default_value(0), "Specify window position increment (use with -w)")
	#endif
	("file,f", po::value< std::vector<std::string> >(), "Read sequence(s) from specified file")
	//SEQ ("sequence,s", po::value< std::vector<std::string> >(), "Read sequence from specified sequence")
	("type,t", po::value<unsigned int>(&hishape_type)->default_value(2), "Specify hishape abstract type")
	("kbest,k", po::value<unsigned int>(&kbest)->default_value(3), "Choose the k-best classes") 
	("Related,R", po::value<std::string>(&match_str) /*->default_value("4,14.5")*/, "Output related hishapes regarding the given helix indices")
	("exact,e", po::value<unsigned int>(&exact)->implicit_value(1), "Faster computation of exact RNA hishapes")
        ("proba,p", po::value<unsigned int>(&proba)->implicit_value(1), "Calculate accumulated hishape probabilities")
	//("rapid,q", po::value<unsigned int>(&rapid)->implicit_value(1), "Faster computation of exact RNA hishape probabilities")
	("partition,Z", po::value<unsigned int>(&partition)->implicit_value(1), "Assign partition functions of hishape classes")
	//("Partition,Z", po::value<unsigned int>(&Partition)->implicit_value(1), "Assign exact partition functions of hishape classes")
        ("thresh,x", po::value<double>(&thresh)->implicit_value(10000.0), "Specify a threshold value for initstem")
        ("Temperature,T", po::value<double>(&T)->implicit_value(37.0), "Specify the temperature for energy parameters in degree Celsius.") 
	("Parameter,P", po::value<std::string>(&P), "Specify the path of energy parameter file.")
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
	//conflicting_options(vm, "align", "convert");
	//IMPROTANT, FROM LOGIC, TWO ITEMS ARE CONFLICT, BUT ONLY for simplied calculation: conflicting_options(vm, "kbest", "proba"); => using rapidshapes
	//conflicting_options(vm, "proba", "rapid");
	conflicting_options(vm, "proba", "partition");
	//conflicting_options(vm, "proba", "Partition");
	//conflicting_options(vm, "rapid", "partition");
	//conflicting_options(vm, "rapid", "Partition");
	//conflicting_options(vm, "partition", "Partition");
	option_dependency(vm, "exact", "Related");

	
        if (vm.count("help") || vm.size()==2) {
            std::cout << "Usage: RNAlihishapes (-s INPUT|-f INPUT-file) [options]\n";  // [<sequence>|<first hishape>] [<second hishape>]
	    std::cout << "Examples:" << std::endl;
	    //SEQ std::cout << "  ./RNAalihishapes -s GAGGGCUG_UCUGACCCAGCCCG_____AUGAG___CUGAACCUGCCCUGUGACCUACCGCGUGGCAUUUC#GUGGUCUUUAC_GCCCAAGGGAGGGGUAAUGGGCCUGAGACCGGGCCC_UCGGGC______GUUGAUAACC#GAGGGUAGCUCUGACUCUGCCCG_____AUGAC___ACCAACCUGCCCUGUGACCUAUCGUGUUGAAUUCC# -k 10" << std::endl;
	    std::cout << "  ./RNAlihishapes ../examples/test2.aln -P ./librna/vienna/rna_turner1999.par" << std::endl;
	    std::cout << "  ./RNAlihishapes ../examples/test2.aln -R 27,38 -k 20" << std::endl;
	    std::cout << "  ./RNAlihishapes -f ../examples/test2.aln -k 10 -R 36.5m,41.5,27 -e -t 2 -x0" << std::endl;
            std::cout << desc;
            return 0;
        }

        if (vm.count("version")) {
            std::cout << "RNAlihishapes 1.2.1 (June 04, 2020)\n";  // TODO: use the C++ buildin date function 
            return 0;
        }
        
        if (vm.count("file")) {
	    std::vector<std::string> files = vm["file"].as< std::vector<std::string> >();   
	    if (!files.empty()) {    
		if (files.size()==1)
		{
		    int mis=0;
		    FILE     *clust_file = fopen ( files[0].c_str(), "r");  //stdin;
		    
		    int n_seq, length, i;
		    char     *AS[MAX_NUM_NAMES];          /* aligned sequences, a array of sequences */
		    char     *names[MAX_NUM_NAMES];       /* sequence names */

		    n_seq = read_clustal(clust_file, AS, names);  // the output of AS contains '-'

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
		    
		    for (i=0; AS[i]; i++) {
			free(AS[i]); free(names[i]);
		    }
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
	
	if (vm.count("Temperature"))
	{
#ifdef RNALIB_H
            temperature = T;
#endif
	}

	if (vm.count("Parameter"))
	{
#ifdef RNALIB_H
            librna_read_param_file(/*par_filename*/P.c_str());
#endif
	}
	else
	{
#ifdef RNALIB_H
	    librna_read_param_file(0);
#endif
	}
  
    } catch (std::exception &e) {
        std::cerr << "Exception: " << e.what() << '\n';
        std::exit(1);
    }


    double pfunc = 0.0;
    if (proba==1 /*|| rapid==1*/)
    {
  // ########################
  // ## calculate pfunc    ##
  // ########################
  gapc::pfunc_cls pfunc_obj;
  try {
    pfunc_obj.init(inputs);
  } catch (std::exception &e) {
    std::cerr << "Exception: " << e.what() << '\n';
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
  unsigned n = pfunc_obj.t_0_seq.size();
  for (unsigned int i = 0; ; i+=opts.window_increment) {
    unsigned int right = std::min(n, i+opts.window_size);
    gapc::return_type res = pfunc_obj.run();
    std::cout << "Answer ("
      << i << ", " << right << ") :\n";
    pfunc_obj.print_result(std::cout, res);
    for (unsigned int j = 0; j<opts.repeats; ++j)
      pfunc_obj.print_backtrack(std::cout, res);
    if (i+opts.window_size >= n)
      break;
    pfunc_obj.window_increment();
  }
#else
  gapc::add_event("start");

  pfunc_obj.cyk();
  gapc::pfunc_ret pfunc_res = pfunc_obj.run();

  gapc::add_event("end_computation");

  //std::cout << "Answer: \n";
  //pfunc_obj.print_result(std::cout, res);
  pfunc = boost::lexical_cast<double>(pfunc_res);
  //std::cout << pfunc;

  gapc::add_event("end_result_pp");

#ifdef TRACE
  std::cerr << "start backtrack\n";
#endif
  //for (unsigned int i = 0; i<opts.repeats; ++i)
  //  pfunc_obj.print_backtrack(std::cout, res);
  //pfunc_obj.print_subopt(std::cout, opts.delta);

  gapc::add_event("end");
#endif
    }
    
    
    // #########################
    // ## consensus structure ##
    // #########################
    gapc::cons_cls cons_obj;
    try {
      cons_obj.init(inputs);
    } catch (std::exception &e) {
      std::cerr << "Exception: " << e.what() << '\n';
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
    cons_obj.cyk();
    gapc::cons_ret cons_res = cons_obj.run();
    //std::string consensus = boost::lexical_cast<std::string>(cons_res);
  


    if (!inputs.empty())
    {
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
	    std::vector<std::pair<const char*, unsigned> >::iterator it;
	    for ( it=inputs.begin() ; it < inputs.end(); it++ ) {
//		std::cout << "length = " << (*it).second << std::endl;
// 		for (i = 0; i < (*it).second; i++)
// 		{
// 		    putchar(toupper((*it).first[i]));
// 		}
		//std::cout << "concensus sequence " << std::endl;
		std::stringstream cons_ss;
		cons_obj.print_result(cons_ss, cons_res);
		std::string cons_str = cons_ss.str();  
		std::transform(cons_str.begin(), cons_str.end(),cons_str.begin(), ::toupper);
		cons_str.erase(std::remove(cons_str.begin(), cons_str.end(), '\n'), cons_str.end());
		std::cout << cons_str << "    #consensus" << std::endl;
	    }
	    

	    if (proba == 0 /*&& rapid == 0*/ && partition == 0 /*&& Partition == 0*/) 
	    {
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
		    std::cout << tokens[25*i+20] << (boost::format(format_pattern) % ((boost::lexical_cast<double>(tokens[25*i+6]))/100.0f) % ("["+tokens[25*i+2]+"]")).str() << std::endl;
		} 
	    }
	    else if (proba == 1) 
	    {
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
		    //std::cout << tokens[25*i+2] << "," << tokens[25*i+6] << "," << tokens[25*i+15] << "," << tokens[25*i+20] << std::endl;
		}      
				      
		// prepare format_pattern and 2 maps to output
		std::string format_pattern = ( boost::format("%%%d.2f %%%ds %%11.6f") % (length_longest_mfe+5) % (length_longest_hishape+5) ).str();
		i = 0;
		for (i=0; 25*i < tokens.size(); i++) {
															  // #### wrap hishape with "[...]" ####
		    std::cout << tokens[25*i+20] << (boost::format(format_pattern) % ((boost::lexical_cast<double>(tokens[25*i+6]))/100.0f) % ("["+tokens[25*i+2]+"]") % ((boost::lexical_cast<double>(tokens[25*i+15]))/pfunc)).str() << std::endl;
		}
	    }
	    else if (partition == 1) 
	    {
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
		     
		// prepare format_pattern and 2 maps to output
		std::string format_pattern = ( boost::format("%%%d.2f %%%ds %%15.10f") % (length_longest_mfe+5) % (length_longest_hishape+5) ).str();
		i = 0;
		for (i=0; 25*i < tokens.size(); i++) {
															  // #### wrap hishape with "[...]" ####
		    std::cout << tokens[25*i+20] << (boost::format(format_pattern) % ((boost::lexical_cast<double>(tokens[25*i+6]))/100.0) % ("["+tokens[25*i+2]+"]") % (boost::lexical_cast<double>(tokens[25*i+15]))).str() << std::endl;
		}
	    }                
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
  
  /* delete all allocated memory */
  //for (inputs_t::iterator i = inputs.begin(); i != inputs.end(); ++i)
  //  delete (*i).first;
  if (input != NULL)
    free(input);


  return 0;
}
