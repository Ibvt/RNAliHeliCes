/*
 * main.cc
 *
 *  Created on: Feb 3, 2011
 *      Author: jhuang
 * 
 * Install the core libraries:
 * :~/RNABarrier-distribution/trunk$ sudo apt-get install libboost-all-dev
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


//#include "pf_answer.hh"
//#include "p_func.hh"
//#include "kinetic_hishapem_mfe_pp.hh"
//#include "kinetic_hishapem_mfe_pp.hh"
#include "probs_h_microstate.hh"
#include "probs_hplus_microstate.hh"
#include "probs_m_microstate.hh"
#include "probs_a_microstate.hh"
#include "pfall_microstate.hh"


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


/*
  extent the main.cc with a new option 'P' for probability calculation
*/
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
	("sequence,s", po::value< std::vector<std::string> >(), "Read sequence from specified sequence")
	("type,t", po::value<unsigned int>(&hishape_type)->default_value(2), "Specify hishape abstract type")
	("kbest,k", po::value<unsigned int>(&kbest)->default_value(3), "Choose the k-best classes") 
	("Related,R", po::value<std::string>(&match_str) /*->default_value("4,14.5")*/, "Output related hishapes regarding the given helix indices")
	("exact,e", po::value<unsigned int>(&exact)->implicit_value(1), "Faster computation of exact RNA hishapes")
        ("proba,p", po::value<unsigned int>(&proba)->implicit_value(1), "Calculate accumulated hishape probabilities")
	//("rapid,q", po::value<unsigned int>(&rapid)->implicit_value(1), "Faster computation of exact RNA hishape probabilities")
	("partition,z", po::value<unsigned int>(&partition)->implicit_value(1), "Assign partition functions of hishape classes")
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
	conflicting_options(vm, "file", "sequence");
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
            std::cout << "Usage: RNAalihishapes (-s INPUT|-f INPUT-file) [options]\n";  // [<sequence>|<first hishape>] [<second hishape>]
	    std::cout << "Examples:" << std::endl;
	    std::cout << "  ./RNAalihishapes -s GAGGGCUG_UCUGACCCAGCCCG_____AUGAG___CUGAACCUGCCCUGUGACCUACCGCGUGGCAUUUC#GUGGUCUUUAC_GCCCAAGGGAGGGGUAAUGGGCCUGAGACCGGGCCC_UCGGGC______GUUGAUAACC#GAGGGUAGCUCUGACUCUGCCCG_____AUGAC___ACCAACCUGCCCUGUGACCUAUCGUGUUGAAUUCC# -k 10" << std::endl;
	    std::cout << "  ./RNAalihishapes ../examples/collosoma_slrna.seq -P ./librna/vienna/rna_turner2004.par" << std::endl;
	    std::cout << "  ./RNAalihishapes ../examples/collosoma_slrna.seq -R 27,38 -k 20" << std::endl;
	    std::cout << "  ./RNAalihishapes -f ../examples/collosoma_slrna.seq -k 10 -R 36.5m,41.5,27 -e -t 2 -x0" << std::endl;
            std::cout << desc;
            return 0;
        }

        if (vm.count("version")) {
            std::cout << "RNAalihishapes 1.0 (May 16, 2013)\n";  // TODO: use the C++ buildin date function 
            return 0;
        }
        
        
        //std::string strErrors;
        if (vm.count("file"))
        {
	    std::vector<std::string> files = vm["file"].as< std::vector<std::string> >();
            if (!files.empty()) {
	        // ###################
	        // ## standard mode ##
		if (files.size()==1)
		{
		    std::ifstream file((vm["file"].as< std::vector<std::string> >())[0].c_str());
		    file.exceptions(std::ios_base::badbit |
			std::ios_base::failbit |
			std::ios_base::eofbit);
		    std::filebuf *buffer = file.rdbuf();
		    size_t size = buffer->pubseekoff(0, std::ios::end, std::ios::in);
		    buffer->pubseekpos(0, std::ios::in);
		    input = new char[size+1];
		    assert(input);
		    buffer->sgetn(input, size);
		    input[size] = 0;
		    
		    char *end = input+size;
		    for (char *i = input; i != end; ) {
		      char *s = std::strchr(i, '\n');
		      if (s)
			*s = 0;
		      size_t x = std::strlen(i)+1;
		      char *j = new char[x];    
		      std::strncpy(j, i, x);
		      // convert to small letters
		      std::transform(j,j+x-1,j,ToLower());
		      // t -> u
		      std::replace(j,j+x-1,'t','u');
		      // delete '.'  and '-' from alignment files
		      std::remove(j,j+x-1,'.');
		      std::remove(j,j+x-1,'-');
		      inputs.push_back(std::make_pair(j, x-1));
		      if (s)
			i = s + 1;
		      else
			break;
		    }
		    
		    delete[] input;  //input is different as inputs, inputs are used in each path and it will be deleted in the last line of main()
		    // delete fileinput_char;
		}
		else
		{
		    std::cout << "Status: files.size()=" << files.size() << std::endl;
		    throw gapc::OptException("ERROR: parse options");
		}
	    }
        }

        
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
	    // delete[] input;
	}
	
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
  gapc::class_name obj;
  try {
    obj.init(inputs);
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
  unsigned n = obj.t_0_seq.size();
  for (unsigned int i = 0; ; i+=opts.window_increment) {
    unsigned int right = std::min(n, i+opts.window_size);
    gapc::return_type res = obj.run();
    std::cout << "Answer ("
      << i << ", " << right << ") :\n";
    obj.print_result(std::cout, res);
    for (unsigned int j = 0; j<opts.repeats; ++j)
      obj.print_backtrack(std::cout, res);
    if (i+opts.window_size >= n)
      break;
    obj.window_increment();
  }
#else
  gapc::add_event("start");

  obj.cyk();
  gapc::return_type res = obj.run();

  gapc::add_event("end_computation");

  //std::cout << "Answer: \n";
  //obj.print_result(std::cout, res);
  pfunc = boost::lexical_cast<double>(res);
  //std::cout << pfunc;

  gapc::add_event("end_result_pp");

#ifdef TRACE
  std::cerr << "start backtrack\n";
#endif
  //for (unsigned int i = 0; i<opts.repeats; ++i)
  //  obj.print_backtrack(std::cout, res);
  //obj.print_subopt(std::cout, opts.delta);

  gapc::add_event("end");
#endif
    }



    if (!inputs.empty())
    {
	// ########################
	// ## calculate hishapes ##
	// ########################

  // actual performance gains like 20%
  // see also http://www.ddj.com/cpp/184401305

  // workaround stupid Sun CC std::cout to fd0 after sync_with_stdio 
  // with -m64 and stlport4 bug:
  // http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=6804239
//#if defined(__SUNPRO_CC) && __SUNPRO_CC <= 0x5100
//  #warning Enable sync_with_stdio because of Sun CC 12 Compiler Bug
//#else
//  std::ios_base::sync_with_stdio(false);
//#endif
//  std::cin.tie(0);  
  
	try {
	    std::vector<std::string> tokens;
// 	    if ((proba == 0 && rapid == 0 && partition == 0 && Partition == 0) || rapid==1 || Partition==1) {  // proba == 0 || 
// 		if (hishape_type==4) {
// 		    gapc::hishapeh_pp_cls obj;
// 		    
// 		    obj.init(inputs, /*(int)kbest*1.1*/kbest, match_str, exact, thresh*100, locmin, true);
// 		    
//                     // placeholder
// 		    obj.cyk();  // do nothing
// 		    gapc::hishapeh_pp_ret res = obj.run();
// 		    //## gapc::add_event("end_computation");
// 		    
// 		    obj.print_result(std::cout, res);  // do nothing
// 		    //##gapc::add_event("end_result_pp");
// 		    
// 		    // extends the kbest to 10% 
// 		    //res.set_k(kbest);
// 
// 		    //List_Ref<std::pair<std::pair<Rope, mfeanswer> , intrusive_ptr<Backtrace<String, unsigned int> > > > btp  = obj.bt_proxy_nt_struct();
// 		    //intrusive_ptr<Backtrace<String, unsigned int> >  bt = execute_backtrack_k(btp);
// 		    intrusive_ptr<Backtrace<String, unsigned int> >  bt   = obj.backtrack(obj.t_0_left_most);
// 		    intrusive_ptr<Backtrace_List<String, unsigned int> > l =
// 		      boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
// 		    assert(!bt || (bt && l));
// 		    if (l) {
// 			/* save the data in a vector 'tokens' */
// 			for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
// 			    std::stringstream ss;
// 			    (*i)->print(ss);
// 
// 			    std::copy(std::istream_iterator<std::string>(ss),
// 				      std::istream_iterator<std::string>(),
// 				      std::back_inserter<std::vector<std::string> >(tokens) );
// 			}	 
// 		    }
// 		} else if (hishape_type==3) {
// 		    gapc::hishapehplus_pp_cls obj;
// 		    obj.init(inputs, kbest, match_str, exact, thresh*100, locmin, true);
// 		    obj.cyk();
// 		    gapc::hishapehplus_pp_ret res = obj.run();
// 		    obj.print_result(std::cout, res);   
// 		    //List_Ref<std::pair<std::pair<Rope, mfeanswer> , intrusive_ptr<Backtrace<String, unsigned int> > > > btp  = obj.bt_proxy_nt_struct();
// 		    //intrusive_ptr<Backtrace<String, unsigned int> >  bt = execute_backtrack_k(btp);
// 		    intrusive_ptr<Backtrace<String, unsigned int> >  bt   = obj.backtrack(obj.t_0_left_most);
// 		    intrusive_ptr<Backtrace_List<String, unsigned int> > l =
// 		      boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
// 		    assert(!bt || (bt && l));
// 		    if (l) {
// 			/* save the data in a vector 'tokens' */
// 			for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
// 			    std::stringstream ss;
// 			    (*i)->print(ss);
// 
// 			    std::copy(std::istream_iterator<std::string>(ss),
// 				      std::istream_iterator<std::string>(),
// 				      std::back_inserter<std::vector<std::string> >(tokens) );
// 			}	 
// 		    }
// 		} else if (hishape_type==2) {
// 		    gapc::hishapem_pp_cls obj;
// 		    obj.init(inputs, kbest, match_str, exact, thresh*100, locmin, true);
// 		    obj.cyk();
// 		    gapc::hishapem_pp_ret res = obj.run();
// 		    obj.print_result(std::cout, res);  
// 		    //List_Ref<std::pair<std::pair<Rope, mfeanswer> , intrusive_ptr<Backtrace<String, unsigned int> > > > btp  = obj.bt_proxy_nt_struct();
// 		    //intrusive_ptr<Backtrace<String, unsigned int> >  bt = execute_backtrack_k(btp);
// 		    intrusive_ptr<Backtrace<String, unsigned int> >  bt   = obj.backtrack(obj.t_0_left_most);
// 		    intrusive_ptr<Backtrace_List<String, unsigned int> > l =
// 		      boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
// 		    assert(!bt || (bt && l));
// 		    if (l) {
// 			/* save the data in a vector 'tokens' */
// 			for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
// 			    std::stringstream ss;
// 			    (*i)->print(ss);
// 
// 			    std::copy(std::istream_iterator<std::string>(ss),
// 				      std::istream_iterator<std::string>(),
// 				      std::back_inserter<std::vector<std::string> >(tokens) );
// 			}	 
// 		    } 
// 		} else if (hishape_type==1) {
// 		    gapc::hishapeb_pp_cls obj;
// 		    obj.init(inputs, kbest, match_str, exact, thresh*100, locmin, true);
// 		    obj.cyk();
// 		    gapc::hishapeb_pp_ret res = obj.run();
// 		    obj.print_result(std::cout, res); 
// 		    //List_Ref<std::pair<std::pair<Rope, mfeanswer> , intrusive_ptr<Backtrace<String, unsigned int> > > > btp  = obj.bt_proxy_nt_struct();
// 		    //intrusive_ptr<Backtrace<String, unsigned int> >  bt = execute_backtrack_k(btp);
// 		    intrusive_ptr<Backtrace<String, unsigned int> >  bt   = obj.backtrack(obj.t_0_left_most);
// 		    intrusive_ptr<Backtrace_List<String, unsigned int> > l =
// 		      boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
// 		    assert(!bt || (bt && l));
// 		    if (l) {
// 			/* save the data in a vector 'tokens' */
// 			for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
// 			    std::stringstream ss;
// 			    (*i)->print(ss);
// 
// 			    std::copy(std::istream_iterator<std::string>(ss),
// 				      std::istream_iterator<std::string>(),
// 				      std::back_inserter<std::vector<std::string> >(tokens) );
// 			}	 
// 		    }	 
// 		} else {
// 		    std::cout << "It does not exist hishape abstract type " << hishape_type << "." << std::endl;
// 		}
// 	    }
// 	    else 
	    
	    if (proba==1 || partition == 1) 
	    {
		if (hishape_type==4) {
		    
		    gapc::hishapeh_pfx_cls obj;
		    //TODO: implement pfx with obj.init(inputs, kbest, match_str, exact, thresh*100, locmin, true);
		    obj.init(inputs, kbest, match_str, exact);
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
                        std::cout << ss.str();
			std::copy(std::istream_iterator<std::string>(ss),
				  std::istream_iterator<std::string>(),
				  std::back_inserter<std::vector<std::string> >(tokens) );
		    }
		    
// 		    List_Ref<std::pair<std::pair<Rope, std::pair<mfeanswer, pfanswer> > , intrusive_ptr<Backtrace<String, unsigned int> > > > btp;
// 		    #ifdef WINDOW_MODE
// 			btp  = obj.bt_proxy_nt_struct(obj.t_0_left_most, obj.t_0_right_most);
// 		    #else
// 			btp  = obj.bt_proxy_nt_struct();
// 		    #endif
// 		    intrusive_ptr<Backtrace<String, unsigned int> >  bt   = execute_backtrack_k(btp);
// 		    intrusive_ptr<Backtrace_List<String, unsigned int> > l = boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
// 		    assert(l);
// 		    /// save the data in a vector 'tokens' /
// 		    for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
// 			std::stringstream ss;
// 			(*i)->print(ss);
// 
// 			std::copy(std::istream_iterator<std::string>(ss),
// 				  std::istream_iterator<std::string>(),
// 				  std::back_inserter<std::vector<std::string> >(tokens) );
// 		    }
// 		} 		
// 		else if (hishape_type==3) {
// 		    gapc::hishapehplus_pfx_cls obj;
// 		    obj.init(inputs, kbest, match_str, exact);
// 		    obj.cyk();
// 		    gapc::hishapehplus_pfx_ret res = obj.run();
// 		    obj.print_result(std::cout, res);   
// 		    //======for (unsigned int i = 0; i<opts.repeats; ++i) obj.print_backtrack(std::cout, res);======
// 		    List_Ref<std::pair<std::pair<Rope, std::pair<mfeanswer, pfanswer> > , intrusive_ptr<Backtrace<String, unsigned int> > > > btp;
// intrusive_ptr<Backtrace<std::pair<String, double> , unsigned int> >  bt    = backtrack(t_0_left_most);
// intrusive_ptr<Backtrace_List<std::pair<String, double> , unsigned int> > l =
//   boost::dynamic_pointer_cast<Backtrace_List<std::pair<String, double> , unsigned int> > (bt);
// assert(!bt || (bt && l));
// if (l) {
// for (Backtrace_List<std::pair<String, double> , unsigned int>::iterator i = l->begin();
//      i != l->end(); ++i)
//   (*i)->print(out);
// }
// 		    #ifdef WINDOW_MODE
// 			btp  = obj.bt_proxy_nt_struct(obj.t_0_left_most, obj.t_0_right_most);
// 		    #else
// 			btp  = obj.bt_proxy_nt_struct();
// 		    #endif
// 		    intrusive_ptr<Backtrace<String, unsigned int> >  bt   = execute_backtrack_k(btp);
// 		    intrusive_ptr<Backtrace_List<String, unsigned int> > l = boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
// 		    assert(l);
// 		    /// save the data in a vector 'tokens' /
// 		    for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
// 			std::stringstream ss;
// 			(*i)->print(ss);
// 
// 			std::copy(std::istream_iterator<std::string>(ss),
// 				  std::istream_iterator<std::string>(),
// 				  std::back_inserter<std::vector<std::string> >(tokens) );
// 		    }
// 		} else if (hishape_type==2) {
// 		    gapc::hishapem_pfx_cls obj;
// 		    obj.init(inputs, kbest, match_str, exact);
// 		    obj.cyk();
// 		    gapc::hishapem_pfx_ret res = obj.run();
// 		    obj.print_result(std::cout, res);  
// 		    //======for (unsigned int i = 0; i<opts.repeats; ++i) obj.print_backtrack(std::cout, res);======
// 		    List_Ref<std::pair<std::pair<Rope, std::pair<mfeanswer, pfanswer> > , intrusive_ptr<Backtrace<String, unsigned int> > > > btp;
// 		    #ifdef WINDOW_MODE
// 			btp  = obj.bt_proxy_nt_struct(obj.t_0_left_most, obj.t_0_right_most);
// 		    #else
// 			btp  = obj.bt_proxy_nt_struct();
// 		    #endif
// 		    intrusive_ptr<Backtrace<String, unsigned int> >  bt   = execute_backtrack_k(btp);
// 		    intrusive_ptr<Backtrace_List<String, unsigned int> > l = boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
// 		    assert(l);
// 		    /// save the data in a vector 'tokens' /
// 		    for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
// 			std::stringstream ss;
// 			(*i)->print(ss);
// 
// 			std::copy(std::istream_iterator<std::string>(ss),
// 				  std::istream_iterator<std::string>(),
// 				  std::back_inserter<std::vector<std::string> >(tokens) );
// 		    } 
// 		} else if (hishape_type==1) {
// 		    gapc::hishapeb_pfx_cls obj;
// 		    obj.init(inputs, kbest, match_str, exact);
// 		    obj.cyk();
// 		    gapc::hishapeb_pfx_ret res = obj.run();
// 		    obj.print_result(std::cout, res); 
// 		    //======for (unsigned int i = 0; i<opts.repeats; ++i) obj.print_backtrack(std::cout, res);======
// 		    List_Ref<std::pair<std::pair<Rope, std::pair<mfeanswer, pfanswer> > , intrusive_ptr<Backtrace<String, unsigned int> > > > btp;
// 		    #ifdef WINDOW_MODE
// 			btp  = obj.bt_proxy_nt_struct(obj.t_0_left_most, obj.t_0_right_most);
// 		    #else
// 			btp  = obj.bt_proxy_nt_struct();
// 		    #endif
// 		    intrusive_ptr<Backtrace<String, unsigned int> >  bt   = execute_backtrack_k(btp);
// 		    intrusive_ptr<Backtrace_List<String, unsigned int> > l = boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
// 		    assert(l);
// 		    /// save the data in a vector 'tokens' /
// 		    for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
// 			std::stringstream ss;
// 			(*i)->print(ss);
// 
// 			std::copy(std::istream_iterator<std::string>(ss),
// 				  std::istream_iterator<std::string>(),
// 				  std::back_inserter<std::vector<std::string> >(tokens) );
// 		    }
		} else {
		    std::cout << "It does not exist hishape abstract type " << hishape_type << "." << std::endl;
		}
	    }
	    
	    
	    // #########################
	    // ## processing commonly ##
	    unsigned int i;
	    // print sequence as well as its length
	    std::vector<std::pair<const char*, unsigned> >::iterator it;
	    for ( it=inputs.begin() ; it < inputs.end(); it++ ) {
		std::cout << "length = " << (*it).second << std::endl;
		for (i = 0; i < (*it).second; i++)
		{
		    putchar(toupper((*it).first[i]));
		}
		std::cout << std::endl;
	    }
	    

	    if (proba == 0 /*&& rapid == 0*/ && partition == 0 /*&& Partition == 0*/) 
	    {
		// calculate the longest hishape and mfe
		unsigned int length_longest_hishape = 0, length_hishape = 0, length_longest_mfe = 0, length_mfe = 0;
		for (i=0; 9*i < tokens.size(); i++) {
		    if (tokens[9*i+2]!="_")
			tokens[9*i+2].erase (tokens[9*i+2].end()-1); 
		    
		    length_hishape = tokens[9*i+2].length();
		    length_mfe = tokens[9*i+4].length();
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
		for (i=0; 9*i < tokens.size(); i++) {  
															// #### wrap hishape with "[...]" ####
		    std::cout << tokens[9*i+7] << (boost::format(format_pattern) % ((boost::lexical_cast<int>(tokens[9*i+4]))/100.0f) % ("["+tokens[9*i+2]+"]")).str() << std::endl;
		} 
	    }
	    else if (proba == 1) 
	    {
		// calculate the longest hishape and mfe
		unsigned int length_longest_hishape = 0, length_hishape = 0, length_longest_mfe = 0, length_mfe = 0;
		for (i=0; 13*i < tokens.size(); i++) {
		    if (tokens[13*i+2]!="_")
		        tokens[13*i+2].erase (tokens[13*i+2].end()-1); 
		    
		    length_hishape = tokens[13*i+2].length();
		    length_mfe = tokens[13*i+5].length();
		    if (length_longest_hishape < length_hishape) {
			length_longest_hishape = length_hishape;
		    }
		    if (length_longest_mfe < length_mfe) {
			length_longest_mfe = length_mfe;
		    }
		}      
				      
		// prepare format_pattern and 2 maps to output
		std::string format_pattern = ( boost::format("%%%d.2f %%%ds %%11.6f") % (length_longest_mfe+5) % (length_longest_hishape+5) ).str();
		i = 0;
		for (i=0; 13*i < tokens.size(); i++) {
															  // #### wrap hishape with "[...]" ####
		    std::cout << tokens[13*i+11] << (boost::format(format_pattern) % ((boost::lexical_cast<int>(tokens[13*i+5]))/100.0) % ("["+tokens[13*i+2]+"]") % ((boost::lexical_cast<double>(tokens[13*i+7]))/pfunc)).str() << std::endl;
		}
	    }
	    else if (partition == 1) 
	    {
		// calculate the longest hishape and mfe
		unsigned int length_longest_hishape = 0, length_hishape = 0, length_longest_mfe = 0, length_mfe = 0;
		for (i=0; 13*i < tokens.size(); i++) {
		    if (tokens[13*i+2]!="_")
		        tokens[13*i+2].erase (tokens[13*i+2].end()-1); 
		    
		    length_hishape = tokens[13*i+2].length();
		    length_mfe = tokens[13*i+5].length();
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
		for (i=0; 13*i < tokens.size(); i++) {
															  // #### wrap hishape with "[...]" ####
		    std::cout << tokens[13*i+11] << (boost::format(format_pattern) % ((boost::lexical_cast<int>(tokens[13*i+5]))/100.0) % ("["+tokens[13*i+2]+"]") % (boost::lexical_cast<double>(tokens[13*i+7]))).str() << std::endl;
		}
	    }
// 	    else if (rapid == 1)
// 	    {
// 		// calculate the longest hishape and mfe
// 		unsigned int length_longest_hishape = 0, length_hishape = 0, length_longest_mfe = 0, length_mfe = 0;
// 		for (i=0; 9*i < tokens.size(); i++) {  
// 		    length_hishape = tokens[9*i+2].length();
// 		    length_mfe = tokens[9*i+4].length();
// 		    if (length_longest_hishape < length_hishape) {
// 			length_longest_hishape = length_hishape;
// 		    }
// 		    if (length_longest_mfe < length_mfe) {
// 			length_longest_mfe = length_mfe;
// 		    }
// 		}
// 		
// 		
// 		std::map<std::string, double> hishape_pf;
// 		std::string token_9_2 = "";
// 		// prepare format_pattern and output the results
// 		std::string format_pattern = ( boost::format("%%%d.2f %%%ds %%11.6f") % (length_longest_mfe+5) % (length_longest_hishape+4) ).str();
// 		i = 0;
// 		for (i=0; 9*i < tokens.size(); i++) {   
//                     //std::string hishape = tokens[9*i+2];
//                     double pf = 0.0;
// 		    if (hishape_pf.find(tokens[9*i+2]) != hishape_pf.end()) {
// 			pf = hishape_pf[tokens[9*i+2]];  
// 		    } else {
//                         //match_str = tokens[9*i+2]
//                         std::string match_str_for_function = "";    
// 			std::vector<std::string> helixIndices;
// 			tokenize(tokens[9*i+2], helixIndices);
// 			std::vector<std::string>::iterator it;
// 			for ( it=helixIndices.begin() ; it != helixIndices.end(); it++ )
// 			{
// 			    std::string currentHi = *it;
// 			    if ( (currentHi[currentHi.size()-1] != 'b') && (currentHi[currentHi.size()-1] != 'i') && \
// 				  (currentHi[currentHi.size()-1] != 'm') && (currentHi != "(") && (currentHi != ")") )  // in case h or _
// 			    {
// 				match_str_for_function = match_str_for_function + "," + currentHi;
// 			    }
// 			}
// 
// 		        calculateHishapeProbabilities(inputs,1000000000,hishape_type, match_str_for_function, true, &hishape_pf);
// 			pf = hishape_pf[tokens[9*i+2]]; 
// 			//std::cout << "calculateHishapeProbabilities("<<","<<1000000000<<","<<hishape_type<<","<<match_str_for_function<<")" << std::endl;
// 		    }
// 		    token_9_2 = tokens[9*i+2];
// 		    if (token_9_2 != "_")
// 		        token_9_2.erase (token_9_2.end()-1); 
// 		    std::cout << tokens[9*i+7] << (boost::format(format_pattern) % ((boost::lexical_cast<int>(tokens[9*i+4]))/100.0f) % ("["+token_9_2+"]") % (pf/pfunc)).str() << std::endl;
// 		}
// 	    }
// 	    // output free energies of hishapes
// 	    else if (Partition == 1)
// 	    {
// 	        //_kT = -0.00198717*(273.15 + T);
// 		//std::cout << "T=" << T << ",_kT=" << _kT << std::endl;
// 		// calculate the longest hishape and mfe
// 		unsigned int length_longest_hishape = 0, length_hishape = 0, length_longest_mfe = 0, length_mfe = 0;
// 		for (i=0; 9*i < tokens.size(); i++) {  
// 		    length_hishape = tokens[9*i+2].length();
// 		    length_mfe = tokens[9*i+4].length();
// 		    if (length_longest_hishape < length_hishape) {
// 			length_longest_hishape = length_hishape;
// 		    }
// 		    if (length_longest_mfe < length_mfe) {
// 			length_longest_mfe = length_mfe;
// 		    }
// 		}
// 		
// 		
// 		std::map<std::string, double> hishape_pf;
// 		std::string token_9_2 = "";
// 		// prepare format_pattern and output the results
// 		std::string format_pattern = ( boost::format("%%%d.2f %%%ds %%15.10f") % (length_longest_mfe+5) % (length_longest_hishape+4) ).str();
// 		i = 0;
// 		for (i=0; 9*i < tokens.size(); i++) {   
//                     //std::string hishape = tokens[9*i+2];
//                     double pf = 0.0;
// 		    if (hishape_pf.find(tokens[9*i+2]) != hishape_pf.end()) {
// 			pf = hishape_pf[tokens[9*i+2]];  
// 		    } else {
//                         //match_str = tokens[9*i+2]
//                         std::string match_str_for_function = "";    
// 			std::vector<std::string> helixIndices;
// 			tokenize(tokens[9*i+2], helixIndices);
// 			std::vector<std::string>::iterator it;
// 			for ( it=helixIndices.begin() ; it != helixIndices.end(); it++ )
// 			{
// 			    std::string currentHi = *it;
// 			    if ( (currentHi[currentHi.size()-1] != 'b') && (currentHi[currentHi.size()-1] != 'i') && \
// 				  (currentHi[currentHi.size()-1] != 'm') && (currentHi != "(") && (currentHi != ")") )  // in case h or _
// 			    {
// 				match_str_for_function = match_str_for_function + "," + currentHi;
// 			    }
// 			}
// 
// 		        calculateHishapeProbabilities(inputs,1000000000,hishape_type, match_str_for_function, true, &hishape_pf);
// 			pf = hishape_pf[tokens[9*i+2]]; 
// 			//std::cout << "calculateHishapeProbabilities("<<","<<1000000000<<","<<hishape_type<<","<<match_str_for_function<<")" << std::endl;
// 		    }
// 		    token_9_2 = tokens[9*i+2];
// 		    if (token_9_2 != "_")
// 		        token_9_2.erase (token_9_2.end()-1); 
// 		    std::cout << tokens[9*i+7] << (boost::format(format_pattern) % ((boost::lexical_cast<int>(tokens[9*i+4]))/100.0f) % ("["+token_9_2+"]") % pf).str() << std::endl;
// 		}
// 	    }


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


      /* delete all allocated memory */
      for (inputs_t::iterator i = inputs.begin(); i != inputs.end(); ++i)
	delete[] (*i).first;
  }
  
  //if (input != NULL)
  //  delete[] input;


  return 0;
}
