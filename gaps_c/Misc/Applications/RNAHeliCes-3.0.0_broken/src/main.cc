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
#include "p_func.hh"
//#include "kinetic_hishapem_mfe_pp.hh"
//#include "kinetic_hishapem_mfe_pp.hh"
#include "hishapeh_mfeV4_pp.hh"
#include "hishapehplus_mfeV4_pp.hh"
#include "hishapem_mfeV4_pp.hh"
#include "hishapeb_mfeV4_pp.hh"
#include "hishapeh_mfe_pfx.hh"
#include "hishapehplus_mfe_pfx.hh"
#include "hishapem_mfe_pfx.hh"
#include "hishapeb_mfe_pfx.hh"

#if 0
#include "p_func_04.hh"
#include "hishapeh_mfe_pp_04.hh"
#include "hishapehplus_mfe_pp_04.hh"
#include "hishapem_mfe_pp_04.hh"
#include "hishapeb_mfe_pp_04.hh"
#include "hishapeh_mfe_pfx_04.hh"
#include "hishapehplus_mfe_pfx_04.hh"
#include "hishapem_mfe_pfx_04.hh"
#include "hishapeb_mfe_pfx_04.hh"
#endif


#include "rtlib/string.hh"
#include "rtlib/list.hh"
#include "rtlib/hash.hh"
#include "rtlib/asymptotics.hh"
#include "rtlib/generic_opts.hh"

#include "comp/DebugTools.hpp"
#include "comp/Hishape.hpp"
#include "comp/HelixIndex.hpp"
#include "comp/ShapiroAligner.hpp"
#include "comp/NASecondaryStructure.hpp"
#include "comp/NASSDotBracketLoader.hpp"
#include "comp/TreeNode.hpp"
#include "comp/HelixIndex.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/replace.hpp>

namespace po = boost::program_options;

//using namespace std;

///////////////////////////////////////////////////////////////////
/////////////// post processing of stable local minima ////////////
#include "RNAStackFold.hh"
#include "clusterRS.hh"
#include "neighborjoin.hh"
//#include <sys/time.h>
//#include "ComputeE.cc"


//static double   _kT = -1.;
//std::string infile, outfile;//, logfile;
//std::string configfile("Config/config");
Config cf;

//int knum = 5;
bool is_barriercutoff = false;
double barriercutoff = 0;
//bool is_barriercutoff_percent = false;
//double barriercutoff_percent = 0.38;
//bool isdelta = false;
//bool isdelta_p = false;
//int delta = 1.0;
//int delta_p = 20;
//bool reqLocOpt = true;
TurnerE TE;

//int MINMAXOPT = 1;//by default, MIN.
//bool is_onlyLOpt = false;//false=output both LOpt and SLOpt stack configurations, true=only output LOpt stack configurations
//vvd is a lower triangular matrix, each entry (i,j) of which represents the conflicting degree between i and j
//cutoff is the new cut off, for each entry (i,j) with value less or equal to cutoff, we will merge column i and j to i, since the free energy of Si is always less or equal to Sj
//tree[i] = j means structrue Si is merged to Sj, if j==-1, Si is not merged. 
//cluster contains the indices of structures to consider. 
void applycutoff(std::vector <std::vector <double> > & vvd, double & cutoff, std::vector<TNode> & tree, std::vector<int> & cluster) {
	int n = sz(vvd);
	assert(n == sz(cluster));
	std::vector<int> newcluster;
	std::vector<std::vector<double> > newvvd;
	int m = 0;
	FOR(i,0,n) {
		bool merge = false;
		for(int j = sz(vvd[i])-1; j>=0; j--) {
			if(tree[j].pt == -1) {
				if(vvd[i][j] <= cutoff) {
					//merge j to i
					merge = true;
					tree[i].pt = j;
					break;
				}
			}
		}
		if(!merge) {
			m++; // the size of newcluster
			newcluster.PB(cluster[i]);
			assert(m==sz(newcluster));
			std::vector<double> vd;
			newvvd.PB(vd);
			assert(sz(newvvd)==m);
			FOR(j,0,m-1) {
				int k = newcluster[j];
				newvvd[m-1].PB(vvd[i][k]);
			}

		}
	}
	cluster = newcluster;
	vvd = newvvd;
}

void print_phylipltm(const std::vector<std::vector<double> > & vvd) {
//DEBUG	std::cout << "   " <<sz(vvd)<<std::endl;
	FOR(i,0,sz(vvd)) {
		std::cout<<"S"/*<<left<<std::setw(12)<<i*/;
		printf("%10d",i);
		assert(sz(vvd[i]) == i);
// 		FOR(j,0,sz(vvd[i])) {
// 			//std::cout/*<<left*/<<std::setw(6)<<std::setprecision(3)<<vvd[i][j];
// 			printf("%10.4g",vvd[i][j]);
// 		}
		FOR(j,0,sz(vvd)) {
		    if (j<sz(vvd[i])) {
			printf("%10.4g",vvd[i][j]);
		    } else if (j==sz(vvd[i])) {
		        printf("%10.4g",0.0);
	            } else {
		        printf("%10.4g",vvd[j][i]);
		    }
		}
		std::cout<<std::endl;
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


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

// ################ functions for hishape comparison #################
// the third parameter of from_string() should be 
// one of std::hex, std::dec or std::oct
template <class T>
bool from_string(T& t, 
                 const std::string& s, 
                 std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}

typedef val::biocpp::tools::CTreeNode<hishape::comp::CHelixIndex> CHelixIndexTreeNode;
//! A tree node that contains a pair of nucleotides.
CHelixIndexTreeNode::CTreeNodePointer HishapeToTree(std::string strHishape)
{
	strHishape = strHishape.substr(1,strHishape.size()-2);
	//std::cout << "processed strHishape=" << strHishape << std::endl;
	

	// ###############################
	// #### prepare mmBracketPair ####
	std::multimap<int,int> mmBracketPair;
	std::multimap<int,int>::iterator itBracketPair;
	
	std::vector<std::string> helixIndices;
	tokenize(strHishape, helixIndices);  
	std::stack<int> helixIndexStack;
	int iHishapeLength = (signed)helixIndices.size();
	//BUG_TRACK  std::cout << iHishapeLength << std::endl;
	int i = 0;
	while (i < iHishapeLength)
	{
		if ("(" == helixIndices[i])
		{
			helixIndexStack.push(i++);
		}
		else if (")" == helixIndices[i])
		{
			if (true == helixIndexStack.empty())
			{
				throw std::invalid_argument("SNASSDotBracketLoader: The given sequence does not contain matching brackets!");
			}
			mmBracketPair.insert(std::pair<int,int>(helixIndexStack.top()-1,i));  // for m ... )
			//mmBracketPair.insert(std::pair<int,int>(helixIndexStack.top(),i));    // for ( ... )
			helixIndexStack.pop();
			i++;
		}
		else
		{
		        i++;
			//mmBracketPair.insert(std::pair<int,int>(i++,-1));
		}
	}    
	std::vector<std::string>::iterator it;
	for ( it=helixIndices.begin() ; it != helixIndices.end(); it++ )
	{
	    std::string currentHi = *it;
	    if ( (currentHi[currentHi.size()-1] != 'b') && (currentHi[currentHi.size()-1] != 'i') && \
		 (currentHi[currentHi.size()-1] != 'm') && (currentHi != "_") && (currentHi != "(") && (currentHi != ")") )  // in case h
	    {
	        //BUG_TRACK  std::cout << " " << *it;
                (*it) = (*it) + "h";
	    }
	}
	    


	// show content: 算出所有的 bracket对的位置
	//BUG_TRACK  for ( itBracketPair=mmBracketPair.begin() ; itBracketPair != mmBracketPair.end(); itBracketPair++ )
	//BUG_TRACK      std::cout << (*itBracketPair).first << " => " << (*itBracketPair).second << std::endl;
	



	
	//CHelixIndexTreeNode::CTreeNodePointer ptNode = NULL;
	//std::stack<std::pair<int, CHelixIndexTreeNode::CTreeNodePointer> > tStack;
	//int iHishapeLength = tHishape.GetLength();
	//int iIsolatedIndex = 0;
	
	
	

	// node root
	hishape::comp::CHelixIndex tHelixIndex(0.0f, 'N');
	CHelixIndexTreeNode::CTreeNodePointer rootNode = new CHelixIndexTreeNode(tHelixIndex);
	CHelixIndexTreeNode::CTreeNodePointer actualNode = rootNode;
	
	std::deque<int> mlPositionDeque;
	std::deque<CHelixIndexTreeNode::CTreeNodePointer> mlPointerDeque;
	
	// add the first level
	for (i = 0; i < iHishapeLength; i++)
	{
	    //std::cout << helixIndices[i] << std::endl;
	    std::string currentHi = helixIndices[i];
	    //std::cout << currentHi.substr (0,currentHi.size()-1) << std::endl;
	    char lastChar = currentHi.at(currentHi.size()-1);
	    float f;
	    if (lastChar == 'b' || lastChar == 'i')
	    {
	        //assert (from_string<float>(f, currentHi.substr (0,currentHi.size()-1), std::dec));
	        if(from_string<float>(f, currentHi.substr (0,currentHi.size()-1), std::dec))
                {
	            hishape::comp::CHelixIndex tHelixIndex(f, lastChar);
	            actualNode = actualNode->AddChild(tHelixIndex);
		}
		else
		{   
		    throw std::invalid_argument("Error occurs during transformation from string to float!");
		}
	    }
	    else if (lastChar == 'h')
	    {
	        if(from_string<float>(f, currentHi.substr (0,currentHi.size()-1), std::dec))
                {
	            hishape::comp::CHelixIndex tHelixIndex(f, lastChar);
	            actualNode->AddChild(tHelixIndex);
		}
		else
		{
		    throw std::invalid_argument("Error occurs during transformation from string to float!");
		}
		actualNode = rootNode;
	    }
	    else if (lastChar == 'm')
	    {
	        if(from_string<float>(f, currentHi.substr (0,currentHi.size()-1), std::dec))
                {
	            hishape::comp::CHelixIndex tHelixIndex(f, lastChar);
	            actualNode = actualNode->AddChild(tHelixIndex);
		}
		else
		{
		    throw std::invalid_argument("Error occurs during transformation from string to float!");
		}
		mlPositionDeque.push_back(i);
		mlPointerDeque.push_back(actualNode);
		actualNode = rootNode;
		i = mmBracketPair.find(i)->second;
	    }
	    
	}
	
	/*
	int mlPosition = mlPositionDeque.front();
	mlPositionDeque.pop_front();
	std::cout << mlPosition << std::endl;
	CHelixIndexTreeNode::CTreeNodePointer mlPointer = mlPointerDeque.front();
	mlPointerDeque.pop_front();
	hishape::comp::CHelixIndex tempHelixIndex(100.0f, 'z');
	mlPointer->AddChild(tempHelixIndex);
	*/

	// for further recursive calculation
	while (false == mlPositionDeque.empty())
	{
	    int mlPosition = mlPositionDeque.front();
	    mlPositionDeque.pop_front();
	    CHelixIndexTreeNode::CTreeNodePointer localRootNode = mlPointerDeque.front();
	    actualNode = localRootNode;
	    mlPointerDeque.pop_front();
	    int start = mlPosition + 1;                          // position of bracket open
	    int end = mmBracketPair.find(mlPosition)->second;    // position of bracket close
	    
	    for (i = (start+1); i < end; i++)
	    {
		//std::cout << helixIndices[i] << std::endl;
		std::string currentHi = helixIndices[i];
		//std::cout << currentHi.substr (0,currentHi.size()-1) << std::endl;
		char lastChar = currentHi.at(currentHi.size()-1);
		float f;
		if (lastChar == 'b' || lastChar == 'i')
		{
		    //assert (from_string<float>(f, currentHi.substr (0,currentHi.size()-1), std::dec));
		    if(from_string<float>(f, currentHi.substr (0,currentHi.size()-1), std::dec))
		    {
			hishape::comp::CHelixIndex tHelixIndex(f, lastChar);
			actualNode = actualNode->AddChild(tHelixIndex);
		    }
		    else
		    {   
			throw std::invalid_argument("Error occurs during transformation from string to float!");
		    }
		}
		else if (lastChar == 'h')
		{
		    if(from_string<float>(f, currentHi.substr (0,currentHi.size()-1), std::dec))
		    {
			hishape::comp::CHelixIndex tHelixIndex(f, lastChar);
			actualNode->AddChild(tHelixIndex);
		    }
		    else
		    {
			throw std::invalid_argument("Error occurs during transformation from string to float!");
		    }
		    actualNode = localRootNode;
		}
		else if (lastChar == 'm')
		{
		    if(from_string<float>(f, currentHi.substr (0,currentHi.size()-1), std::dec))
		    {
			hishape::comp::CHelixIndex tHelixIndex(f, lastChar);
			actualNode = actualNode->AddChild(tHelixIndex);
		    }
		    else
		    {
			throw std::invalid_argument("Error occurs during transformation from string to float!");
		    }
		    mlPositionDeque.push_back(i);
		    mlPointerDeque.push_back(actualNode);
		    actualNode = localRootNode;
		    i = mmBracketPair.find(i)->second;
		}
	    }
	}
	
	
	CHelixIndexTreeNode::CTreeNodePointerList tSuffixNodesList1 = rootNode->GetSuffixList();
	////std::cout << "Suffix order for tree 1:" << std::endl;
	std::string strSuffix1("[");
	for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tSuffixIter1 = tSuffixNodesList1.begin(); tSuffixIter1 != tSuffixNodesList1.end(); tSuffixIter1++)
	{
		strSuffix1 += (*tSuffixIter1)->GetObject().ToString() /*GetPairBaseCodes()*/ + ", ";
	}
	strSuffix1.resize((signed)strSuffix1.size() - 2);
	strSuffix1 += "]";
	//BUG_TRACK  std::cout << strSuffix1 << std::endl;

	return rootNode;
	
}


std::string PrintStemLoop(const val::biocpp::CNASecondaryStructure::SSubstructureData &tRetData, const val::biocpp::CNASecondaryStructure &tSecStructure, unsigned int hishape_type)
{
        std::ostringstream outs;
	int iPosition=tRetData.iStartIndex1;
	for ( ; tSecStructure.GetNucleotideConfig(iPosition+1)!=val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP; iPosition++)
	{       
                // because the hairpin.size() >= 3, the position (iPosition+2) is always safe
	        //if (tSecStructure.GetNucleotideConfig(iPosition+2)==val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP)
	        //{
		//    break;
		//}
		if ( tSecStructure.IsPaired(iPosition) && (!tSecStructure.IsPaired(iPosition+1) || !tSecStructure.IsPaired(tSecStructure.GetPairedBaseIndex(iPosition)-1)) ) 
		{
		      float helixIndex = (iPosition + tSecStructure.GetPairedBaseIndex(iPosition))/2.0f;
		      int checkPosition = iPosition+1;
		      if (tSecStructure.IsPaired(checkPosition))
		      {
			  checkPosition = tSecStructure.GetPairedBaseIndex(iPosition)-1;
		      }
		      assert (!tSecStructure.IsPaired(checkPosition));

		      switch (tSecStructure.GetNucleotideConfig(checkPosition))
		      {
			    case val::biocpp::na_secondary_structure::I_EXTERNAL_BASE:
			    {
				    outs << "Error_I_EXTERNAL_BASE";
				    break;
			    }
			    case val::biocpp::na_secondary_structure::I_STACKED_PAIR:
			    {
				    outs << "Error_I_STACKED_PAIR";
				    break;
			    }
			    case val::biocpp::na_secondary_structure::I_BULGE:
			    {
			            if (hishape_type == HISHAPE_B)
				        outs << helixIndex << "b,";
				    break;
			    }
			    case val::biocpp::na_secondary_structure::I_INTERNAL_LOOP:
			    {
			            if (hishape_type == HISHAPE_B)
				        outs << helixIndex << "i,";
				    break;
			    }
			    case val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP:
			    {
				    outs << "Error_I_HAIRPIN_LOOP";
				    break;
			    }
			    case val::biocpp::na_secondary_structure::I_MULTIBRANCHED_LOOP:
			    {
				    outs << "Error_I_MULTIBRANCHED_LOOP";
				    break;
			    }
			    default:
			    {
				    outs << "Error_" << tSecStructure.GetNucleotideConfig(iPosition+1);
			    }
		      }
		}
	}
	assert (tSecStructure.IsPaired(iPosition));
	outs << (iPosition+tSecStructure.GetPairedBaseIndex(iPosition))/2.0f << ","; 
	return outs.str();
}


std::string PrintStem(const val::biocpp::CNASecondaryStructure::SSubstructureData &tRetData, const val::biocpp::CNASecondaryStructure &tSecStructure, unsigned int hishape_type)
{
        std::ostringstream outs;
	for (int iPosition=tRetData.iStartIndex1; iPosition<tRetData.iEndIndex1; iPosition++)
	{
		if ( tSecStructure.IsPaired(iPosition) && (!tSecStructure.IsPaired(iPosition+1) || !tSecStructure.IsPaired(tSecStructure.GetPairedBaseIndex(iPosition)-1)) ) 
		{
		      float helixIndex = (iPosition + tSecStructure.GetPairedBaseIndex(iPosition))/2.0f;
		      int checkPosition = iPosition+1;
		      if (tSecStructure.IsPaired(checkPosition))
		      {
			  checkPosition = tSecStructure.GetPairedBaseIndex(iPosition)-1;
		      }
		      assert (!tSecStructure.IsPaired(checkPosition));


		      switch (tSecStructure.GetNucleotideConfig(checkPosition))
		      {
			    case val::biocpp::na_secondary_structure::I_EXTERNAL_BASE:
			    {
				    outs << "Error_I_EXTERNAL_BASE";
				    break;
			    }
			    case val::biocpp::na_secondary_structure::I_STACKED_PAIR:
			    {
				    outs << "Error_I_STACKED_PAIR";
				    break;
			    }
			    case val::biocpp::na_secondary_structure::I_BULGE:
			    {
			            if (hishape_type == HISHAPE_B)
				        outs << helixIndex << "b,";
				    break;
			    }
			    case val::biocpp::na_secondary_structure::I_INTERNAL_LOOP:
			    {
			            if (hishape_type == HISHAPE_B)
				        outs << helixIndex << "i,";
				    break;
			    }
			    case val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP:
			    {
				    outs << "Error_I_HAIRPIN_LOOP";
				    break;
			    }
			    case val::biocpp::na_secondary_structure::I_MULTIBRANCHED_LOOP:
			    {
				    outs << "Error_I_MULTIBRANCHED_LOOP";
				    break;
			    }
			    default:
			    {
				    outs << "Error_" << tSecStructure.GetNucleotideConfig(iPosition+1);
			    }
		      }
		}
	}
	if (hishape_type <= HISHAPE_M)
	{
	    outs << (tRetData.iEndIndex1+tRetData.iStartIndex2)/2.0f << 'm' << ",";
	}
	return outs.str();
}
	

typedef val::biocpp::tools::CTreeNode<hishape::comp::CHelixIndex> CHelixIndexTreeNode;
std::string PrintRecursive( const CHelixIndexTreeNode::CTreeNodePointer &treeNode, const val::biocpp::CNASecondaryStructure &tSecStructure, const std::multimap<int,int> &mmFirstPosJ, unsigned int hishape_type)
{
	CHelixIndexTreeNode::CTreeNodePointerList suffixList_Children = treeNode->GetDirectChildrenList();
	CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter = suffixList_Children.begin();
        std::ostringstream outs;
	
	while (tIter != suffixList_Children.end())
        {
		
		int j = mmFirstPosJ.find((*tIter)->GetObject().GetHiNumber())->second;
		val::biocpp::CNASecondaryStructure::SSubstructureData tRetData = tSecStructure.GetStemOrStemLoopData(j);
		if (tRetData.iStructureNature == val::biocpp::na_secondary_structure::I_STEM_TYPE)  // knot m
		{
		    outs << PrintStem(tRetData, tSecStructure, hishape_type);
		    if (hishape_type <= HISHAPE_H_PLUS)    outs << "(,";
		    outs << PrintRecursive(*tIter, tSecStructure, mmFirstPosJ, hishape_type);
		    if (hishape_type <= HISHAPE_H_PLUS)    outs << "),";
		}   
		else    // leaf
		{
		    outs << PrintStemLoop(tRetData, tSecStructure, hishape_type);
		    outs << PrintRecursive(*tIter, tSecStructure, mmFirstPosJ, hishape_type);
		}

		tIter++;
	}
	return outs.str();
}

// transform secondary structure to hishape
std::string sSToHishape(const std::string &strStructure1, int hishape_type)
{
	val::biocpp::loaders::SNASSDotBracketLoader tDotBracketLoader;
	// TODO: use class CStem instead of class CHelixIndex in Tree
	// ## reuse following data structures
	std::multimap<int,int> mmFirstPosJ;
	std::multimap<int,int>::iterator itFirstPosJ;
	
	
	// ###################################
	// #### convert SS1 to Hishape1 ####
	// ###################################
	val::biocpp::CNASecondaryStructure tSecStructure1 = tDotBracketLoader(strStructure1);
	
	// #####################################################
	// #### get the tree structure of StemOrStemLoopData
	
	CHelixIndexTreeNode::CTreeNodePointer ptRoot1 = tSecStructure1.ParseStructureWithTree();
	CHelixIndexTreeNode::CTreeNodePointerList rootSuffix1 = ptRoot1->GetSuffixList();
	
	// #####################################################
	// #### prepare a map(firstPosition, j)
	int iSecondStructureStemsCount1 = tSecStructure1.GetStemsAndStemLoopsCount();
	val::biocpp::CNASecondaryStructure::SSubstructureData tRetData1;
	
	for (int j1=1; j1 <= iSecondStructureStemsCount1; ++j1)
	{
		tRetData1 = tSecStructure1.GetStemOrStemLoopData(j1);
		mmFirstPosJ.insert (std::pair<int,int>(tRetData1.iStartIndex1,j1));
	}
	
	// #####################################################
	// #### first round and further recursive
	std::ostringstream stmHishape;
	//if (hishape_type <= HISHAPE_H_PLUS)    std::cout << "(,";
	stmHishape << "[";
	CHelixIndexTreeNode::CTreeNodePointerList ptRootChildren1 = ptRoot1->GetDirectChildrenList();
	for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter1 = ptRootChildren1.begin(); tIter1 != ptRootChildren1.end(); tIter1++)
	{
	    //if (hishape_type <= HISHAPE_H_PLUS)    std::cout << "(,";
	    int j1 = mmFirstPosJ.find((*tIter1)->GetObject().GetHiNumber())->second;
	    val::biocpp::CNASecondaryStructure::SSubstructureData tRetData1 = tSecStructure1.GetStemOrStemLoopData(j1);
	    if (tRetData1.iStructureNature == val::biocpp::na_secondary_structure::I_STEM_TYPE)  // knot m
	    {
		    stmHishape << PrintStem(tRetData1, tSecStructure1, hishape_type);
		    if (hishape_type <= HISHAPE_H_PLUS)    stmHishape << "(,";
		    stmHishape << PrintRecursive(*tIter1, tSecStructure1, mmFirstPosJ, hishape_type);
		    if (hishape_type <= HISHAPE_H_PLUS)    stmHishape << "),";
	    }   
	    else    // leaf
	    {
		stmHishape << PrintStemLoop(tRetData1, tSecStructure1, hishape_type);
		stmHishape << PrintRecursive(*tIter1, tSecStructure1, mmFirstPosJ, hishape_type);
	    }
	    //if (hishape_type <= HISHAPE_H_PLUS)    std::cout << "),";
	}
	//if (hishape_type <= HISHAPE_H_PLUS)    std::cout << "),";
	std::string strHishape = stmHishape.str();
	if (!strHishape.empty())
	{
	    strHishape = strHishape.substr(0,strHishape.size()-1);
	}
	
	strHishape = strHishape + "]";
	delete ptRoot1;
	return strHishape; 
}

int pairAlign(std::string strInput1, std::string strInput2, int hishape_type, float ratio)
{
	// Code for Trim trailing Spaces only
	size_t endpos = strInput1.find_last_not_of("\n"); // Find the first character position from reverse
	if ( std::string::npos != endpos )
	    strInput1 = strInput1.substr( 0, endpos+1 );
	endpos = strInput2.find_last_not_of("\n");
	if ( std::string::npos != endpos )
	    strInput2 = strInput2.substr( 0, endpos+1 );
	
	try
	{
		// #####################################################################
		// #### transfer a hishape type h_plus, m, b into the tree ####
		//BUG_TRACK  std::cout << "strInput1=" << strInput1 << std::endl;
		//BUG_TRACK  std::cout << "strInput2=" << strInput2 << std::endl;
		//BUG_TRACK  std::cout << "hishape_type=" << hishape_type << std::endl;
		if (hishape_type==HISHAPE_H_PLUS) {
		    boost::replace_all(strInput1, "(", "-1m,(");
		    boost::replace_all(strInput2, "(", "-1m,(");
		}
		//BUG_TRACK  std::cout << "strInput1=" << strInput1 << std::endl;
		//BUG_TRACK  std::cout << "strInput2=" << strInput2 << std::endl;
		CHelixIndexTreeNode::CTreeNodePointer ptTree1 = HishapeToTree(strInput1);
		CHelixIndexTreeNode::CTreeNodePointer ptTree2 = HishapeToTree(strInput2);
		
		if (ratio < 0.0f || ratio > 10000.0f) 
		{
		    std::cout << "The ratio value should be specified between 0.0 and 10000.0!" << std::endl;
		    return 1;
		}
		hishape::comp::algorithms::CShapiroAligner tShapiroAligner(ptTree1, ptTree2, ratio);
		hishape::comp::CScoreScheme tScoreScheme;
		tShapiroAligner.SetScoreScheme(tScoreScheme);
		tShapiroAligner.Align();
		// TODO: fitting the score calculation
		//tShapiroAligner.SetFirstTree(HishapeToTree("[27]"));
		//tShapiroAligner.SetSecondTree(HishapeToTree("[36.5m,(,27,38,)]"));
		//tShapiroAligner.Align(false,false);
		
		//BUG_TRACK  std::cout << "aligning finished!" << std::endl;
		std::cout << tShapiroAligner.GetScore() << std::endl;                     // delete tShapiroAligner?? 
		return 0;
	}
	catch (std::exception &ex)
	{
		std::cout << "ERROR: " << ex.what() << std::endl;
	}
	catch (...)
	{
		std::cout << "ERROR: unexpected exception!" << std::endl;
	}
}

// method for multialign calculation
// TODO: simplify the following function
// TODO: return a map means this is a deep copy, potential improvement
void calculateHishapes(const std::vector<std::pair<const char*, unsigned> > &inputs, uint32_t kbest, unsigned int hishape_type, std::map<std::string,std::string> *hishape_ppmfehishape, bool exact, double thresh, bool locmin)
{
	try {
	    std::vector<std::string> tokens;

	    if (hishape_type==4) {
		gapc::hishapeh_pp_cls obj;
		obj.init(inputs, kbest, std::string(""), exact, thresh*100, locmin);
		obj.cyk();
		gapc::hishapeh_pp_ret res = obj.run();
		obj.print_result(std::cout, res);
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
		obj.init(inputs, kbest, std::string(""), exact, thresh*100, locmin);
		obj.cyk();
		gapc::hishapehplus_pp_ret res = obj.run();
		obj.print_result(std::cout, res);   
		//List_Ref<std::pair<std::pair<Rope, mfeanswer> , intrusive_ptr<Backtrace<String, unsigned int> > > > btp  = obj.bt_proxy_nt_struct();
		//intrusive_ptr<Backtrace<String, unsigned int> >  bt = execute_backtrack_k(btp);
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
		obj.init(inputs, kbest, std::string(""), exact, thresh*100, locmin);
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
//## gapc::add_event("start");
		obj.cyk();
		gapc::hishapem_pp_ret res = obj.run();
		obj.print_result(std::cout, res);  
		//List_Ref<std::pair<std::pair<Rope, mfeanswer> , intrusive_ptr<Backtrace<String, unsigned int> > > > btp  = obj.bt_proxy_nt_struct();
		//intrusive_ptr<Backtrace<String, unsigned int> >  bt = execute_backtrack_k(btp);
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
		obj.init(inputs, kbest, std::string(""), exact, thresh*100, locmin);
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
//## gapc::add_event("start");
		obj.cyk();
		gapc::hishapeb_pp_ret res = obj.run();
		obj.print_result(std::cout, res); 
		//List_Ref<std::pair<std::pair<Rope, mfeanswer> , intrusive_ptr<Backtrace<String, unsigned int> > > > btp  = obj.bt_proxy_nt_struct();
		//intrusive_ptr<Backtrace<String, unsigned int> >  bt = execute_backtrack_k(btp);
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
		std::cout << "It is not allowed to align hishapes with abstract type " << hishape_type << "." << std::endl;
	    }
		
	        
	    // #########################
	    // ## processing commonly ##
	    unsigned int i;
	    // print sequence as well as its length
	    /*
	    std::vector<std::pair<const char*, unsigned> >::iterator it;
	    for ( it=inputs.begin() ; it < inputs.end(); it++ ) {
		std::cout << "length = " << (*it).second << std::endl;
		for (i = 0; i < (*it).second; i++)
		{
		    putchar(toupper((*it).first[i]));
		}
		std::cout << std::endl;
	    }*/
	    

	    for (i=0; 9*i < tokens.size(); i++) {
		tokens[9*i+2].erase(tokens[9*i+2].end()-1); 
	    }
	    
	    // calculate the longest hishape and mfe
	    unsigned int length_longest_hishape = 0, length_hishape = 0, length_longest_mfe = 0, length_mfe = 0;
	    for (i=0; 9*i < tokens.size(); i++) {
		length_hishape = tokens[9*i+2].length();
		length_mfe = tokens[9*i+4].length();
		if (length_longest_hishape < length_hishape) {
		    length_longest_hishape = length_hishape;
		    length_longest_mfe = length_mfe;
		}
	    }
					
	    // prepare format_pattern and 2 maps to output
	    std::string format_pattern; 
	    std::map<std::string,std::string> hishape_mfehishape, mfehishape_pp;
	    i = 0;
	    for (i=0; 9*i < tokens.size(); i++) {
		format_pattern = ( boost::format("%%%d.2f %%%ds") % (length_longest_mfe+5) % (length_longest_hishape+5) ).str();
														      // #### wrap hishape with "[...]" ####
		std::string mfehishape = (boost::format(format_pattern) % ((boost::lexical_cast<int>(tokens[9*i+4]))/100.0f) % ("["+tokens[9*i+2]+"]")).str();
		std::string pp = tokens[9*i+7];                                       
		hishape_mfehishape[tokens[9*i+2]] = mfehishape;
		mfehishape_pp[mfehishape] = pp;
	    }
											
	    // to delete the duplicate hishape records, for example double [46.5], we need hishape_mfehishape
	    // mfehishape combine with mfehishape_pp to output records
	    std::map<std::string,std::string> mfehishape_hishape = flip_map(hishape_mfehishape);
	    std::map<std::string,std::string>::reverse_iterator rit;
	    size_t found;
	    for ( rit=mfehishape_hishape.rbegin(); rit != mfehishape_hishape.rend(); rit++ ) {
	      found = (*rit).first.find_first_not_of(" ");
		  if ( (*rit).first[found] == '-' ) {
		      std::cout << mfehishape_pp[(*rit).first] << (*rit).first << std::endl;
		      std::string hishape = mfehishape_hishape[(*rit).first];
		      // process hishape in type 3 in the keys of the map: (,27,38,) --> -1m,(,27,38,)
		      if (hishape_type==HISHAPE_H_PLUS) {
			  boost::replace_all(hishape, "(", "-1m,(");
		      }    
		      (*hishape_ppmfehishape)["["+hishape+"]"] = mfehishape_pp[(*rit).first] + (*rit).first;
		  }    
	    }
	    std::map<std::string,std::string>::iterator it2;
	    for ( it2=mfehishape_hishape.begin(); it2 != mfehishape_hishape.end(); it2++ ) {
	      found = (*it2).first.find_first_not_of(" ");
		  if ( (*it2).first[found] != '-' ) {
		      std::cout << mfehishape_pp[(*it2).first] << (*it2).first << std::endl;
		      std::string hishape = mfehishape_hishape[(*it2).first];
		      // process hishape in type 3 in the keys of the map: (,27,38,) --> -1m,(,27,38,)
		      if (hishape_type==HISHAPE_H_PLUS) {
			  boost::replace_all(hishape, "(", "-1m,(");
		      } 
		      (*hishape_ppmfehishape)["["+hishape+"]"] = mfehishape_pp[(*rit).first] + (*rit).first;
		  }    
	    }    
	    //return hishape_ppmfehishape;
	    //======for (unsigned int i = 0; i<opts.repeats; ++i) obj.print_backtrack(std::cout, res);======                      
	} catch (std::exception &e) {
	    std::cerr << "Exception: " << e.what() << '\n';
	    std::exit(1);
	}
	//return NULL;
}


void calculateHishapeProbabilities(const std::vector<std::pair<const char*, unsigned> > &inputs, uint32_t kbest, unsigned int hishape_type, const std::string &match_str, unsigned int exact, std::map<std::string,double> *hishape_pf)
{
  	    std::vector<std::string> tokens;
	    if (hishape_type==4) {
		
		gapc::hishapeh_pfx_cls obj;
		obj.init(inputs, kbest, match_str, exact);
		obj.cyk();
		gapc::hishapeh_pfx_ret res = obj.run();
		obj.print_result(std::cout, res);
		
		//======for (unsigned int i = 0; i<opts.repeats; ++i) obj.print_backtrack(std::cout, res);======
		List_Ref<std::pair<std::pair<Rope, std::pair<mfeanswer, pfanswer> > , intrusive_ptr<Backtrace<String, unsigned int> > > > btp;
		#ifdef WINDOW_MODE
		    btp  = obj.bt_proxy_nt_struct(obj.t_0_left_most, obj.t_0_right_most);
		#else
		    btp  = obj.bt_proxy_nt_struct();
		#endif
		intrusive_ptr<Backtrace<String, unsigned int> >  bt   = execute_backtrack_k(btp);
		intrusive_ptr<Backtrace_List<String, unsigned int> > l = boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
		assert(l);
		/// save the data in a vector 'tokens' /
		for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
		    std::stringstream ss;
		    (*i)->print(ss);

		    std::copy(std::istream_iterator<std::string>(ss),
			      std::istream_iterator<std::string>(),
			      std::back_inserter<std::vector<std::string> >(tokens) );
		}
	    } else if (hishape_type==3) {
		gapc::hishapehplus_pfx_cls obj;
		obj.init(inputs, kbest, match_str, exact);
		obj.cyk();
		gapc::hishapehplus_pfx_ret res = obj.run();
		obj.print_result(std::cout, res);   
		//======for (unsigned int i = 0; i<opts.repeats; ++i) obj.print_backtrack(std::cout, res);======
		List_Ref<std::pair<std::pair<Rope, std::pair<mfeanswer, pfanswer> > , intrusive_ptr<Backtrace<String, unsigned int> > > > btp;
		#ifdef WINDOW_MODE
		    btp  = obj.bt_proxy_nt_struct(obj.t_0_left_most, obj.t_0_right_most);
		#else
		    btp  = obj.bt_proxy_nt_struct();
		#endif
		intrusive_ptr<Backtrace<String, unsigned int> >  bt   = execute_backtrack_k(btp);
		intrusive_ptr<Backtrace_List<String, unsigned int> > l = boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
		assert(l);
		/// save the data in a vector 'tokens' /
		for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
		    std::stringstream ss;
		    (*i)->print(ss);

		    std::copy(std::istream_iterator<std::string>(ss),
			      std::istream_iterator<std::string>(),
			      std::back_inserter<std::vector<std::string> >(tokens) );
		}
	    } else if (hishape_type==2) {
		gapc::hishapem_pfx_cls obj;
		obj.init(inputs, kbest, match_str, exact);
		obj.cyk();
		gapc::hishapem_pfx_ret res = obj.run();
		obj.print_result(std::cout, res);  
		//======for (unsigned int i = 0; i<opts.repeats; ++i) obj.print_backtrack(std::cout, res);======
		List_Ref<std::pair<std::pair<Rope, std::pair<mfeanswer, pfanswer> > , intrusive_ptr<Backtrace<String, unsigned int> > > > btp;
		#ifdef WINDOW_MODE
		    btp  = obj.bt_proxy_nt_struct(obj.t_0_left_most, obj.t_0_right_most);
		#else
		    btp  = obj.bt_proxy_nt_struct();
		#endif
		intrusive_ptr<Backtrace<String, unsigned int> >  bt   = execute_backtrack_k(btp);
		intrusive_ptr<Backtrace_List<String, unsigned int> > l = boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
		assert(l);
		/// save the data in a vector 'tokens' /
		for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
		    std::stringstream ss;
		    (*i)->print(ss);

		    std::copy(std::istream_iterator<std::string>(ss),
			      std::istream_iterator<std::string>(),
			      std::back_inserter<std::vector<std::string> >(tokens) );
		} 
	    } else if (hishape_type==1) {
		gapc::hishapeb_pfx_cls obj;
		obj.init(inputs, kbest, match_str, exact);
		obj.cyk();
		gapc::hishapeb_pfx_ret res = obj.run();
		obj.print_result(std::cout, res); 
		//======for (unsigned int i = 0; i<opts.repeats; ++i) obj.print_backtrack(std::cout, res);======
		List_Ref<std::pair<std::pair<Rope, std::pair<mfeanswer, pfanswer> > , intrusive_ptr<Backtrace<String, unsigned int> > > > btp;
		#ifdef WINDOW_MODE
		    btp  = obj.bt_proxy_nt_struct(obj.t_0_left_most, obj.t_0_right_most);
		#else
		    btp  = obj.bt_proxy_nt_struct();
		#endif
		intrusive_ptr<Backtrace<String, unsigned int> >  bt   = execute_backtrack_k(btp);
		intrusive_ptr<Backtrace_List<String, unsigned int> > l = boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
		assert(l);
		/// save the data in a vector 'tokens' /
		for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
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
	    unsigned int i = 0;
	    for (i=0; 13*i < tokens.size(); i++) {
	        (*hishape_pf)[tokens[13*i+2]] = boost::lexical_cast<double>(tokens[13*i+7]);
	    }
}



// ostream* stmOutput = &std::cout;
// method for multialign calculation
void multiAlign(const std::vector<std::vector<std::string> > &allHishapes, size_t vecIndex, const std::string &hishapeSoFar, /*std::vector<CHelixIndexTreeNode::CTreeNodePointer> treeVec, hishape::comp::algorithms::CShapiroAligner tShapiroAligner,*/ float scoreSoFar, float ratio, float& scoreMin, std::string& hishapeMin/*, std::ostream& stmOutput*/)
{
    /*stmOutput << ".";*/
    if (vecIndex >= allHishapes.size())
    {
        //std::cout << "hishapeSoFar=" << hishapeSoFar << std::endl;
	//std::cout << "scoreSoFar=" << scoreSoFar << std::endl;
	if (scoreMin > scoreSoFar)
	{
	    scoreMin = scoreSoFar;
	    hishapeMin = hishapeSoFar;
	}
        return;
    }
    for (size_t i=0; i<allHishapes[vecIndex].size(); i++)
    {
        //std::cout << vecIndex << ":" << hishapeSoFar << ":" << allHishapes[vecIndex][i] << std::endl;

	//##CHelixIndexTreeNode::CTreeNodePointer ptTree2 = HishapeToTree(allHishapes[vecIndex][i]);
	//##tShapiroAligner.SetSecondTree(ptTree2);
	
	std::vector<std::string> hishapeSoFarVec;
	tokenize(hishapeSoFar, hishapeSoFarVec, "_");
	float scoreThisTurn = 0.0f;
	for (size_t j=0; j<hishapeSoFarVec.size(); j++)
	{
	    //std::cout << hishapeSoFarVec[j] << " vs " << allHishapes[vecIndex][i] << " = ";
	    
	    //##CHelixIndexTreeNode::CTreeNodePointer ptTree1 = HishapeToTree(hishapeSoFarVec[j]);    
	    //##tShapiroAligner.SetFirstTree(ptTree1);
	    //##tShapiroAligner.Align(false,false);
	    CHelixIndexTreeNode::CTreeNodePointer ptTree2 = HishapeToTree(allHishapes[vecIndex][i]);
	    CHelixIndexTreeNode::CTreeNodePointer ptTree1 = HishapeToTree(hishapeSoFarVec[j]); 
	    hishape::comp::algorithms::CShapiroAligner tShapiroAligner(ptTree1, ptTree2, ratio);
	    hishape::comp::CScoreScheme tScoreScheme;
	    tShapiroAligner.SetScoreScheme(tScoreScheme);
	    tShapiroAligner.Align(false,false);
            float scoreThisPair = tShapiroAligner.GetScore(); 
	    //std::cout << "(scoreThisPair=" << scoreThisPair;
	    scoreThisTurn += scoreThisPair; 
	    //std::cout << ",scoreThisTurn=" << scoreThisTurn << ")" << std::endl;
	    //##delete ptTree1;
	}
	//scoreSoFar += scoreThisTurn;
	if ((scoreSoFar+scoreThisTurn) < scoreMin)
	{
            multiAlign(allHishapes, vecIndex+1, hishapeSoFar+"_"+allHishapes[vecIndex][i], /*treeVec, tShapiroAligner,*/ scoreSoFar+scoreThisTurn, ratio, scoreMin, hishapeMin/*, stmOutput*/);
	}
	//##delete ptTree2;
	
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
    unsigned int convert = 0;
    unsigned int align = 0;
    unsigned int proba = 0;
    unsigned int rapid = 0;
    unsigned int partition = 0;
    unsigned int Partition = 0;
    float ratio = 5.0f;
    //unsigned int filter1 = 0;
    double thresh = 10000.0;
    unsigned int locmin = 0;
    float minh = 0.0f;
    double T = 37.0;
    std::string P = "";
    unsigned int multialign = 0;


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
	("file,f", po::value< std::vector<std::string> >(), "Read sequence(s) from specified file")
	("sequence,s", po::value< std::vector<std::string> >(), "Read sequence from specified sequence")
	("type,t", po::value<unsigned int>(&hishape_type)->default_value(2), "Specify hishape abstract type")
	//(4 stands for an algorithm in which helix indices "
	//"from hairpin loops are considered, 3 represents an algorithm in which helix indices from hairpin loops and nested information are considered, "
	//"2 stands for an algorithm in which helix indices from hairpin- and multiloops are in nested form, 1 represents an algorithm in which all helix indices are in nested form). "
	//"The letter 'm' is attached to the end of helix indices from multiloops in type 1 as well as in type 2, " 
	//"while the letters 'b' and 'i' denote helix indices from bulge- and internal loops in type 1, respectively.")
	("kbest,k", po::value<unsigned int>(&kbest)->default_value(3), "Choose the k-best classes") 
	("Related,R", po::value<std::string>(&match_str) /*->default_value("4,14.5")*/, "Output related hishapes regarding the given helix indices")
	("exact,e", po::value<unsigned int>(&exact)->implicit_value(1), "Faster computation of exact RNA hishapes")
        ("convert,c", po::value<unsigned int>(&convert)->implicit_value(1), "Convert secondary structure to hishape")
        ("align,a", po::value<unsigned int>(&align)->implicit_value(1), "Align two Hishapes by minimizing Tree editing distance (HiTed). Given a fas-format file, the option will compare the 1st record with the 2nd one, the 3rd with the 4th, and so on.")
        ("proba,p", po::value<unsigned int>(&proba)->implicit_value(1), "Calculate accumulated hishape probabilities")
	("rapid,q", po::value<unsigned int>(&rapid)->implicit_value(1), "Faster computation of exact RNA hishape probabilities")
	("partition,z", po::value<unsigned int>(&partition)->implicit_value(1), "Assign partition functions of hishape classes")
	("Partition,Z", po::value<unsigned int>(&Partition)->implicit_value(1), "Assign exact partition functions of hishape classes")
        ("ratio,r", po::value<float>(&ratio)->default_value(5.0f), "Specify the scaling coefficient of the score calculated from helix index distance "
        "=> Score_all=Score_conversion+ratio*Score_distance. " 
        "When ratio is equal to 0.0, it takes only the score from helix index type conversion.")  // into account
	//("filter1,i", po::value<unsigned int>(&filter1)->implicit_value(1), "A substructure A will not be added in the external loop or multiloop if the free energy of A > x kcal/mol. x value can be given with option thresh and its default value is 0")
        ("thresh,x", po::value<double>(&thresh)->implicit_value(10000.0), "Specify a threshold value for initstem")
        ("locmin,l", po::value<unsigned int>(&locmin)->implicit_value(1), "Print only locally minimal hishreps")
	("minh,m", po::value<float>(&minh)->default_value(0.0f), "Print only minima with energy barrier greater than the specified delta value")
        ("Temperature,T", po::value<double>(&T)->implicit_value(37.0), "Specify the temperature for energy parameters in degree Celsius.") 
	("Parameter,P", po::value<std::string>(&P), "Specify the path of energy parameter file.")
        ("multialign,u", po::value<unsigned int>(&multialign)->implicit_value(1), "Calculate the minimum Sum-Of-Pairs (SP) score of a multiple hishape alignment.") 
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
	conflicting_options(vm, "align", "convert");
	//conflicting_options(vm, "kbest", "proba");
	conflicting_options(vm, "proba", "rapid");
	conflicting_options(vm, "proba", "partition");
	conflicting_options(vm, "proba", "Partition");
	conflicting_options(vm, "rapid", "partition");
	conflicting_options(vm, "rapid", "Partition");
	conflicting_options(vm, "partition", "Partition");
	//option_dependency(vm, "ratio", "align");
	//option_dependency(vm, "ratio", "multialign");
	//option_dependency(vm, "locmin", "threshold");
	option_dependency(vm, "exact", "Related");
	//option_dependency(vm, "minh", "locmin");
	

	
        //std::cout << "vm.size()=" << vm.size() << std::endl;
        if (vm.count("help") || vm.size()==4) {
            std::cout << "Usage: RNAHeliCes (-s INPUT|-f INPUT-file) [options]\n";  // [<sequence>|<first hishape>] [<second hishape>]
	    std::cout << "Examples:" << std::endl;
	    std::cout << "  ./RNAHeliCes ../examples/riboswitches.fas -c" << std::endl;            // convert every record 
	    std::cout << "  ./RNAHeliCes '[27]' '[36.5m,(,27,38,)]' -a -r 1" << std::endl;
	    std::cout << "  ./RNAHeliCes ../examples/riboswitches.fas -a -t 1 -r 1" << std::endl;  // compare the 1st and 2nd, 3rd and 4nd
	    std::cout << "  ./RNAHeliCes ../examples/ires_picornaviridae.fa -u -k 2 -t 2 -r 2" << std::endl;
	    std::cout << "  ./RNAHeliCes -s GGGGGGCCCCCC" << std::endl;
	    std::cout << "  ./RNAHeliCes ../examples/collosoma_slrna.seq -P ./librna/vienna/rna_turner2004.par" << std::endl;
	    std::cout << "  ./RNAHeliCes ../examples/collosoma_slrna.seq -R 27,38 -k 20" << std::endl;
	    //std::cout << "  ./RNAHeliCes -f ../examples/multiple_sequences.fa -u -k 2 -t 2 -r 2" << std::endl;
	    std::cout << "  ./RNAHeliCes -f ../examples/collosoma_slrna.seq -k 10 -R 36.5m,41.5,27 -e -t 2 -x0" << std::endl;
	    std::cout << "  ./RNAHeliCes -f ../examples/collosoma_slrna.seq -k 100 -x0 -l" << std::endl;
            std::cout << desc; //<< '\n';
            //std::cout << "input files are: " << vm["file"].as< std::vector<std::string> >() << '\n';
            //std::cout << "input sequences are: " << vm["sequence"].as< std::vector<std::string> >() << '\n';
            //std::cout << "positional options are: " << vm["sequence"].as< std::vector<std::string> >() << '\n';  // if it is empty, return bad_any_cast error
            return 0;
        }

        if (vm.count("version")) {
            std::cout << "RNAHeliCes 2.0.12 (Feb. 09, 2012)\n";  // TODO: use the C++ buildin date function 
            return 0;
        }
        
//**        if (vm.count("locmin")) {
//**	 	if (thresh > 0) 
//**		{
//**		    std::cout << "The energy threshold should be less or equal to zero if locmin-filter is used!" << std::endl;
//**		    return 1;
//**		} 
//**	}
        
        
        //std::string strErrors;
        if (vm.count("file"))
        {
           // std::string fileinput_str = (vm["file"].as< std::vector<std::string> >())[0];
           // char * fileinput_char = new char[fileinput_str.length()+1];
           // strcpy(fileinput_char, fileinput_str.c_str());


            

	    std::vector<std::string> files = vm["file"].as< std::vector<std::string> >();
	    //std::vector<std::string>::iterator it;
            //for ( it=files.begin() ; it < files.end(); it++ )
            //    std::cout << "input file names: " << *it;
            //std::cout << std::endl;
            if (!files.empty()) {
	        // ##################
	        // ## convert mode ##
	        if (convert==1 && files.size()==1) 
		{
		  	if (hishape_type < 1 || hishape_type > 4) 
			{
				std::cout << "The hishape abstract type can only be specified between 1 and 4." << std::endl;  // For multiple hishape alignment, 
				return 1;
			}

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
					std::string ss = "";
					int line_no = 0;
					while ( getline(fiInputFile, line) )
					{

						if ( (line_no%3 == 0) ) {
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
						    if (line_no%3 == 1)
							seq = line;
						    if (line_no%3 == 2)
							ss = line;
						}


						if (header != "" && seq != "" && ss != "") {
							std::cout << header << std::endl;
							std::cout << seq << std::endl;
							std::cout << ss << std::endl;
                                                        std::cout << sSToHishape(header+"\n"+seq+"\n"+ss+"\n", hishape_type) /*<< " (hishape type = " << hishape_type << ")" */<< std::endl;
							header = "";
							seq = "";
							ss = "";
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
	        // ##################
	        // #### align mode ##
		if (align==1 && files.size()==2) 
		{
                        
		        std::string strInput1,strInput2;
		        if (!files[0].empty() && files[0].at(0) == '[' &&  !files[1].empty() && files[1].at(0) == '[')
			{
				strInput1 = files[0];
				strInput2 = files[1];
				pairAlign(strInput1, strInput2, hishape_type, ratio);
			}
		}
		else
		if (align==1 && files.size()==1) 
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

					std::string header1 = "";
					std::string seq1 = "";
					std::string ss1 = "";
					std::string header2 = "";
					std::string seq2 = "";
					std::string ss2 = "";
					int line_no = 0;
					while ( getline(fiInputFile, line) )
					{

						if ( (line_no%6 == 0) ) {
						    if (line[0] == '>')
						    {
							// get header only
							//header1 = line.substr(1);
							header1 = line;
						    }
						    else
						    {
							throw "ERROR: wrong faa-format in input file (\"" + files[0] + "\")!";
						    }
						}
						else if ( (line_no%6 == 3) ) {
						    if (line[0] == '>')
						    {
							// get header only
							//header2 = line.substr(1);
							header2 = line;
						    }
						    else
						    {
							throw "ERROR: wrong faa-format in input file (\"" + files[0] + "\")!";
						    }
						}
						else 
						{
						    if (line_no%6 == 1)
							seq1 = line;
						    if (line_no%6 == 2)
							ss1 = line;
						    if (line_no%6 == 4)
							seq2 = line;
						    if (line_no%6 == 5)
							ss2 = line;
						}


						if (header1 != "" && header2 != "" && seq1 != "" && seq2 != "" && ss1 != "" && ss2 != "") {
						        std::cout << "#### " << header1 + " vs. " + header2 + " ####" << std::endl;  
                                                        pairAlign(sSToHishape(header1+"\n"+seq1+"\n"+ss1+"\n", hishape_type), sSToHishape(header2+"\n"+seq2+"\n"+ss2+"\n", hishape_type), hishape_type, ratio);
							header1 = "";
							header2 = "";
							seq1 = "";
							seq2 = "";
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
		// ###############################
	        // #### multiple alignment mode ##
		if (multialign==1 && files.size()==1) 
		{
			//std::list< std::pair<std::pair<std::string, std::string>, std::vector<std::string> >  > hishapesOfMultiSequences;
			std::string comment,sequence;
			//unsigned int rowNumber;
			//unsigned int i,j,m,n;


			std::fstream fiInputFile1;
			if (false == fiInputFile1.is_open())
			{       
				try
				{
					// get first input file
					fiInputFile1.open(files[0].c_str(), std::ios::in);
					// checks if the file was opened
					if ((false == fiInputFile1) || (false == fiInputFile1.is_open()))
					{
						// not opened, throws an exception
						throw "Error loading first input file (\"" + files[0] + "\")!";
					}
				}
				catch (...)
				{
					std::cout << "ERROR: unable to open first input file (\"" + files[0] + "\")!\n";
				}
			}       
			if (true == fiInputFile1.is_open())
			{
				try
				{       std::string comment,sequence;
				        std::vector<std::string> allComments;
				        std::vector<std::vector<std::string> > allHishapes;
					std::vector<std::map<std::string, std::string> > allHishapePpmfehishape;

					if (hishape_type < 1 || hishape_type > 4) 
					{
					    std::cout << "The hishape abstract type can only be specified between 1 and 4." << std::endl;  // For multiple hishape alignment, 
					    return 1;
					}
					while (fiInputFile1.eof() == false)  // fiInputFile1.good()
					{    
					      getline (fiInputFile1,comment);
					      getline (fiInputFile1,sequence);
					      if (!(comment.empty()) && !(sequence.empty())) 
					      {
						    allComments.push_back(comment);
						    std::cout << sequence << std::endl;
						    inputs.clear();
						    inputs.push_back(std::make_pair(sequence.c_str(), (unsigned int)sequence.size()));
						    std::map<std::string, std::string> hishape_ppmfehishape;
						    calculateHishapes(inputs,kbest,hishape_type, &hishape_ppmfehishape, exact, thresh, locmin);

                                                    std::vector<std::string> keys;
						    // get the keys (and not the values) from the map 'hishape_ppmfehishape' and save them into the vector 'keys': boost::transform_iterator have the same function 
						    std::transform(hishape_ppmfehishape.begin(), hishape_ppmfehishape.end(), std::inserter(keys, keys.begin()), GetKey<std::map<std::string, std::string>::value_type>());
						    //std::copy(s.begin(), s.end(), std::ostream_iterator<std::string>(std::cout, " "));

						    allHishapes.push_back(keys);
						    allHishapePpmfehishape.push_back(hishape_ppmfehishape);
					      }
					}



					//**CHelixIndexTreeNode::CTreeNodePointer ptTree1 = HishapeToTree("[]");
					//**CHelixIndexTreeNode::CTreeNodePointer ptTree2 = HishapeToTree("[]");
					
					if (ratio < 0.0f || ratio > 10000.0f) 
					{
					    std::cout << "The ratio value should be specified between 0.0 and 10000.0!" << std::endl;
					    return 1;
					}
					//**hishape::comp::algorithms::CShapiroAligner tShapiroAligner(ptTree1, ptTree2, ratio);
					//**hishape::comp::CScoreScheme tScoreScheme;
					//**tShapiroAligner.SetScoreScheme(tScoreScheme);
					//**tShapiroAligner.Align(false,false);
					float scoreMin = std::numeric_limits<float>::infinity();
					std::string hishapeMin = "";
					multiAlign(allHishapes, 0, "", /*tShapiroAligner,*/ 0.0f, ratio, scoreMin, hishapeMin/*, std::cout*/);
					std::cout << "The minimum Sum-Of-Pairs (SP) score = "<< scoreMin << std::endl;
					//std::cout << hishapeMin << std::endl;
					std::vector<std::string> tokens;
					tokenize(hishapeMin,tokens,"_");
					
					for (size_t i=0; i<allHishapePpmfehishape.size(); i++)
					{
					    std::cout << allComments[i] << std::endl;  // output comment
					    
					    std::string foundPpmfehishape = allHishapePpmfehishape[i].find(tokens[i])->second;
					    std::vector<std::string> tokens2;
					    tokenize(foundPpmfehishape,tokens2," "); 
                                            for (size_t j=0; j<tokens2.size(); j++)
					    {
					        std::cout << tokens2[j] << "    ";
					    }
					    std::cout << std::endl;
					}
				}
				catch (std::exception &ex)
				{
					std::cout << "ERROR: " << ex.what() << std::endl;
				}
				catch (...)
				{
					std::cout << "ERROR: unexpected exception!" << std::endl;
				}
			}
			if (true == fiInputFile1.is_open())
			{
				fiInputFile1.close();
			}
		}
		else  
	        // ###################
	        // ## standard mode ##
		if (convert==0 && align==0 && multialign==0 && files.size()==1)
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
		    std::cout << "Status: align=" << align << ",convert=" << convert << ",multialign=" << multialign << ",files.size()=" << files.size() << std::endl;
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

//     TE.init(LEfile, SEfile, TerminalSMHEfile, TerminalSMIEfile, SingleMIEfile, 
//     OneTwoMIEfile, TetraLoopEfile, TwoTwoMIEfile, d3Efile, d5Efile);
    TE.init_SE(SEfile);

    double pfunc = 0.0;
    if (proba==1 || rapid==1)
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



    if (convert==0 && align==0 && multialign==0 && !inputs.empty())
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
	    if ((proba == 0 && rapid == 0 && partition == 0 && Partition == 0) || rapid==1 || Partition==1) {  // proba == 0 || 
		if (hishape_type==4) {
		    gapc::hishapeh_pp_cls obj;
		    
		    obj.init(inputs, /*(int)kbest*1.1*/kbest, match_str, exact, thresh*100, locmin, true);
		    
                    // placeholder
		    obj.cyk();  // do nothing
		    gapc::hishapeh_pp_ret res = obj.run();
		    //## gapc::add_event("end_computation");
		    
		    obj.print_result(std::cout, res);  // do nothing
		    //##gapc::add_event("end_result_pp");
		    
		    // extends the kbest to 10% 
		    //res.set_k(kbest);

		    //List_Ref<std::pair<std::pair<Rope, mfeanswer> , intrusive_ptr<Backtrace<String, unsigned int> > > > btp  = obj.bt_proxy_nt_struct();
		    //intrusive_ptr<Backtrace<String, unsigned int> >  bt = execute_backtrack_k(btp);
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
		    obj.init(inputs, kbest, match_str, exact, thresh*100, locmin, true);
		    obj.cyk();
		    gapc::hishapehplus_pp_ret res = obj.run();
		    obj.print_result(std::cout, res);   
		    //List_Ref<std::pair<std::pair<Rope, mfeanswer> , intrusive_ptr<Backtrace<String, unsigned int> > > > btp  = obj.bt_proxy_nt_struct();
		    //intrusive_ptr<Backtrace<String, unsigned int> >  bt = execute_backtrack_k(btp);
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
		    obj.init(inputs, kbest, match_str, exact, thresh*100, locmin, true);
		    obj.cyk();
		    gapc::hishapem_pp_ret res = obj.run();
		    obj.print_result(std::cout, res);  
		    //List_Ref<std::pair<std::pair<Rope, mfeanswer> , intrusive_ptr<Backtrace<String, unsigned int> > > > btp  = obj.bt_proxy_nt_struct();
		    //intrusive_ptr<Backtrace<String, unsigned int> >  bt = execute_backtrack_k(btp);
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
		    obj.init(inputs, kbest, match_str, exact, thresh*100, locmin, true);
		    obj.cyk();
		    gapc::hishapeb_pp_ret res = obj.run();
		    obj.print_result(std::cout, res); 
		    //List_Ref<std::pair<std::pair<Rope, mfeanswer> , intrusive_ptr<Backtrace<String, unsigned int> > > > btp  = obj.bt_proxy_nt_struct();
		    //intrusive_ptr<Backtrace<String, unsigned int> >  bt = execute_backtrack_k(btp);
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
	    }
	    else if (proba==1 || partition == 1) 
	    {
		if (hishape_type==4) {
		    
		    gapc::hishapeh_pfx_cls obj;
		    //TODO: implement pfx with obj.init(inputs, kbest, match_str, exact, thresh*100, locmin, true);
		    obj.init(inputs, kbest, match_str, exact);
		    obj.cyk();
		    gapc::hishapeh_pfx_ret res = obj.run();
		    obj.print_result(std::cout, res);
		    
		    //======for (unsigned int i = 0; i<opts.repeats; ++i) obj.print_backtrack(std::cout, res);======
		    List_Ref<std::pair<std::pair<Rope, std::pair<mfeanswer, pfanswer> > , intrusive_ptr<Backtrace<String, unsigned int> > > > btp;
		    #ifdef WINDOW_MODE
			btp  = obj.bt_proxy_nt_struct(obj.t_0_left_most, obj.t_0_right_most);
		    #else
			btp  = obj.bt_proxy_nt_struct();
		    #endif
		    intrusive_ptr<Backtrace<String, unsigned int> >  bt   = execute_backtrack_k(btp);
		    intrusive_ptr<Backtrace_List<String, unsigned int> > l = boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
		    assert(l);
		    /// save the data in a vector 'tokens' /
		    for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
			std::stringstream ss;
			(*i)->print(ss);

			std::copy(std::istream_iterator<std::string>(ss),
				  std::istream_iterator<std::string>(),
				  std::back_inserter<std::vector<std::string> >(tokens) );
		    }
		} else if (hishape_type==3) {
		    gapc::hishapehplus_pfx_cls obj;
		    obj.init(inputs, kbest, match_str, exact);
		    obj.cyk();
		    gapc::hishapehplus_pfx_ret res = obj.run();
		    obj.print_result(std::cout, res);   
		    //======for (unsigned int i = 0; i<opts.repeats; ++i) obj.print_backtrack(std::cout, res);======
		    List_Ref<std::pair<std::pair<Rope, std::pair<mfeanswer, pfanswer> > , intrusive_ptr<Backtrace<String, unsigned int> > > > btp;
		    #ifdef WINDOW_MODE
			btp  = obj.bt_proxy_nt_struct(obj.t_0_left_most, obj.t_0_right_most);
		    #else
			btp  = obj.bt_proxy_nt_struct();
		    #endif
		    intrusive_ptr<Backtrace<String, unsigned int> >  bt   = execute_backtrack_k(btp);
		    intrusive_ptr<Backtrace_List<String, unsigned int> > l = boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
		    assert(l);
		    /// save the data in a vector 'tokens' /
		    for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
			std::stringstream ss;
			(*i)->print(ss);

			std::copy(std::istream_iterator<std::string>(ss),
				  std::istream_iterator<std::string>(),
				  std::back_inserter<std::vector<std::string> >(tokens) );
		    }
		} else if (hishape_type==2) {
		    gapc::hishapem_pfx_cls obj;
		    obj.init(inputs, kbest, match_str, exact);
		    obj.cyk();
		    gapc::hishapem_pfx_ret res = obj.run();
		    obj.print_result(std::cout, res);  
		    //======for (unsigned int i = 0; i<opts.repeats; ++i) obj.print_backtrack(std::cout, res);======
		    List_Ref<std::pair<std::pair<Rope, std::pair<mfeanswer, pfanswer> > , intrusive_ptr<Backtrace<String, unsigned int> > > > btp;
		    #ifdef WINDOW_MODE
			btp  = obj.bt_proxy_nt_struct(obj.t_0_left_most, obj.t_0_right_most);
		    #else
			btp  = obj.bt_proxy_nt_struct();
		    #endif
		    intrusive_ptr<Backtrace<String, unsigned int> >  bt   = execute_backtrack_k(btp);
		    intrusive_ptr<Backtrace_List<String, unsigned int> > l = boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
		    assert(l);
		    /// save the data in a vector 'tokens' /
		    for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
			std::stringstream ss;
			(*i)->print(ss);

			std::copy(std::istream_iterator<std::string>(ss),
				  std::istream_iterator<std::string>(),
				  std::back_inserter<std::vector<std::string> >(tokens) );
		    } 
		} else if (hishape_type==1) {
		    gapc::hishapeb_pfx_cls obj;
		    obj.init(inputs, kbest, match_str, exact);
		    obj.cyk();
		    gapc::hishapeb_pfx_ret res = obj.run();
		    obj.print_result(std::cout, res); 
		    //======for (unsigned int i = 0; i<opts.repeats; ++i) obj.print_backtrack(std::cout, res);======
		    List_Ref<std::pair<std::pair<Rope, std::pair<mfeanswer, pfanswer> > , intrusive_ptr<Backtrace<String, unsigned int> > > > btp;
		    #ifdef WINDOW_MODE
			btp  = obj.bt_proxy_nt_struct(obj.t_0_left_most, obj.t_0_right_most);
		    #else
			btp  = obj.bt_proxy_nt_struct();
		    #endif
		    intrusive_ptr<Backtrace<String, unsigned int> >  bt   = execute_backtrack_k(btp);
		    intrusive_ptr<Backtrace_List<String, unsigned int> > l = boost::dynamic_pointer_cast<Backtrace_List<String, unsigned int> >(bt);
		    assert(l);
		    /// save the data in a vector 'tokens' /
		    for (Backtrace_List<String, unsigned int>::iterator i = l->begin(); i != l->end(); ++i) {
			std::stringstream ss;
			(*i)->print(ss);

			std::copy(std::istream_iterator<std::string>(ss),
				  std::istream_iterator<std::string>(),
				  std::back_inserter<std::vector<std::string> >(tokens) );
		    }
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
	    

	    if (proba == 0 && rapid == 0 && partition == 0 && Partition == 0) 
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


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////// post processing of stable local minima by checking barrier energies pairwise ///////////////

		time_t initial = time(NULL), first = time(NULL), second = time(NULL);
		if(minh>0) {
		        //cf.setConfig(configfile);
			//std::cout << "TwoTwoMIEfile=" << TwoTwoMIEfile << std::endl;

			//std::cout << "STACKING=" << TE.STACKING << std::endl;
		        std::vector<std::vector<double> > vvd;
			
		        std::string seq(inputs[0].first, inputs[0].second);
                        for(unsigned int i=0;i<seq.length();i++)
			{
			    seq[i] = toupper(seq[i]);
			}
			
			std::vector<RNA> vrna;
			std::set<stem> stem_set;
			for (i=0; 9*i < tokens.size(); i++) {  
			    vrna.push_back(RNA(tokens[9*i+7], boost::lexical_cast<int>(tokens[9*i+4]), "["+tokens[9*i+2]+"]"));//, "["+tokens[9*i+2]+"]"));
                            RStack rs_hishrep(tokens[9*i+7]);
			    FOR(t,0,sz(rs_hishrep.vstems)) {
			      stem_set.insert(rs_hishrep.vstems[t]);
// 			      if (!st.containST(rs_hishrep.vstems[t])) {
// 				st.add(rs_hishrep.vstems[t]);
// 			      }
			    }
			}
			stemlist st(seq, stem_set);  // st refers to stemlist
			

			first = second ;
	//		gettimeofday(&ms_first, NULL);
	//		cerr<<"The elapsed time is "<<(ms_first.tv_sec- ms_initial.tv_sec)<<" seconds "<<(ms_first.tv_usec- ms_initial.tv_usec)/1000<<" (miliseconds) "<<std::endl;


// 			std::cout << "te:";
// 			for (std::map<std::string, int>::iterator it = TE.STACKING.begin() ; it != TE.STACKING.end(); ++it)
// 			  std::cout << (*it).first << "," << (*it).second;
// 			std::cout << "========\n";
			
			//compute free energy of all possible stems in st.vstems.
			std::map<stem, int> stemEmap = computeStemE(TE, seq, st.vstems); 

	
// 			std::map<stem,int>::iterator it;
//                         for ( it=stemEmap.begin() ; it != stemEmap.end(); it++ ) {
// 			  std::cout<<"("<<(*it).first.outer.lp<<","<<(*it).first.outer.rp<<") ~ ("<<(*it).first.inner.lp<<","<<(*it).first.inner.rp<<")" <<" => "<< (*it).second << std::endl;
// 			}
		
			//Heuristic
//TREE			std::vector<TNode> tree;
			std::vector<int> vrna_cluster; //for saving indices of structures considered. 
			std::map<stemstem, std::vector<int> > CDSTmap;
			FOR(t,0,sz(st.vstems)) {
				FOR(k,0,sz(st.vstems)) {
					if(st.vstems[t].conflictwith(st.vstems[k]) && st.vstems[t]!=st.vstems[k]) {
					        //st.vstems[t].show();
						get_conflictingdegree_between_stems(CDSTmap, st.vstems[t], st.vstems[k]);
					}
				}
			}
//DEBUG std::cout << "minh=" << minh << std::endl;
//DEBUG std::cout << "cf.BARRIERCUTOFF=" << cf.BARRIERCUTOFF << std::endl;
//DEBUG std::cout << "is_barriercutoff=" << is_barriercutoff << std::endl;
//cf.BARRIERCUTOFF = 0;
			//if(is_barriercutoff) {
//DEBUG			        std::cout << barriercutoff << std::endl;
				cf.BARRIERCUTOFF = minh;
			//}		
			


			
			FOR(t,0,sz(vrna)) {
//DEBUG
// 				if(t%500==0 && t!=0) {
// 					std::cout<<"The algorithm has processed ["<<t<<"/"<<sz(vrna)<<"] local optimal stack configurations."<<std::endl;
// 				}
				RStack rs_t(vrna[t].str);
				//compute the distance between vrna[t] and each structure in vrna_cluster, 
				//if vrna[t] is close to vrna[k=vrna_cluster[x]], then do nothing
				//otherwise, add t to vrna_cluster.
				//assume the minimum conflicting degree is greater than conflicting degree cut off.
				double min_e = cf.BARRIERCUTOFF+1;//min_e could be any value greater than cf.BARRIERCUTOFF
				std::vector<double> reverse_vd;
				int indx = -1;
				for(int x = sz(vrna_cluster)-1; x>=0; x--) {
					int k = vrna_cluster[x];  // record k is last record of already stored ones, t is comming one
					RStack rs_k(vrna[k].str);
					double cdsum_tk = - compute_pairwise_conflictingdegree(TE, seq, rs_t, rs_k, stemEmap, CDSTmap);
					double cdsum_kt = - compute_pairwise_conflictingdegree(TE, seq, rs_k, rs_t, stemEmap, CDSTmap);
					double cdsum = min(cdsum_tk, cdsum_kt);
					reverse_vd.PB(cdsum);
					if(min_e > cdsum) {
						min_e = cdsum;
						indx = k;
					}
//TREE 					if(min_e < cf.BARRIERCUTOFF) {  // cut前留影
// 						//the conflicting degree between vrna[t] and vrna[indx] is less than the cut off, then node[t] is merged to node[indx]
// 						tree.PB(indx);
// 						break;
// 					} else {
// 					//	std::cout<<"Initially, the conflicting degree of vrna["<<k<<"] and vrna["<<t<<"] is "<<cdsum<<std::endl;
// 					//	std::cout<<vrna[k].str<<std::endl;
// 					//	std::cout<<vrna[t].str<<std::endl;
// 					}
				}
				
				if(min_e >= cf.BARRIERCUTOFF) {
					//the least conflicting degree between vrna[t] and any structure with indices vrna_cluster[.] is greater than the cut off 
					//t should be preserved.
					vrna_cluster.PB(t);
					//reverse_vd[x] contains the conflicting degree between vrna[t] and vrna[vrna_cluster[sz(reverse_vd)-1-x]]
					std::vector<double> vd;
					for(int x=sz(reverse_vd)-1; x>=0; x--) {
						vd.PB(reverse_vd[x]);
					}
					vvd.PB(vd);
//TREE					tree.PB(-1);//vrna[t] is not merged with any structure
	//				cerr<<"t="<<t<<std::endl;
	//				cerr<<"sz(vrna_cluster)="<<sz(vrna_cluster)<<", vrna_cluster.end() = "<<vrna_cluster[sz(vrna_cluster)-1]<<std::endl;
	//				cerr<<"sz(vd)="<<sz(vd)<<", sz(vvd)="<<sz(vvd)<<std::endl;
					
					assert(sz(vrna_cluster) == sz(vd)+1 && sz(vvd)==sz(vd)+1);
					FOR(x,0,sz(vd)) {
						int k = vrna_cluster[x];
					//	std::cout<<"Then, the conflicting degree of vrna["<<k<<"] and vrna["<<t<<"] is "<<vd[x]<<std::endl;
					//	std::cout<<vrna[k].str<<std::endl;
					//	std::cout<<vrna[t].str<<std::endl;
					}
				}
				
// 				if(sz(vrna_cluster) > cf.MAXMATRIXSZ) {
// 					cf.increaseBARRIERCUTOFF();
// 					std::cout<<"The number of stable local optimal stack configurations generated is "<<sz(vrna_cluster)<<", while we can process no more than MAXMATRIXSZ("<<cf.MAXMATRIXSZ<<")."<<std::endl;
// 					std::cout<<"We will increase the energy barrier cutoff to "<<cf.BARRIERCUTOFF<<" kcal/mol in order to decrease the number of stable local optimal structures."<<std::endl;
// 					applycutoff(vvd, cf.BARRIERCUTOFF, tree, vrna_cluster);
// 				}
			}
			
		
			std::vector<RNA> vrna_tmp;
			FOR(t,0,sz(vrna_cluster)) {
				vrna_tmp.PB(vrna[vrna_cluster[t]]);
				printf("%6d    ",t+1);
				std::cout<<vrna[vrna_cluster[t]].str<<"\t";
				printf("%6.2f",vrna[vrna_cluster[t]].e/100.0f);
				std::cout<<"\t"<<vrna[vrna_cluster[t]].hishape<<std::endl;
			}
			
//DEBUG
// 			std::cout<<"The algorithm has obtained "<<sz(vrna_cluster)<<" stable local optimal stack configurations from "<<sz(vrna)<<" local optimal stack configurations, and ranked them according to the minimal energy barrier associated."<<std::endl;
// 			std::cout<<"Note that numbers in the first column indicate ranks of structures according to their free energies."<<std::endl;
// 			std::cout<<"The algorithm has obtained "<<sz(vrna_cluster)<<" stable local optimal stack configurations from "<<sz(vrna)<<" local optimal stack configurations."<<std::endl;
// 			second = time(NULL);
// 			std::cout<<"The algorithm has spent "<<(second - first) << " seconds for computing pairwise energy barriers among stable local optimal stack configurations."<<std::endl;
// 			first = second ;

			print_phylipltm(vvd);

			mytable mt (vvd);
			std::vector<int> vorder;
			std::vector<int> vi = mt.DJ(vrna_cluster.size()/*knum*/, vrna_tmp, std::cout, vorder);
//DEBUG
// 			mt.show();
// 			
// 			second = time(NULL);
// 			std::cout<<"The algorithm has spent "<<(second - first) << " seconds for neighbor joining."<<std::endl;
// 			first = second ;

//DEBUG
// 			int t = sz(vorder)-1;
// 			while(t>=0) {
// 			        //TODO: potential errors with vorder
// 				std::cout<<vorder[t]<<"\t"<<vrna_tmp[vorder[t]].str<<"\t";
// 				printf("%6.2f",vrna_tmp[vorder[t]].e/100.0f);
// 				std::cout<<"\t"<<vrna_tmp[vorder[t]].hishape<<std::endl;
// 				t--;
// 			}
//DEBUG			
// 			std::cout<<std::endl<<std::endl<<"The top "<<sz(vi)<<" stable local optimal stack configurations ranked according to their associated minimal energy barrier are as follows:"<<std::endl;
// 			FOR(t,0,sz(vi)) {
// 				std::cout<<vi[t]+1<<"\t"<<vrna_tmp[vi[t]].str<<"\t";
// 				printf("%6.2f",vrna_tmp[vi[t]].e/100.0f);
// 				std::cout<<"\t"<<vrna_tmp[vi[t]].hishape<<std::endl;
// 			}    
//
// 			second = time(NULL);
// 			std::cout<<"The algorithm has spent "<<(second - first) << " seconds for writing the top "<<sz(vi)<<" stable local optimal stack configurations according to their associated minimal energy barrier to outfile."<<std::endl;
		}
		//std::cout<<"The algorithm for filtering out stable local minima and clustering them has spent "<<(second - initial)<< " seconds in total."<<std::endl;

		//	gettimeofday(&ms_second, NULL);
		//	cerr<<"The elapsed time is "<<(ms_second.tv_sec- ms_first.tv_sec)<<" seconds "<<(ms_second.tv_usec- ms_first.tv_usec)/1000<<" (miliseconds) "<<std::endl;
		//	cerr<<"The overall elapsed time is "<<((ms_second.tv_sec- ms_initial.tv_sec)*1000 + (ms_second.tv_usec- ms_initial.tv_usec)/1000)<<" (miliseconds) "<<std::endl;
		else {

		    // prepare format_pattern and output the results
		    std::string format_pattern = ( boost::format("%%%d.2f %%%ds") % (length_longest_mfe+5) % (length_longest_hishape+5) ).str();
		    i = 0;
		    for (i=0; 9*i < tokens.size(); i++) {  
															    // #### wrap hishape with "[...]" ####
			std::cout << tokens[9*i+7] << (boost::format(format_pattern) % ((boost::lexical_cast<int>(tokens[9*i+4]))/100.0f) % ("["+tokens[9*i+2]+"]")).str() << std::endl;
		    } 
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
	    else if (rapid == 1)
	    {
		// calculate the longest hishape and mfe
		unsigned int length_longest_hishape = 0, length_hishape = 0, length_longest_mfe = 0, length_mfe = 0;
		for (i=0; 9*i < tokens.size(); i++) {  
		    length_hishape = tokens[9*i+2].length();
		    length_mfe = tokens[9*i+4].length();
		    if (length_longest_hishape < length_hishape) {
			length_longest_hishape = length_hishape;
		    }
		    if (length_longest_mfe < length_mfe) {
			length_longest_mfe = length_mfe;
		    }
		}
		
		
		std::map<std::string, double> hishape_pf;
		std::string token_9_2 = "";
		// prepare format_pattern and output the results
		std::string format_pattern = ( boost::format("%%%d.2f %%%ds %%11.6f") % (length_longest_mfe+5) % (length_longest_hishape+4) ).str();
		i = 0;
		for (i=0; 9*i < tokens.size(); i++) {   
                    //std::string hishape = tokens[9*i+2];
                    double pf = 0.0;
		    if (hishape_pf.find(tokens[9*i+2]) != hishape_pf.end()) {
			pf = hishape_pf[tokens[9*i+2]];  
		    } else {
                        //match_str = tokens[9*i+2]
                        std::string match_str_for_function = "";    
			std::vector<std::string> helixIndices;
			tokenize(tokens[9*i+2], helixIndices);
			std::vector<std::string>::iterator it;
			for ( it=helixIndices.begin() ; it != helixIndices.end(); it++ )
			{
			    std::string currentHi = *it;
			    if ( (currentHi[currentHi.size()-1] != 'b') && (currentHi[currentHi.size()-1] != 'i') && \
				  (currentHi[currentHi.size()-1] != 'm') && (currentHi != "(") && (currentHi != ")") )  // in case h or _
			    {
				match_str_for_function = match_str_for_function + "," + currentHi;
			    }
			}

		        calculateHishapeProbabilities(inputs,1000000000,hishape_type, match_str_for_function, true, &hishape_pf);
			pf = hishape_pf[tokens[9*i+2]]; 
			//std::cout << "calculateHishapeProbabilities("<<","<<1000000000<<","<<hishape_type<<","<<match_str_for_function<<")" << std::endl;
		    }
		    token_9_2 = tokens[9*i+2];
		    if (token_9_2 != "_")
		        token_9_2.erase (token_9_2.end()-1); 
		    std::cout << tokens[9*i+7] << (boost::format(format_pattern) % ((boost::lexical_cast<int>(tokens[9*i+4]))/100.0f) % ("["+token_9_2+"]") % (pf/pfunc)).str() << std::endl;
		}
	    }
	    // output free energies of hishapes
	    else if (Partition == 1)
	    {
	        //_kT = -0.00198717*(273.15 + T);
		//std::cout << "T=" << T << ",_kT=" << _kT << std::endl;
		// calculate the longest hishape and mfe
		unsigned int length_longest_hishape = 0, length_hishape = 0, length_longest_mfe = 0, length_mfe = 0;
		for (i=0; 9*i < tokens.size(); i++) {  
		    length_hishape = tokens[9*i+2].length();
		    length_mfe = tokens[9*i+4].length();
		    if (length_longest_hishape < length_hishape) {
			length_longest_hishape = length_hishape;
		    }
		    if (length_longest_mfe < length_mfe) {
			length_longest_mfe = length_mfe;
		    }
		}
		
		
		std::map<std::string, double> hishape_pf;
		std::string token_9_2 = "";
		// prepare format_pattern and output the results
		std::string format_pattern = ( boost::format("%%%d.2f %%%ds %%15.10f") % (length_longest_mfe+5) % (length_longest_hishape+4) ).str();
		i = 0;
		for (i=0; 9*i < tokens.size(); i++) {   
                    //std::string hishape = tokens[9*i+2];
                    double pf = 0.0;
		    if (hishape_pf.find(tokens[9*i+2]) != hishape_pf.end()) {
			pf = hishape_pf[tokens[9*i+2]];  
		    } else {
                        //match_str = tokens[9*i+2]
                        std::string match_str_for_function = "";    
			std::vector<std::string> helixIndices;
			tokenize(tokens[9*i+2], helixIndices);
			std::vector<std::string>::iterator it;
			for ( it=helixIndices.begin() ; it != helixIndices.end(); it++ )
			{
			    std::string currentHi = *it;
			    if ( (currentHi[currentHi.size()-1] != 'b') && (currentHi[currentHi.size()-1] != 'i') && \
				  (currentHi[currentHi.size()-1] != 'm') && (currentHi != "(") && (currentHi != ")") )  // in case h or _
			    {
				match_str_for_function = match_str_for_function + "," + currentHi;
			    }
			}

		        calculateHishapeProbabilities(inputs,1000000000,hishape_type, match_str_for_function, true, &hishape_pf);
			pf = hishape_pf[tokens[9*i+2]]; 
			//std::cout << "calculateHishapeProbabilities("<<","<<1000000000<<","<<hishape_type<<","<<match_str_for_function<<")" << std::endl;
		    }
		    token_9_2 = tokens[9*i+2];
		    if (token_9_2 != "_")
		        token_9_2.erase (token_9_2.end()-1); 
		    std::cout << tokens[9*i+7] << (boost::format(format_pattern) % ((boost::lexical_cast<int>(tokens[9*i+4]))/100.0f) % ("["+token_9_2+"]") % pf).str() << std::endl;
		}
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


      /* delete all allocated memory */
      for (inputs_t::iterator i = inputs.begin(); i != inputs.end(); ++i)
	delete[] (*i).first;
  }
  
  //if (input != NULL)
  //  delete[] input;


  return 0;
}
