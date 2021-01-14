#ifndef _HITED_UTIL_HH_
#define _HITED_UTIL_HH_


/*
 * hited.cc
 *
 *  Created on: June 1, 2012
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
#include <cassert>

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

/*
  extent the main.cc with a new option 'P' for probability calculation
*/
// hishape abstraction type
const unsigned int HISHAPE_B                = 1;
const unsigned int HISHAPE_M                = 2;
const unsigned int HISHAPE_H_PLUS           = 3;
const unsigned int HISHAPE_H                = 4;

/* functions used in common */
void conflicting_options(const po::variables_map& vm,
                         const char* opt1, const char* opt2);

/* Function used to check that of 'for_what' is specified, then
   'required_option' is specified too. */
void option_dependency(const po::variables_map& vm,
                        const char* for_what, const char* required_option);

void tokenize(const std::string& str,
                      std::vector<std::string>& tokens,
                      const std::string& delimiters);

template <class T>
bool from_string(T& t, 
                 const std::string& s, 
                 std::ios_base& (*f)(std::ios_base&));

typedef val::biocpp::tools::CTreeNode<hishape::comp::CHelixIndex> CHelixIndexTreeNode;
//! A tree node that contains a pair of nucleotides.
CHelixIndexTreeNode::CTreeNodePointer HishapeToTree(std::string &strHishape);

std::string PrintStemLoop(const val::biocpp::CNASecondaryStructure::SSubstructureData &tRetData, const val::biocpp::CNASecondaryStructure &tSecStructure, unsigned int hishape_type);


std::string PrintStem(const val::biocpp::CNASecondaryStructure::SSubstructureData &tRetData, const val::biocpp::CNASecondaryStructure &tSecStructure, unsigned int hishape_type);

typedef val::biocpp::tools::CTreeNode<hishape::comp::CHelixIndex> CHelixIndexTreeNode;
std::string PrintRecursive( const CHelixIndexTreeNode::CTreeNodePointer &treeNode, const val::biocpp::CNASecondaryStructure &tSecStructure, const std::multimap<int,int> &mmFirstPosJ, unsigned int hishape_type);

// transform secondary structure to hishape
std::string sSToHishape( const std::string &strStructure1, int hishape_type);

int pairAlign(std::string strInput1, std::string strInput2, int hishape_type, float ratio);

#endif

