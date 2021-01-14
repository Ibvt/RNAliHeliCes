#include "hited_util.hh"


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
CHelixIndexTreeNode::CTreeNodePointer HishapeToTree(std::string &strHishape)
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
std::string sSToHishape( const std::string &strStructure1, int hishape_type)
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

// END FUNCTIONS
//******************************************************************************