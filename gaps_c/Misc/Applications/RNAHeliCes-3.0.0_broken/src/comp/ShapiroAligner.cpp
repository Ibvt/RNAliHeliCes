//*****************************************************************************
#include "ShapiroAligner.hpp"
#include "DebugTools.hpp"
//#include "HishapeToTreeConverter.hpp"
#include <stack>
#include <functional>
#include <vector>
#include <map>
#include <sstream>
//#include <stdlib.h>
#include <math.h>

#include <boost/algorithm/string/replace.hpp>

namespace hishape
{
namespace comp
{
namespace algorithms
{
/**
Standard constructor
********************
Creates a new uninitialized instance of the algorithm. Once created, the user
should set the 2 structures to compare.
*/
CShapiroAligner::CShapiroAligner(void)
{
	m_bInited = false;
}
CShapiroAligner::CShapiroAligner(const CHelixIndexTreeNode::CTreeNodePointer &ptTree1, const CHelixIndexTreeNode::CTreeNodePointer &ptTree2, float ratio/*, int hishape_type*/)
{
  /*
        m_tFirst.m_iLength = hishape1.m_iLength;
	m_tFirst.m_tNucleotidesPairs = hishape1.m_tNucleotidesPairs;
	m_tFirst.m_helixIndices = hishape1.m_helixIndices;	
	m_tFirst.m_tNucleotidesConfig = hishape1.m_tNucleotidesConfig;
	std::cout << hishape1.GetPairedBaseIndex(0) << std::endl;
	std::cout << m_tFirst.GetPairedBaseIndex(0) << std::endl;
	
	m_tSecond.m_iLength = hishape2.m_iLength;
	m_tSecond.m_tNucleotidesPairs = hishape2.m_tNucleotidesPairs;
	m_tSecond.m_helixIndices = hishape2.m_helixIndices;	
	m_tSecond.m_tNucleotidesConfig = hishape2.m_tNucleotidesConfig;
	std::cout << m_tSecond.GetPairedBaseIndex(0) << std::endl;
	*/
  
  
  
	m_tFirst = ptTree1;
	m_tSecond = ptTree2;
	m_ratio = ratio;
	/*m_hishapeType = hishape_type;*/
	m_bInited = false;
}


/**
Destructor
**********
A non-virtual destructor (this class shouldn't be derivated).
*/
CShapiroAligner::~CShapiroAligner(void)
{
}


/**
@fn const float& GetMin(const float& fVal1, const float& fVal2)
***************************************************************
Returns the minimum value between the 2 given floats.

@param const float& fFirstValue:
  the first float value.

@param const float& fSecondValue:
  the second float value.

@return const float&:
  fFirstValue if fFirstValue is less than fSecondValue, fSecondValue otherwise.

*/
const float& CShapiroAligner::GetMin(const float& fFirstValue, const float& fSecondValue)
{
	return (fFirstValue < fSecondValue)? fFirstValue:fSecondValue;
}


/**
@fn const float& GetMax(const float& fVal1, const float& fVal2)
***************************************************************
Returns the maximum value between the 2 given floats.

@param const float& fFirstValue:
  the first float value.

@param const float& fSecondValue:
  the second float value.

@return const float&:
  fFirstValue if fFirstValue is greater than fSecondValue, fSecondValue otherwise.

*/
const float& CShapiroAligner::GetMax(const float& fFirstValue, const float& fSecondValue)
{
	return (fFirstValue > fSecondValue)? fFirstValue:fSecondValue;
} 
    

/*
@fn CShapiroAligner::CHelixIndexTreeNode::CTreeNodePointerList CShapiroAligner::GetLRKeyRoots(CHelixIndexTreeNode::CTreeNodePointer ptNode) const
*********************************************************************************************************************************************************************
Returns the list of keyroots as they are described in Zhang-Shasha algorithm.

@param CHelixIndexTreeNode& tNode:\n
 tree node to start from.

@return CHelixIndexTreeNode::CTreeNodePointerList:\n
 the list of key roots (node pointers in a list).
*/
// though GetLRKeyRoots is private, it doesn't need to set private keyword in cpp file
CShapiroAligner::CHelixIndexTreeNode::CTreeNodePointerList CShapiroAligner::GetLRKeyRoots(CHelixIndexTreeNode::CTreeNodePointer ptNode) const
{
	CHelixIndexTreeNode::CTreeNodePointerList tRootNodesList;
	CHelixIndexTreeNode::CTreeNodePointerList tSuffixList = ptNode->GetSuffixList();
	//ptNode->PrintSuffixList();
	
	for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter = tSuffixList.begin(); tIter != tSuffixList.end(); tIter++)
	{
		if (true == (*tIter)->HasParent())
		{
			// if the node is not the left child of its parent, add it
			if (false == (*tIter)->IsDirectLeftChildOf( (*tIter)->GetParentPointer() ))
			{
				tRootNodesList.push_back(*tIter);
			}
		}
		else
		{
			// tree root node
			tRootNodesList.push_back(*tIter);
		}
	}
	return tRootNodesList;
}



/**
@fn int CShapiroAligner::Align(bool bStoreDistMatrix, bool bComputeBacktracking)
**********************************************************************************************
Returns the alignment score between two secondary structures.

@param bool bStoreDistMatrix:
 enable distance matrix storing if true, otherwise save memory.

@param bool bComputeBacktracking:
 enable backtracking computation if true, otherwise save memory and speed.

@return int:
 the alignement score.
*/
float CShapiroAligner::Align(bool bStoreDistMatrix, bool bComputeBacktracking)
{
	if (true == bComputeBacktracking)
	{
		return AlignWithBacktracking(bStoreDistMatrix);
	}
	else
	{
		return AlignNoBacktracking(bStoreDistMatrix);
	}
}


/**
@fn int CShapiroAligner::AlignNoBacktracking(bool bStoreDistMatrix)
******************************************************************
Returns the alignment score between two secondary structures.

@param bool bStoreDistMatrix:
 enable distance matrix storing if true, otherwise save memory.

@return int:
 the alignement score.
*/
float CShapiroAligner::AlignNoBacktracking(bool bStoreDistMatrix)
{
	m_iCyclesCount = 0;
	std::ostringstream tOutputStream;

	float fDistance = -1.0f;
	
	
	// ## Hishapes ==> Trees ##
	CHelixIndexTreeNode::CTreeNodePointer ptTree1 = m_tFirst;
	CHelixIndexTreeNode::CTreeNodePointer ptTree2= m_tSecond;
	
	
	// ## Trees ==> GetLRKeyRoots ##
	/*
	    CHelixIndexTreeNode::CTreeNodePointerList tGetLRKeyRoots1 = GetLRKeyRoots(ptTree1);
	    CHelixIndexTreeNode::CTreeNodePointerList tGetLRKeyRoots2 = GetLRKeyRoots(ptTree2);
	    CHelixIndexTreeNode::CTreeNodePointerList tSuffixNodesList1 = ptTree1->GetSuffixList();
	    CHelixIndexTreeNode::CTreeNodePointerList tSuffixNodesList2 = ptTree2->GetSuffixList();
	*/
//	if (m_bInited == false)
//	{
		// init stuff
		// compute GetLRKeyRoots for each tree
		CHelixIndexTreeNode::CTreeNodePointerList tGetLRKeyRoots1 = GetLRKeyRoots(ptTree1);
//#ifdef _DEBUG
/*
		TRACED("GetLRKeyRoots for tree 1 returned " << (signed) tGetLRKeyRoots1.size() << " nodes:");
		std::cout << "GetLRKeyRoots for tree 1 returned " << (signed) tGetLRKeyRoots1.size() << " nodes:" << std::endl;
		std::string strKeyRoots1("[");
		for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tRootIter1 = tGetLRKeyRoots1.begin(); tRootIter1 != tGetLRKeyRoots1.end(); tRootIter1++)
		{
			//strKeyRoots1 += (*tRootIter1)->GetObject().GetPairBaseCodes() + ", ";
			strKeyRoots1 += (*tRootIter1)->GetObject().ToString() + ", ";
		}
		strKeyRoots1.resize((signed)strKeyRoots1.size() - 2);
		strKeyRoots1 += "]";
		TRACED(strKeyRoots1);
		std::cout << strKeyRoots1 << std::endl;
		*/
//#else
//#endif
		CHelixIndexTreeNode::CTreeNodePointerList tGetLRKeyRoots2 = GetLRKeyRoots(ptTree2);
//#ifdef _DEBUG
/*
		TRACED("GetLRKeyRoots for tree 2 returned " << (signed) tGetLRKeyRoots2.size() << " nodes:");
		std::cout << "GetLRKeyRoots for tree 2 returned " << (signed) tGetLRKeyRoots2.size() << " nodes:" << std::endl;
		std::string strKeyRoots2("[");
		for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tRootIter2 = tGetLRKeyRoots2.begin(); tRootIter2 != tGetLRKeyRoots2.end(); tRootIter2++)
		{
			//strKeyRoots2 += (*tRootIter2)->GetObject().GetPairBaseCodes() + ", ";
			strKeyRoots2 += (*tRootIter2)->GetObject().ToString() + ", ";
		}
		strKeyRoots2.resize((signed)strKeyRoots2.size() - 2);
		strKeyRoots2 += "]";
		TRACED(strKeyRoots2);
		std::cout << strKeyRoots2 << std::endl;
		*/
//#else
//#endif
		CHelixIndexTreeNode::CTreeNodePointerList tSuffixNodesList1 = ptTree1->GetSuffixList();
//#ifdef _DEBUG
/*
		TRACED("Suffix order for tree 1:");
		std::cout << "Suffix order for tree 1:" << std::endl;
		std::string strSuffix1("[");
		for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tSuffixIter1 = tSuffixNodesList1.begin(); tSuffixIter1 != tSuffixNodesList1.end(); tSuffixIter1++)
		{
			//strSuffix1 += (*tSuffixIter1)->GetObject().GetPairBaseCodes() + ", ";
			strSuffix1 += (*tSuffixIter1)->GetObject().ToString() + ", ";
		}
		strSuffix1.resize((signed)strSuffix1.size() - 2);
		strSuffix1 += "]";
		TRACED(strSuffix1);
		std::cout << strSuffix1 << std::endl;
		*/
//#else
//#endif
		CHelixIndexTreeNode::CTreeNodePointerList tSuffixNodesList2 = ptTree2->GetSuffixList();
//#ifdef _DEBUG
/*
		TRACED("Suffix order for tree 2:");
		std::cout << "Suffix order for tree 2:" << std::endl;
		std::string strSuffix2("[");
		for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tSuffixIter2 = tSuffixNodesList2.begin(); tSuffixIter2 != tSuffixNodesList2.end(); tSuffixIter2++)
		{
			//strSuffix2 += (*tSuffixIter2)->GetObject().GetPairBaseCodes() + ", ";
			strSuffix2 += (*tSuffixIter2)->GetObject().ToString() + ", ";
		}
		strSuffix2.resize((signed)strSuffix2.size() - 2);
		strSuffix2 += "]";
		TRACED(strSuffix2);
		std::cout << strSuffix2 << std::endl;
		*/
//#else
//#endif
//		m_bInited = true;
//	}


	// ## comparing Trees ##
	
	// Compute lookup tables to assign a node its index in the initial tree
	std::map<CHelixIndexTreeNode::CTreeNodePointer, int, SNodePtrCompare> tNodeIndex1;
	std::map<CHelixIndexTreeNode::CTreeNodePointer, int, SNodePtrCompare> tNodeIndex2;
	// initialize map 1 (lookup table for tree 1)
	int iNodeIndex = 0;
	for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter1 = tSuffixNodesList1.begin(); tIter1 != tSuffixNodesList1.end(); tIter1++)
	{
		tNodeIndex1[(*tIter1)] = iNodeIndex++;
		//std::cout << "ShapiroAligner:" << (*tIter1)->GetObject().ToString() << "," << iNodeIndex << std::endl;
	}
	// initialize map 2 (lookup table for tree 2)
	iNodeIndex = 0;
	for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter2 = tSuffixNodesList2.begin(); tIter2 != tSuffixNodesList2.end(); tIter2++)
	{
		tNodeIndex2[(*tIter2)] = iNodeIndex++;
	}
	
	
	// local tree index node lookup tables
	// used while computing forestdist
	// note: these lookup tables do not contain the real index of a node in a
	//       local forest but its real index minus one for optimisation in the
	//       way the lookup table is used.
	std::map<CHelixIndexTreeNode::CTreeNodePointer, int, SNodePtrCompare> tLocalTreeNodeIndex1;
	std::map<CHelixIndexTreeNode::CTreeNodePointer, int, SNodePtrCompare> tLocalTreeNodeIndex2;

	// treedist table (permanent array)
	
	int iTree1NodeCount = ptTree1->GetChildrenCount() + 1; // "+1" for root node
	int iTree2NodeCount = ptTree2->GetChildrenCount() + 1; // "+1" for root node
	// allocate memory for treedist table
	// declares that aafTreeDist[i] is a pointer to an array of iTree1NodeCount floats.
	// aafTreeDist is the set of such pointers
	// treeDist needs only
	float **aafTreeDist = NULL;
	aafTreeDist = new float*[iTree1NodeCount];
	for (int i=0; i < iTree1NodeCount; i++)
	{
		aafTreeDist[i] = new float[iTree2NodeCount];
	}
	// allocate memory for forestdist tables (max space)
	float **aafForestDist = new float*[iTree1NodeCount+1];
	for (int i=0; i<=iTree1NodeCount; i++)
	{
		aafForestDist[i] = new float[iTree2NodeCount+1];
	}


	// first line and column initialisation iterator
	CHelixIndexTreeNode::CTreeNodePointerList::iterator tInitIter;
	// loops on keyroots
	// GetLRKeyRoots for tree 1 returned 3 nodes: [[4,c], [4,e], [0,N]]
        // GetLRKeyRoots for tree 2 returned 3 nodes: [[4,b], [4,e], [0,N]]
	for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter1 = tGetLRKeyRoots1.begin(); tIter1 != tGetLRKeyRoots1.end(); tIter1++)
	{
	        // GetLRKeyRoots for tree 2 returned 1 nodes: [[0,N]]
		for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter2 = tGetLRKeyRoots2.begin(); tIter2 != tGetLRKeyRoots2.end(); tIter2++)
		{

			tSuffixNodesList1 = (*tIter1)->GetSuffixList();

			//std::cout << "ShapiroAligner:tSuffixNodesList1 of " << (*tIter1)->GetObject().ToString() << " is " << tSuffixNodesList1.ToString() << "." << std::endl;
			// gets a [l(j) to j] list
			tSuffixNodesList2 = (*tIter2)->GetSuffixList();

			// gets table dimensions
			int iTableHeight = (*tIter1)->GetChildrenCount() + 2; // "+1" for local root node "+1" for del/insert
			int iTableWidth = (*tIter2)->GetChildrenCount() + 2; // "+1" for local root node "+1" for del/insert

			float fDeletion = 0;
			float fInsertion = 0;
			float fOtherNodeDistance = 0;
			int iForestDistIndex1, iForestDistIndex2;
			
			// ###########################################################################
			// ######## here starts treedist(i,j) computation for following pairs ########
			// ###########################################################################
			/*
			[[4,b], [4,c]]                                 vs [[4,b]]
			[[4,b], [4,c]]                                 vs [[4,e]]
			[[4,b], [4,c]]                                 vs [[4,a], [4,b], [4,m], [4,c], [4,e], [0,N]]
			[[4,e]]                                        vs [[4,b]]
			[[4,e]]                                        vs [[4,e]]
			[[4,e]]                                        vs [[4,a], [4,b], [4,m], [4,c], [4,e], [0,N]]
			[[4,a], [4,b], [4,c], [4,m], [4,e], [0,N]]     vs [[4,b]]
			[[4,a], [4,b], [4,c], [4,m], [4,e], [0,N]]     vs [[4,e]]
			[[4,a], [4,b], [4,c], [4,m], [4,e], [0,N]]     vs [[4,a], [4,b], [4,m], [4,c], [4,e], [0,N]]
			*/
			// SCORE CALCULATION 0
			// forestdist(0,0)=0;
			aafForestDist[0][0] = 0;
			tInitIter = tSuffixNodesList1.begin();
			// for i1:=l(i) to i
			for (int i=1; i < iTableHeight; i++)
			{
			        // SCORE CALCULATION 1
				// forestdist(T1[l(i)..i1],0)=forestdist(T1[l(i)..i1-1],0)+@(T1[i1]->^)    DELETE
				aafForestDist[i][0] = aafForestDist[i-1][0] + m_ptScoreScheme->GetScore((*tInitIter)->GetObject().GetIntCode(), (*tInitIter)->GetObject().GetGapCode());
				tInitIter++;
			}
			tInitIter = tSuffixNodesList2.begin();
			// for j1:=l(j) to j
			for (int j=1; j < iTableWidth; j++)
			{
			        // SCORE CALCULATION 2
				// forestdist(0,T2[l(j)..j1])=forestdist(0,T2[l(j)..j1-1],0)+@(^->T2[j1])
				aafForestDist[0][j] = aafForestDist[0][j-1]  + m_ptScoreScheme->GetScore((*tInitIter)->GetObject().GetGapCode(), (*tInitIter)->GetObject().GetIntCode());
				tInitIter++;
			}
			// prepares local tree nodes index lookup table
			tLocalTreeNodeIndex1.clear();
			// for i1:=l(i) to i
			int i1=0;
			for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter3 = tSuffixNodesList1.begin(); tIter3 != tSuffixNodesList1.end(); tIter3++)
			{
				tLocalTreeNodeIndex1[(*tIter3)] = i1++; // note: post-increment for optimization (no further need to do a "tLocalTreeNodeIndex1[theNode] - 1")
				// prepares local tree nodes index lookup table
				tLocalTreeNodeIndex2.clear();
				int j1=0;
				// for j1:=l(j) to j
				for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter4 = tSuffixNodesList2.begin(); tIter4 != tSuffixNodesList2.end(); tIter4++)
				{
					tLocalTreeNodeIndex2[(*tIter4)] = j1++; // note: post-increment for optimization (no further need to do a "tLocalTreeNodeIndex1[theNode] - 1")
					iForestDistIndex1 = i1-1;
					// SCORE CALCULATION 3
					// delete on seconde sequence
					fDeletion = aafForestDist[iForestDistIndex1][j1] + m_ptScoreScheme->GetScore((*tIter3)->GetObject().GetIntCode(), (*tIter3)->GetObject().GetGapCode());
					iForestDistIndex2 = j1-1;
					// SCORE CALCULATION 4
					fInsertion = aafForestDist[i1][iForestDistIndex2] + m_ptScoreScheme->GetScore((*tIter4)->GetObject().GetGapCode(), (*tIter4)->GetObject().GetIntCode());
					m_iCyclesCount++;
					// if l(i1)=l(i) and l(j1)=l(j) then
					if (((*tIter3)->GetLeftLeafPointer() == (*tIter1)->GetLeftLeafPointer())
						&& ((*tIter4)->GetLeftLeafPointer() == (*tIter2)->GetLeftLeafPointer()))
					{
					        // SCORE CALCULATION 5
						// forestdist(T1[l(i)..i1],T2[l(j)..j1]) = min {
						//    forestdist(T1[l(i)..i1-1],T2[l(j)..j1])+@(T1[i1]->^),
						//    forestdist(T1[l(i)..i1],T2[l(j)..j1-1])+@(^->T2[j1]),
						//    forestdist(T1[l(i)..i1-1],T2[l(j)..j1-1])+@(T1[i1]->T2[j1])}
						if ((*tIter3)->GetObject().GetHiNumber()!=-1.0 && (*tIter4)->GetObject().GetHiNumber()!=-1.0)
						{
						    fOtherNodeDistance = aafForestDist[iForestDistIndex1][iForestDistIndex2] + \
									m_ptScoreScheme->GetScore((*tIter3)->GetObject().GetIntCode(), (*tIter4)->GetObject().GetIntCode()) + \
									fabs((*tIter3)->GetObject().GetHiNumber() - (*tIter4)->GetObject().GetHiNumber()) * m_ratio;
						}
						else
						{
						    fOtherNodeDistance = aafForestDist[iForestDistIndex1][iForestDistIndex2] + \
						    m_ptScoreScheme->GetScore((*tIter3)->GetObject().GetIntCode(), (*tIter4)->GetObject().GetIntCode());  
						}
						//## select the smallest value of three
						aafForestDist[i1][j1] = GetMin(fDeletion, GetMin(fInsertion, fOtherNodeDistance));
						// treedist(i1,j1)=forestdist(T1[l(i1)..i1],T2[l(j)..j1]) // put in permanent array
						aafTreeDist[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]] = aafForestDist[i1][j1];
					}
					else
					{
						// forestdist(T1[l(i)..i1],T2[l(j)..j1]) = min {
						//    forestdist(T1[l(i)..i1-1],T2[l(j)..j1])+@(T1[i1]->^),
						//    forestdist(T1[l(i)..i1],T2[l(j)..j1-1])+@(^->T2[j1]),
						//    forestdist(T1[l(i)..l(i1)-1],T2[l(j)..l(j1)-1]) + treedist(i1,j1)}
						if (tNodeIndex1[(*tIter3)->GetLeftLeafPointer()]-1 < tNodeIndex1[(*tIter1)->GetLeftLeafPointer()])
						{
							iForestDistIndex1 = 0;
						}
						else
						{
							iForestDistIndex1 = tLocalTreeNodeIndex1[(*tIter3)->GetLeftLeafPointer()];
						}
						if (tNodeIndex2[(*tIter4)->GetLeftLeafPointer()]-1 < tNodeIndex2[(*tIter2)->GetLeftLeafPointer()])
						{
							iForestDistIndex2 = 0;
						}
						else
						{
							iForestDistIndex2 = tLocalTreeNodeIndex2[(*tIter4)->GetLeftLeafPointer()];
						}
						fOtherNodeDistance = aafForestDist[iForestDistIndex1][iForestDistIndex2] + aafTreeDist[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]];
						aafForestDist[i1][j1] = GetMin(fDeletion, GetMin(fInsertion, fOtherNodeDistance));
					}
				}
			}
			// #### here ends treedist(i,j) computation ####
			// dumps current forestdist matrix
			
			
			//std::cout << "forest_dist (" << (*tIter1)->GetObject().ToString()/*GetPairBaseCodes()*/ << "," << (*tIter2)->GetObject().ToString()/*GetPairBaseCodes()*/ << ")" << std::endl;
			/*for (int i=0; i<iTableHeight; i++)
			{
				for (int j=0; j<iTableWidth; j++)
				{
					std::cout << aafForestDist[i][j] << " ";
				}
				std::cout << std::endl;
			}
			//tOutputStream << std::endl;
			std::cout << std::endl;
			*/
		}
	}

	fDistance = aafTreeDist[iTree1NodeCount-1][iTree2NodeCount-1];
	// dumps current treedist matrix
	/*
	std::cout << "tree_dist" << std::endl;
	for (int i=0; i<iTree1NodeCount; i++)
	{
		for (int j=0; j<iTree2NodeCount; j++)
		{
			std::cout << aafTreeDist[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	*/
	
	// due to the trees were constructed in own data structure, they have to be freed after using. 
	// in compare with this, every data contained in c++ container, for example vector. they don't need to be freed.
	
	// frees forestdist table memory
	for (int i=0; i<=iTree1NodeCount; i++)
	{
		delete [] aafForestDist[i];
	}
	delete [] aafForestDist;

	// frees treedist table memory
	for (int i=0; i<iTree1NodeCount; i++)
	{
		delete [] aafTreeDist[i];
	}
	delete [] aafTreeDist;

	// free trees
	delete ptTree1;
	delete ptTree2;

	m_strOutput = tOutputStream.str();
	m_fLastScore = fDistance;

	
	return fDistance;
}


/**
@fn int CShapiroAligner::AlignWithBacktracking(bool bStoreDistMatrix)
******************************************************************
Returns the alignment score between two secondary structures.

@param bool bStoreDistMatrix:
 enable distance matrix storing if true, otherwise save memory.

@return int:
 the alignement score.
*/
float CShapiroAligner::AlignWithBacktracking(bool bStoreDistMatrix)
{
	m_iCyclesCount = 0;
	std::ostringstream tOutputStream;

	float fDistance = -1;
	
	
	// ## Hishapes ==> Trees ##
	CHelixIndexTreeNode::CTreeNodePointer ptTree1 = m_tFirst;
	CHelixIndexTreeNode::CTreeNodePointer ptTree2= m_tSecond;
	
	
	
	// ## Trees ==> GetLRKeyRoots ##
	/*
	    CHelixIndexTreeNode::CTreeNodePointerList tGetLRKeyRoots1 = GetLRKeyRoots(ptTree1);
	    CHelixIndexTreeNode::CTreeNodePointerList tGetLRKeyRoots2 = GetLRKeyRoots(ptTree2);
	    CHelixIndexTreeNode::CTreeNodePointerList tSuffixNodesList1 = ptTree1->GetSuffixList();
	    CHelixIndexTreeNode::CTreeNodePointerList tSuffixNodesList2 = ptTree2->GetSuffixList();
	*/
//	if (m_bInited == false)
//	{
		// init stuff
		// compute GetLRKeyRoots for each tree
		CHelixIndexTreeNode::CTreeNodePointerList tGetLRKeyRoots1 = GetLRKeyRoots(ptTree1);
		CHelixIndexTreeNode::CTreeNodePointerList tGetLRKeyRoots2 = GetLRKeyRoots(ptTree2);

		CHelixIndexTreeNode::CTreeNodePointerList tSuffixNodesList1 = ptTree1->GetSuffixList();
		
		TRACED("Suffix order for tree 1:");
		std::cout << "Suffix order for tree 1:" << std::endl;
		std::string strSuffix1("[");
		for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tSuffixIter1 = tSuffixNodesList1.begin(); tSuffixIter1 != tSuffixNodesList1.end(); tSuffixIter1++)
		{
			//strSuffix1 += (*tSuffixIter1)->GetObject().GetPairBaseCodes() + ", ";
			strSuffix1 += (*tSuffixIter1)->GetObject().ToString() + ", ";
		}
		strSuffix1.resize((signed)strSuffix1.size() - 2);
		strSuffix1 += "]";
		TRACED(strSuffix1);
		std::cout << strSuffix1 << std::endl;
		
		CHelixIndexTreeNode::CTreeNodePointerList tSuffixNodesList2 = ptTree2->GetSuffixList();
		
		TRACED("Suffix order for tree 2:");
		std::cout << "Suffix order for tree 2:" << std::endl;
		std::string strSuffix2("[");
		for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tSuffixIter2 = tSuffixNodesList2.begin(); tSuffixIter2 != tSuffixNodesList2.end(); tSuffixIter2++)
		{
			//strSuffix2 += (*tSuffixIter2)->GetObject().GetPairBaseCodes() + ", ";
			strSuffix2 += (*tSuffixIter2)->GetObject().ToString() + ", ";
		}
		strSuffix2.resize((signed)strSuffix2.size() - 2);
		strSuffix2 += "]";
		TRACED(strSuffix2);
		std::cout << strSuffix2 << std::endl;
		

//		m_bInited = true;
//	}


	// ## comparing Trees ##
	//###################################################
	//########## prepare lookup maps ####################
	//###################################################
	// Compute lookup tables to assign a node its index in the initial tree
	std::map<CHelixIndexTreeNode::CTreeNodePointer, int, SNodePtrCompare> tNodeIndex1;
	std::map<int, CHelixIndexTreeNode::CTreeNodePointer, IntegerCompare> tIndexNode1;
	std::map<CHelixIndexTreeNode::CTreeNodePointer, int, SNodePtrCompare> tNodeIndex2;
	std::map<int, CHelixIndexTreeNode::CTreeNodePointer, IntegerCompare> tIndexNode2;
	// initialize map 1 (lookup table for tree 1)
	int iNodeIndex = 0;
	for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter1 = tSuffixNodesList1.begin(); tIter1 != tSuffixNodesList1.end(); tIter1++)
	{       
		tNodeIndex1[(*tIter1)] = iNodeIndex;
		tIndexNode1[iNodeIndex] = (*tIter1);
		iNodeIndex++;
		//std::cout << "ShapiroAligner:" << (*tIter1)->GetObject().ToString() << "," << iNodeIndex << std::endl;
	}
	// initialize map 2 (lookup table for tree 2)
	iNodeIndex = 0;
	for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter2 = tSuffixNodesList2.begin(); tIter2 != tSuffixNodesList2.end(); tIter2++)
	{
		tNodeIndex2[(*tIter2)] = iNodeIndex;
		tIndexNode2[iNodeIndex] = (*tIter2);
                iNodeIndex++;
	}
	
	
	// local tree index node lookup tables
	// used while computing forestdist
	// note: these lookup tables do not contain the real index of a node in a
	//       local forest but its real index minus one for optimisation in the
	//       way the lookup table is used.
	std::map<CHelixIndexTreeNode::CTreeNodePointer, int, SNodePtrCompare> tLocalTreeNodeIndex1;
	std::map<CHelixIndexTreeNode::CTreeNodePointer, int, SNodePtrCompare> tLocalTreeNodeIndex2;


	
	
	//###################################################
	//########## global containers ####################
	//###################################################
	int iTree1NodeCount = ptTree1->GetChildrenCount() + 1; // "+1" for root node
	int iTree2NodeCount = ptTree2->GetChildrenCount() + 1; // "+1" for root node
	// allocate memory for treedist table
	// declares that aafTreeDist[i] is a pointer to an array of iTree1NodeCount floats.
	// aafTreeDist is the set of such pointers
	// treeDist needs only
	/*float **aafTreeDist = NULL;
	aafTreeDist = new float*[iTree1NodeCount];
	for (int i=0; i < iTree1NodeCount; i++)
	{
		aafTreeDist[i] = new float[iTree2NodeCount];
	}*/
	// treedist table (permanent array)
	std::vector<std::vector<float> > aafTreeDist(iTree1NodeCount);
	for (int i=0; i < iTree1NodeCount; i++)
	{
		aafTreeDist[i].resize(iTree2NodeCount);
	}
	//@@ initialize aafTreeDistPath
	std::vector<std::vector<  std::list<  std::pair<std::pair<std::pair<CHelixIndexTreeNode::CTreeNodePointer, CHelixIndexTreeNode::CTreeNodePointer>, std::pair<int,int> >, std::pair<char, std::pair<int, int> > >  >  > > aafTreeDistPath(iTree1NodeCount);
	for (int i=0; i < iTree1NodeCount; i++)
	{
		aafTreeDistPath[i].resize(iTree2NodeCount);
	}

	// save all forestdist in a 4 D matrix
	std::vector<std::vector< std::vector<std::vector<float> >  > > aafForestDist(iTree1NodeCount);
	for (int i=0; i<iTree1NodeCount; i++)
	{
	        aafForestDist[i].resize(iTree2NodeCount);
		for (int j=0; j<iTree2NodeCount; j++)
		{
		        aafForestDist[i][j].resize(iTree1NodeCount+1);
		    	for (int k=0; k <= iTree1NodeCount; k++)
			{
				aafForestDist[i][j][k].resize(iTree2NodeCount+1);
			}
		}
	}
	// allocate memory for forestdist tables (max space)
	/*float **aafForestDist = new float*[iTree1NodeCount+1];
	for (int i=0; i<=iTree1NodeCount; i++)
	{
		aafForestDist[i] = new float[iTree2NodeCount+1];
	}*/
	
	
	//## data structure for backtracking
	std::vector<std::vector<std::vector<std::vector<std::list< std::pair<std::pair<std::pair<int, int>, char>, std::pair<int, int> >  >  > > > > tBacktrackingTables(iTree1NodeCount);
	for (int i=0; i<iTree1NodeCount; i++)
	{
	        tBacktrackingTables[i].resize(iTree2NodeCount);
		for (int j=0; j<iTree2NodeCount; j++)
		{
		        tBacktrackingTables[i][j].resize(iTree1NodeCount+1);
			for (int k=0; k<=iTree1NodeCount; k++)
			{
				tBacktrackingTables[i][j][k].resize(iTree2NodeCount+1);
			}
			tBacktrackingTables[i][j][0][0].push_back(make_pair(make_pair(std::pair<int, int>(0, 0), 'N'), std::pair<int,int>(0,0)));
		}
	}

	//@@ initialize aafForestDistPath
	/*
	std::vector<std::vector< std::vector<std::vector<  std::list< std::pair<std::pair<int, int>, char> >  > > > > aafForestDistPath(iTree1NodeCount);
	for (int i=0; i<iTree1NodeCount; i++)
	{
	        aafForestDistPath[i].resize(iTree2NodeCount);
		for (int j=0; j<iTree2NodeCount; j++)
		{
		        aafForestDistPath[i][j].resize(iTree1NodeCount+1);
			for (int k=0; k<=iTree1NodeCount; k++)
			{
				aafForestDistPath[i][j][k].resize(iTree2NodeCount+1);
			}
		}
	}*/
	/*
	std::vector<std::vector< std::vector<std::vector<  std::list< std::pair<std::pair<int, int>, char> >  > > > > aafForestDistPath(iTree1NodeCount);
	for (int i=0; i<iTree1NodeCount; i++)
	{
	        aafForestDistPath[i].resize(iTree2NodeCount);
		for (int j=0; j<iTree2NodeCount; j++)
		{
		        aafForestDistPath[i][j].resize(iTree1NodeCount+1);
			for (int k=0; k<=iTree1NodeCount; k++)
			{
				aafForestDistPath[i][j][k].resize(iTree2NodeCount+1);
			}
		}
	}*/                                                            
	                                                                                                                                       /* previous point */                                  /* replace type */  /* next point */
	std::vector<std::vector<std::vector<std::vector< std::list< std::pair<std::pair<std::pair<CHelixIndexTreeNode::CTreeNodePointer, CHelixIndexTreeNode::CTreeNodePointer>, std::pair<int,int> >, std::pair<char, std::pair<int, int> > >  >  > > > > aafForestDistPath(iTree1NodeCount);
	for (int i=0; i<iTree1NodeCount; i++)
	{
	        aafForestDistPath[i].resize(iTree2NodeCount);
		for (int j=0; j<iTree2NodeCount; j++)
		{
		        aafForestDistPath[i][j].resize(iTree1NodeCount+1);
			for (int k=0; k<=iTree1NodeCount; k++)
			{
				aafForestDistPath[i][j][k].resize(iTree2NodeCount+1);
			}
		}
	}


	//## global loopup table
	//std::vector<std::vector< std::pair<std::pair<CHelixIndexTreeNode::CTreeNodePointer, CHelixIndexTreeNode::CTreeNodePointer>, std::pair<int, int> >    > > tGlobalDistLookupTable(iTree1NodeCount);
	//BUG_TRACK  std::vector<std::vector< std::pair<std::pair<int, int>, std::pair<int, int> >    > > tGlobalDistLookupTable(iTree1NodeCount);  // first pair refers to tabel, second pair refers to the position in table
	//BUG_TRACK  for (int i=0; i<iTree1NodeCount; i++)
	//BUG_TRACK  {
	//BUG_TRACK          tGlobalDistLookupTable[i].resize(iTree2NodeCount);
	//BUG_TRACK  }
	
	
	
	
	//###################################################
	//######## two loops to fill global containers ######
	//###################################################
	// first line and column initialisation iterator
	CHelixIndexTreeNode::CTreeNodePointerList::iterator tInitIter;
	// loops on keyroots
	// GetLRKeyRoots for tree 1 returned 3 nodes: [[4,c], [4,e], [0,N]]
        // GetLRKeyRoots for tree 2 returned 3 nodes: [[4,b], [4,e], [0,N]]
	for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter1 = tGetLRKeyRoots1.begin(); tIter1 != tGetLRKeyRoots1.end(); tIter1++)
	{
	        // GetLRKeyRoots for tree 2 returned 1 nodes: [[0,N]]
		for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter2 = tGetLRKeyRoots2.begin(); tIter2 != tGetLRKeyRoots2.end(); tIter2++)
		{

			tSuffixNodesList1 = (*tIter1)->GetSuffixList();

			//std::cout << "ShapiroAligner:tSuffixNodesList1 of " << (*tIter1)->GetObject().ToString() << " is " << tSuffixNodesList1.ToString() << "." << std::endl;
			// gets a [l(j) to j] list
			tSuffixNodesList2 = (*tIter2)->GetSuffixList();

			// gets table dimensions
			int iTableHeight = (*tIter1)->GetChildrenCount() + 2; // "+1" for local root node "+1" for del/insert
			int iTableWidth = (*tIter2)->GetChildrenCount() + 2; // "+1" for local root node "+1" for del/insert

			float fDeletion = 0;
			float fInsertion = 0;
			float fOtherNodeDistance = 0;
			int iForestDistIndex1, iForestDistIndex2;
			
			// ###########################################################################
			// ######## here starts treedist(i,j) computation for following pairs ########
			// ###########################################################################
			/*
			[[4,b], [4,c]]                                 vs [[4,b]]
			[[4,b], [4,c]]                                 vs [[4,e]]
			[[4,b], [4,c]]                                 vs [[4,a], [4,b], [4,m], [4,c], [4,e], [0,N]]
			[[4,e]]                                        vs [[4,b]]
			[[4,e]]                                        vs [[4,e]]
			[[4,e]]                                        vs [[4,a], [4,b], [4,m], [4,c], [4,e], [0,N]]
			[[4,a], [4,b], [4,c], [4,m], [4,e], [0,N]]     vs [[4,b]]
			[[4,a], [4,b], [4,c], [4,m], [4,e], [0,N]]     vs [[4,e]]
			[[4,a], [4,b], [4,c], [4,m], [4,e], [0,N]]     vs [[4,a], [4,b], [4,m], [4,c], [4,e], [0,N]]
			*/
			// forestdist(0,0)=0;
			// SCORE CALCULATION 0
			aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][0][0] = 0;
			// 因为 tBacktrackingTables 初始化为 make_pair(make_pair(std::pair<int, int>(0, 0), 'N'), std::pair<int,int>(0,0)), 所以这里不需要赋值
			//@@ aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][0][0].push_back(make_pair(make_pair(std::pair<CHelixIndexTreeNode::CTreeNodePointer, CHelixIndexTreeNode::CTreeNodePointer>(*tIter1,*tIter2),std::pair<int, int>(0, 0)),make_pair('N',std::pair<int, int>(0, 0))));
			tInitIter = tSuffixNodesList1.begin();
			// for i1:=l(i) to i
			for (int i=1; i < iTableHeight; i++)
			{
				// forestdist(T1[l(i)..i1],0)=forestdist(T1[l(i)..i1-1],0)+@(T1[i1]->^)    DELETE
				// SCORE CALCULATION 1
				aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i][0] = aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i-1][0] + m_ptScoreScheme->GetScore((*tInitIter)->GetObject().GetIntCode(), (*tInitIter)->GetObject().GetGapCode());
				// ## clear the cell
				tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i][0].clear();
				//tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i][0].push_back(make_pair(std::pair<int, int>(i-1, 0), 'E'));
				tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i][0].push_back(make_pair(make_pair(std::pair<int, int>(i-1, 0), 'D'), std::pair<int,int>(0,0)));
				aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i][0] = aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i-1][0];
				aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i][0].push_back(make_pair(make_pair(std::pair<CHelixIndexTreeNode::CTreeNodePointer,CHelixIndexTreeNode::CTreeNodePointer>(*tIter1,*tIter2), std::pair<int,int>(i-1,0)), make_pair('D',std::pair<int,int>(i,0))));
				tInitIter++;
			}
			tInitIter = tSuffixNodesList2.begin();
			// for j1:=l(j) to j
			for (int j=1; j < iTableWidth; j++)
			{
				// forestdist(0,T2[l(j)..j1])=forestdist(0,T2[l(j)..j1-1],0)+@(^->T2[j1])
				// SCORE CALCULATION 2
				aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][0][j] = aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][0][j-1]  + m_ptScoreScheme->GetScore((*tInitIter)->GetObject().GetGapCode(), (*tInitIter)->GetObject().GetIntCode());
				// ## clear the cell
				tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][0][j].clear();
				//tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][0][j].push_back(make_pair(std::pair<int, int>(0, j-1), 'E'));
				tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][0][j].push_back(make_pair(make_pair(std::pair<int, int>(0, j-1), 'I'), std::pair<int,int>(0,0)));
				aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][0][j] = aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][0][j-1];
				aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][0][j].push_back(make_pair(make_pair(std::pair<CHelixIndexTreeNode::CTreeNodePointer,CHelixIndexTreeNode::CTreeNodePointer>(*tIter1,*tIter2), std::pair<int,int>(0,j-1)), make_pair('I',std::pair<int,int>(0,j))));
				tInitIter++;
			}
			// prepares local tree nodes index lookup table
			tLocalTreeNodeIndex1.clear();
			// for i1:=l(i) to i
			int i1=0;
			for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter3 = tSuffixNodesList1.begin(); tIter3 != tSuffixNodesList1.end(); tIter3++)
			{
				tLocalTreeNodeIndex1[(*tIter3)] = i1++; // note: post-increment for optimization (no further need to do a "tLocalTreeNodeIndex1[theNode] - 1")
				// prepares local tree nodes index lookup table
				tLocalTreeNodeIndex2.clear();
				int j1=0;
				// for j1:=l(j) to j
				for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter4 = tSuffixNodesList2.begin(); tIter4 != tSuffixNodesList2.end(); tIter4++)
				{
					tLocalTreeNodeIndex2[(*tIter4)] = j1++; // note: post-increment for optimization (no further need to do a "tLocalTreeNodeIndex1[theNode] - 1")
					iForestDistIndex1 = i1-1;
					// SCORE CALCULATION 3
					fDeletion = aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][iForestDistIndex1/*important*/][j1] + m_ptScoreScheme->GetScore((*tIter3)->GetObject().GetIntCode(), (*tIter3)->GetObject().GetGapCode());
					iForestDistIndex2 = j1-1;
					// SCORE CALCULATION 4
					fInsertion = aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][iForestDistIndex2/*important*/] + m_ptScoreScheme->GetScore((*tIter4)->GetObject().GetGapCode(), (*tIter4)->GetObject().GetIntCode());
					m_iCyclesCount++;
					
					
					// ## clear the cell
					tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].clear();
					// if l(i1)=l(i) and l(j1)=l(j) then:    in the same corresponding subtree
					if (((*tIter3)->GetLeftLeafPointer() == (*tIter1)->GetLeftLeafPointer())
						&& ((*tIter4)->GetLeftLeafPointer() == (*tIter2)->GetLeftLeafPointer()))
					{
						// forestdist(T1[l(i)..i1],T2[l(j)..j1]) = min {
						//    forestdist(T1[l(i)..i1-1],T2[l(j)..j1])+@(T1[i1]->^),
						//    forestdist(T1[l(i)..i1],T2[l(j)..j1-1])+@(^->T2[j1]),
						//    forestdist(T1[l(i)..i1-1],T2[l(j)..j1-1])+@(T1[i1]->T2[j1])}
						// SCORE CALCULATION 5
						if ((*tIter3)->GetObject().GetHiNumber()!=-1.0 && (*tIter4)->GetObject().GetHiNumber()!=-1.0)
						{
						    fOtherNodeDistance = aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][iForestDistIndex1/*important*/][iForestDistIndex2/*important*/] + \
									m_ptScoreScheme->GetScore((*tIter3)->GetObject().GetIntCode(), (*tIter4)->GetObject().GetIntCode()) + \
									fabs((*tIter3)->GetObject().GetHiNumber() - (*tIter4)->GetObject().GetHiNumber()) * m_ratio;
						}
						else
						{
						    fOtherNodeDistance = aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][iForestDistIndex1/*important*/][iForestDistIndex2/*important*/] + \
									m_ptScoreScheme->GetScore((*tIter3)->GetObject().GetIntCode(), (*tIter4)->GetObject().GetIntCode());
						}
						//## select the smallest value of three
						aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1] = GetMin(fDeletion, GetMin(fInsertion, fOtherNodeDistance));
						

						if (aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1] == fDeletion) {
						    //BUG_TRACK  std::cout <<"C1_D:["<< (*tIter1)->GetObject().ToString()<<","<<(*tIter2)->GetObject().ToString()<<","<<i1<<","<<j1 <<"]="<< (*tIter3)->GetObject().ToString() <<","<< (*tIter4)->GetObject().ToString() <<",("<< fDeletion <<","<< fInsertion <<","<< fOtherNodeDistance <<","<< aafTreeDist[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]]<<")" << std::endl;  // 6,8,0,3=110 + 
						    //tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].push_back(make_pair(std::pair<int, int>(i1-1, j1), 'd'));
						    tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].push_back(make_pair(make_pair(std::pair<int, int>(i1-1, j1), 'D'), std::pair<int,int>(0,0)));
						    aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1] = aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1-1][j1];
				                    aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].push_back(make_pair(make_pair(std::pair<CHelixIndexTreeNode::CTreeNodePointer,CHelixIndexTreeNode::CTreeNodePointer>(*tIter1,*tIter2), std::pair<int,int>(i1-1,j1)), make_pair('D',std::pair<int,int>(i1,j1))));
						    aafTreeDist[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]] = aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1];
						    aafTreeDistPath[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]] = aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1];
						}
						else if (aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1] == fInsertion) {
						    //BUG_TRACK  std::cout <<"C1_I:["<< (*tIter1)->GetObject().ToString()<<","<<(*tIter2)->GetObject().ToString()<<","<<i1<<","<<j1 <<"]="<< (*tIter3)->GetObject().ToString() <<","<< (*tIter4)->GetObject().ToString() <<",("<< fDeletion <<","<< fInsertion <<","<< fOtherNodeDistance <<","<< aafTreeDist[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]]<<")" << std::endl;  // 6,8,0,3=110 + 
						    //tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].push_back(make_pair(std::pair<int, int>(i1, j1-1), 'i'));
						    tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].push_back(make_pair(make_pair(std::pair<int, int>(i1, j1-1), 'I'), std::pair<int,int>(0,0)));
						    aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1] = aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1-1];
				                    aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].push_back(make_pair(make_pair(std::pair<CHelixIndexTreeNode::CTreeNodePointer,CHelixIndexTreeNode::CTreeNodePointer>(*tIter1,*tIter2), std::pair<int,int>(i1,j1-1)), make_pair('I',std::pair<int,int>(i1,j1))));
						    aafTreeDist[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]] = aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1];
						    aafTreeDistPath[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]] = aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1];
						}
						else if (aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1] == fOtherNodeDistance) {
						    //BUG_TRACK  std::cout <<"C1_R:["<< (*tIter1)->GetObject().ToString()<<","<<(*tIter2)->GetObject().ToString()<<","<<i1<<","<<j1 <<"]="<< (*tIter3)->GetObject().ToString() <<","<< (*tIter4)->GetObject().ToString() <<",("<< fDeletion <<","<< fInsertion <<","<< fOtherNodeDistance <<","<< aafTreeDist[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]]<<")" << std::endl;  // 6,8,0,3=110 + 
						    //tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].push_back(make_pair(std::pair<int, int>(i1-1, j1-1), 'r'));
						    //std::cout << fOtherNodeDistance <<" = (R) aafForestDist["<<tNodeIndex1[(*tIter1)]<<"]["<<tNodeIndex2[(*tIter2)]<<"]["<<iForestDistIndex1<<"]["<<iForestDistIndex2<<"]("<<aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][iForestDistIndex1][iForestDistIndex2] <<\
						    ") + "<< m_ptScoreScheme->GetScore((*tIter3)->GetObject().GetIntCode(), (*tIter4)->GetObject().GetIntCode()) <<"+"<< fabs((*tIter3)->GetObject().GetHiNumber() - (*tIter4)->GetObject().GetHiNumber())<< std::endl;
						    tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].push_back(make_pair(make_pair(std::pair<int, int>(i1-1, j1-1), 'R'), std::pair<int,int>(0,0)));
						    aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1] = aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1-1][j1-1];
				                    aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].push_back(make_pair(make_pair(std::pair<CHelixIndexTreeNode::CTreeNodePointer,CHelixIndexTreeNode::CTreeNodePointer>(*tIter1,*tIter2), std::pair<int,int>(i1-1,j1-1)), make_pair('R',std::pair<int,int>(i1,j1))));
						    aafTreeDist[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]] = aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1];
						    
						    
						    aafTreeDistPath[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]] = aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1];
						    // check intermediate results
						    //BUG_TRACK  std::list< std::pair<std::pair<std::pair<CHelixIndexTreeNode::CTreeNodePointer, CHelixIndexTreeNode::CTreeNodePointer>, std::pair<int,int> >, std::pair<char, std::pair<int, int> > >  >::iterator it;
						    //BUG_TRACK  for ( it=aafTreeDistPath[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]].begin() ; it != aafTreeDistPath[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]].end(); it++ )
						    //BUG_TRACK  {
						    //BUG_TRACK  std::cout << "check nach C1_R:[" << ((*it).first.first.first)->GetObject().ToString() << "," << ((*it).first.first.second)->GetObject().ToString() << "," << (*it).first.second.first << "," << (*it).first.second.second << "]" << "=" << (*it).second.first << "=>" << "[" << (*it).second.second.first << "," << (*it).second.second.second << "]" << std::endl;
						    //BUG_TRACK  }
						    //aafTreeDistPath[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]] = aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1-1][j1-1];
				                    //aafTreeDistPath[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]].push_back(make_pair(make_pair(std::pair<CHelixIndexTreeNode::CTreeNodePointer,CHelixIndexTreeNode::CTreeNodePointer>(*tIter1,*tIter2), std::pair<int,int>(i1-1,j1-1)), make_pair('R',std::pair<int,int>(i1,j1))));
						}
						else
						{
						    std::cout << "debug case 1" << std::endl;
						}
						// treedist(i1,j1)=forestdist(T1[l(i1)..i1],T2[l(j)..j1]) // put in permanent array
                                                // ......  
						
						// Routine 1: affTreeDistPath, Routine 2: tGlobalDistLookupTable
						//BUG_TRACK  tGlobalDistLookupTable[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]] = make_pair(std::pair<int, int>(tNodeIndex1[(*tIter1)], tNodeIndex2[(*tIter2)]), std::pair<int,int>(i1,j1));
						//std::cout << "save in permanent array: aafTreeDist["<<tNodeIndex1[(*tIter3)]<<"]["<<tNodeIndex2[(*tIter4)]<<"] =" << aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1] << std::endl;
						//std::cout << "save in affTreeDistPath:";
						/*
						std::list< std::pair<std::pair<std::pair<CHelixIndexTreeNode::CTreeNodePointer, CHelixIndexTreeNode::CTreeNodePointer>, std::pair<int,int> >, std::pair<char, std::pair<int, int> > >  >::iterator it;
						for ( it=aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].begin() ; it != aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].end(); it++ )
						{
						    std::cout << "[" << ((*it).first.first.first)->GetObject().ToString() << "," << ((*it).first.first.second)->GetObject().ToString() << "," << (*it).first.second.first << "," << (*it).first.second.second << "]" << "=" << (*it).second.first << "=>" << "[" << (*it).second.second.first << "," << (*it).second.second.second << "]"  ;
						}
						std::cout << std::endl;
						*/
					}
					else
					{
					        // not in the same corresponding subtree
						// forestdist(T1[l(i)..i1],T2[l(j)..j1]) = min {
						//    forestdist(T1[l(i)..i1-1],T2[l(j)..j1])+@(T1[i1]->^),
						//    forestdist(T1[l(i)..i1],T2[l(j)..j1-1])+@(^->T2[j1]),
						//    forestdist(T1[l(i)..l(i1)-1],T2[l(j)..l(j1)-1]) + treedist(i1,j1)}
						if (tNodeIndex1[(*tIter3)->GetLeftLeafPointer()]-1 < tNodeIndex1[(*tIter1)->GetLeftLeafPointer()])
						{
							iForestDistIndex1 = 0;
						}
						else
						{
							iForestDistIndex1 = tLocalTreeNodeIndex1[(*tIter3)->GetLeftLeafPointer()];
						}
						if (tNodeIndex2[(*tIter4)->GetLeftLeafPointer()]-1 < tNodeIndex2[(*tIter2)->GetLeftLeafPointer()])
						{
							iForestDistIndex2 = 0;
						}
						else
						{
							iForestDistIndex2 = tLocalTreeNodeIndex2[(*tIter4)->GetLeftLeafPointer()];
						}
						                                                                                                                                        /* tree_dist, only using affTreeDist, did not change affTreeDist, so we shouldn't change affTreeDistPath */
						fOtherNodeDistance = aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][iForestDistIndex1][iForestDistIndex2] + aafTreeDist[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]];
/*
C:6,8,0,3,62.5b,94.5h,(46)
C:6,8,0,3,62.5b,94.5i,(51)
C:6,8,0,3,62.5b,93b,(56)
C:6,8,0,3,62.5b,92.5b,(61)
C:6,8,0,3,62.5b,93b,(66)*/
						//BUG_TRACK  std::cout <<"C2:["<< (*tIter1)->GetObject().ToString()<<","<<(*tIter2)->GetObject().ToString()<<","<<i1<<","<<j1 <<"]("<< (*tIter3)->GetObject().ToString() <<","<< (*tIter4)->GetObject().ToString() <<")=("<< fDeletion <<","<< fInsertion <<","<< fOtherNodeDistance <<")from " << "aafForestDist["<<tNodeIndex1[(*tIter1)]<<"]["<<tNodeIndex2[(*tIter2)]<<"]["<<iForestDistIndex1<<"]["<<iForestDistIndex2<<"]("<<aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][iForestDistIndex1][iForestDistIndex2]<<")+"<< aafTreeDist[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]] << std::endl;  // 6,8,0,3=110 + 
						aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1] = GetMin(fDeletion, GetMin(fInsertion, fOtherNodeDistance));
						
						
						if (aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1] == fDeletion) {
						    //tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].push_back(make_pair(std::pair<int, int>(i1-1, j1), 'D'));
						    tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].push_back(make_pair(make_pair(std::pair<int, int>(i1-1, j1), 'D'), std::pair<int,int>(0,0)));
						    aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1] = aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1-1][j1];
				                    aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].push_back(make_pair(make_pair(std::pair<CHelixIndexTreeNode::CTreeNodePointer,CHelixIndexTreeNode::CTreeNodePointer>(*tIter1,*tIter2), std::pair<int,int>(i1-1,j1)), make_pair('D',std::pair<int,int>(i1,j1))));
						    //IMPORTANCE:It is not updated!!!! It wasted me 4 hours!!!! aafTreeDistPath[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]] = aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1];
						}
						else if (aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1] == fInsertion) {  // only consider the first element no matter there is more than one 
						    //tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].push_back(make_pair(std::pair<int, int>(i1, j1-1), 'I'));
						    tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].push_back(make_pair(make_pair(std::pair<int, int>(i1, j1-1), 'I'), std::pair<int,int>(0,0)));
						    aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1] = aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1-1];
				                    aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].push_back(make_pair(make_pair(std::pair<CHelixIndexTreeNode::CTreeNodePointer,CHelixIndexTreeNode::CTreeNodePointer>(*tIter1,*tIter2), std::pair<int,int>(i1,j1-1)), make_pair('I',std::pair<int,int>(i1,j1))));
						    //IMPORTANCE:It is not updated!!!! It wasted me 4 hours!!!! aafTreeDistPath[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]] = aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1];  // TODO: reserve calculation order, at first calculated TreeDistPath, then ForestDistPath = TreeDistPath, as the next case
						}
						else if (aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1] == fOtherNodeDistance) {
						    tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1].push_back(make_pair(make_pair(std::pair<int, int>(iForestDistIndex1, iForestDistIndex2), 'M'), std::pair<int,int>(tNodeIndex1[(*tIter3)], tNodeIndex2[(*tIter4)])));
						
						    
			                            std::list< std::pair<std::pair<std::pair<CHelixIndexTreeNode::CTreeNodePointer, CHelixIndexTreeNode::CTreeNodePointer>, std::pair<int,int> >, std::pair<char, std::pair<int, int> > >  > first = aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][iForestDistIndex1][iForestDistIndex2];  
					            std::list< std::pair<std::pair<std::pair<CHelixIndexTreeNode::CTreeNodePointer, CHelixIndexTreeNode::CTreeNodePointer>, std::pair<int,int> >, std::pair<char, std::pair<int, int> > >  >::iterator it;
                                                    
						    for ( it=aafTreeDistPath[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]].begin() ; it != aafTreeDistPath[tNodeIndex1[(*tIter3)]][tNodeIndex2[(*tIter4)]].end(); it++ )
						    {
							first.push_back(*it);
							//BUG_TRACK  std::cout << "B(tNodeIndex1[(*tIter3)]="<<tNodeIndex1[(*tIter3)]<<","<<"tNodeIndex2[(*tIter4)]="<<tNodeIndex2[(*tIter4)]<<"):[" << ((*it).first.first.first)->GetObject().ToString() << "," << ((*it).first.first.second)->GetObject().ToString() << "," << (*it).first.second.first << "," << (*it).first.second.second << "]" << "=" << (*it).second.first << "=>" << "[" << (*it).second.second.first << "," << (*it).second.second.second << "]" << std::endl;
						    }

						    aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i1][j1] = first;  //==first
						}
						else
						{
						    std::cout << "debug case 2" << std::endl;
						}
					}
				}
			}
			// #### here ends treedist(i,j) computation ####
			// dumps current forestdist matrix
	/*	std::cout << "print loopup table:" << std::endl;
		std::cout << "[";
		for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tSuffixIter2 = tSuffixNodesList2.begin(); tSuffixIter2 != tSuffixNodesList2.end(); tSuffixIter2++)
		{
		        std::cout << tLocalTreeNodeIndex2[(*tSuffixIter2)->GetLeftLeafPointer()];
		}
		std::cout << "]" << std::endl;*/
		
	
		}
	}

	fDistance = aafTreeDist[iTree1NodeCount-1][iTree2NodeCount-1];
	
	
	
	
	//###################################################
	//######### check all materials for backtracking ####
	//###################################################
/*  //BUG_TRACK  
	for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter1 = tGetLRKeyRoots1.begin(); tIter1 != tGetLRKeyRoots1.end(); tIter1++)
	{
	        // GetLRKeyRoots for tree 2 returned 1 nodes: [[0,N]]
		for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter2 = tGetLRKeyRoots2.begin(); tIter2 != tGetLRKeyRoots2.end(); tIter2++)
		{
			std::cout << "forest_dist :" << (*tIter1)->GetObject().ToString() << "," << (*tIter2)->GetObject().ToString() << std::endl;
			for (int i=0; i<=iTree1NodeCount; i++)  // iTableHeight
			{
				for (int j=0; j<=iTree2NodeCount; j++)  // iTableWidth
				{
					std::cout << aafForestDist[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i][j] << " ";
				}
				std::cout << std::endl;
			}
			//tOutputStream << std::endl;
			std::cout << std::endl;
		}
	}
	
	


	for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter1 = tGetLRKeyRoots1.begin(); tIter1 != tGetLRKeyRoots1.end(); tIter1++)
	{
	        // GetLRKeyRoots for tree 2 returned 1 nodes: [[0,N]]
		for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter2 = tGetLRKeyRoots2.begin(); tIter2 != tGetLRKeyRoots2.end(); tIter2++)
		{
			std::cout << "backtracing_table :" << (*tIter1)->GetObject().ToString() << "," << (*tIter2)->GetObject().ToString() << std::endl;
			for (int i=0; i<=iTree1NodeCount; i++)  // iTableHeight
			{
				for (int j=0; j<=iTree2NodeCount; j++)  // iTableWidth
				{ 
				  	//std::cout << "[" << tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][j][i].size() << "]";
				        
				        std::list<std::pair<std::pair<std::pair<int, int>, char>, std::pair<int, int> > >::iterator it;
				        //for (int k=0; k<tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i][j].size(); k++)
					for ( it=tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i][j].begin() ; it != tBacktrackingTables[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i][j].end(); it++ )
				        {
					    // TODO: store some link information into a new multimap = new multimap;
				            std::cout << "[" << (*it).first.first.first << "," << (*it).first.first.second << "]" << "<=" << (*it).first.second << "=" <<"["<<i<<","<<j<<"]";
				        }
				        std::cout << " ";
				}
				std::cout << std::endl;
			}
			//tOutputStream << std::endl;
			std::cout << std::endl;
		}
	}
	
	
	for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter1 = tGetLRKeyRoots1.begin(); tIter1 != tGetLRKeyRoots1.end(); tIter1++)
	{
	        // GetLRKeyRoots for tree 2 returned 1 nodes: [[0,N]]
		for (CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter2 = tGetLRKeyRoots2.begin(); tIter2 != tGetLRKeyRoots2.end(); tIter2++)
		{
			std::cout << "forest_dist_path :" << (*tIter1)->GetObject().ToString() << "," << (*tIter2)->GetObject().ToString() << std::endl;
			for (int i=0; i<=iTree1NodeCount; i++)  // iTableHeight
			{
				for (int j=0; j<=iTree2NodeCount; j++)  // iTableWidth
				{         
				        std::list< std::pair<std::pair<std::pair<CHelixIndexTreeNode::CTreeNodePointer, CHelixIndexTreeNode::CTreeNodePointer>, std::pair<int,int> >, std::pair<char, std::pair<int, int> > >  >::iterator it;
					for ( it=aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i][j].begin() ; it != aafForestDistPath[tNodeIndex1[(*tIter1)]][tNodeIndex2[(*tIter2)]][i][j].end(); it++ )
				        {
				            std::cout << "[" << ((*it).first.first.first)->GetObject().ToString() << "," << ((*it).first.first.second)->GetObject().ToString() << "," << (*it).first.second.first << "," << (*it).first.second.second << "]" << "=" << (*it).second.first << "=>" << "[" << (*it).second.second.first << "," << (*it).second.second.second << "]"  ;
				        }
				        std::cout << std::endl;
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
	}
	
	
	// dumps current treedist matrix
	std::cout << "tree_dist" << std::endl;
	for (int i=0; i<iTree1NodeCount; i++)
	{
		for (int j=0; j<iTree2NodeCount; j++)
		{
			std::cout << aafTreeDist[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	
	std::cout << "tree_dist_path" << std::endl;
	for (int i=0; i<iTree1NodeCount; i++)
	{
		for (int j=0; j<iTree2NodeCount; j++)
		{
			std::list< std::pair<std::pair<std::pair<CHelixIndexTreeNode::CTreeNodePointer, CHelixIndexTreeNode::CTreeNodePointer>, std::pair<int,int> >, std::pair<char, std::pair<int, int> > >  >::iterator it;
			for ( it=aafTreeDistPath[i][j].begin() ; it != aafTreeDistPath[i][j].end(); it++ )
			{
			    std::cout << "[" << ((*it).first.first.first)->GetObject().ToString() << "," << ((*it).first.first.second)->GetObject().ToString() << "," << (*it).first.second.first << "," << (*it).first.second.second << "]" << "=" << (*it).second.first << "=>" << "[" << (*it).second.second.first << "," << (*it).second.second.second << "]"  ;
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	
	std::cout << "tree_dist_lookup_table" << std::endl;
	for (int i=0; i<iTree1NodeCount; i++)
	{
		for (int j=0; j<iTree2NodeCount; j++)
		{
			//std::cout << (*(tGlobalDistLookupTable[i][j].first.first)).GetObject().ToString() << ""
			//          << (*(tGlobalDistLookupTable[i][j].first.second)).GetObject().ToString() << " "
			std::cout << tGlobalDistLookupTable[i][j].first.first << " " \
			          << tGlobalDistLookupTable[i][j].first.second << " " \
			          << tGlobalDistLookupTable[i][j].second.first << " " \
			          << tGlobalDistLookupTable[i][j].second.second << "    ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
*/
	
	
/*	
       TODO: print the tree dist path
             why the last one not complete?
             represent the format
             check the hishape in the alignment have the same order as inputs 
             load a scheme score file from current directory
             integrate the Probabality calculation into program as option -p
*/	

	// backtrackingVector is basis ====using containers====> expression technique
	// ##################### express the backtrackingVector using containers #######################
	// 表达出 这个链: 如果遇到 M， recursive digging
	std::vector<std::string> vecAlignmentRow1, vecAlignmentRow2;
	std::vector<char> vecAlignmentType;

	std::list< std::pair<std::pair<std::pair<CHelixIndexTreeNode::CTreeNodePointer, CHelixIndexTreeNode::CTreeNodePointer>, std::pair<int,int> >, std::pair<char, std::pair<int, int> > >  >::iterator it;
	for ( it=aafForestDistPath[tNodeIndex1[ptTree1]][tNodeIndex2[ptTree2]][ptTree1->GetChildrenCount()+1][ptTree2->GetChildrenCount()+1].begin() ; it != aafForestDistPath[tNodeIndex1[ptTree1]][tNodeIndex2[ptTree2]][ptTree1->GetChildrenCount()+1][ptTree2->GetChildrenCount()+1].end(); it++ )
	{
	    //BUG_TRACK  std::cout << "[" << ((*it).first.first.first)->GetObject().ToString() << "," << ((*it).first.first.second)->GetObject().ToString() << "," << (*it).first.second.first << "," << (*it).first.second.second << "]" << "=" << (*it).second.first << "=>" << "[" << (*it).second.second.first << "," << (*it).second.second.second << "]";
	    
	    CHelixIndexTreeNode::CTreeNodePointerList tempSuffixNodesList1 = ((*it).first.first.first)->GetSuffixList();
	    for (int i=0; i<(*it).second.second.first-1; i++)
	    {
	        tempSuffixNodesList1.pop_front();
	    }
	    //BUG_TRACK  std::cout << tempSuffixNodesList1.front()->GetObject().ToString() << std::endl;
	    
	    CHelixIndexTreeNode::CTreeNodePointerList tempSuffixNodesList2 = ((*it).first.first.second)->GetSuffixList();
	    for (int i=0; i<(*it).second.second.second-1; i++)
	    {
		tempSuffixNodesList2.pop_front();
	    }
	    //BUG_TRACK  std::cout << tempSuffixNodesList2.front()->GetObject().ToString() << std::endl;
	    
	    if ((*it).second.first == 'R')
	    {
	        std::string tempString1 = tempSuffixNodesList1.front()->GetObject().ToString();
		boost::replace_all(tempString1, "-1m", "m");
		std::string tempString2 = tempSuffixNodesList2.front()->GetObject().ToString();
		boost::replace_all(tempString2, "-1m", "m");
		vecAlignmentRow1.push_back(tempString1);
		vecAlignmentType.push_back('R');
		vecAlignmentRow2.push_back(tempString2);
	    }
	    else if ((*it).second.first == 'D')
	    {
	        std::string tempString1 = tempSuffixNodesList1.front()->GetObject().ToString();
		boost::replace_all(tempString1, "-1m", "m");
		vecAlignmentRow1.push_back(tempString1);
		vecAlignmentType.push_back('D');
		vecAlignmentRow2.push_back("-");
	    }
	    else if ((*it).second.first == 'I')
	    {
	        std::string tempString2 = tempSuffixNodesList2.front()->GetObject().ToString();
		boost::replace_all(tempString2, "-1m", "m");
		vecAlignmentRow1.push_back("-");
		vecAlignmentType.push_back('I');
		vecAlignmentRow2.push_back(tempString2);
	    }
	}
	//BUG_TRACK  std::cout << std::endl;



	// present the data structure
	std::ostringstream strAlignmentRow1, strAlignmentType, strAlignmentRow2;  
	for (int i=0; i<(int)(vecAlignmentType.size()); i++)
	{
	    if (vecAlignmentType[i] == 'R')
	    {
		int iHelixIndex1Length = vecAlignmentRow1[i].length();
		int iHelixIndex2Length = vecAlignmentRow2[i].length();
		std::string strMinSpaces (4,' ');   
		if (iHelixIndex1Length > iHelixIndex2Length)
		{
		    std::string strLongerSpaces ((4+iHelixIndex1Length-iHelixIndex2Length), ' ');
		    std::string strTypeSpaces ((4+iHelixIndex1Length-1), ' ');
		    strAlignmentRow1 << strMinSpaces << vecAlignmentRow1[i];
		    strAlignmentType << strTypeSpaces << '|';
		    strAlignmentRow2 << strLongerSpaces << vecAlignmentRow2[i];
		}
		else
		{
		    std::string strLongerSpaces ((4+iHelixIndex2Length-iHelixIndex1Length), ' ');
		    std::string strTypeSpaces ((4+iHelixIndex2Length-1), ' ');
		    strAlignmentRow1 << strLongerSpaces << vecAlignmentRow1[i];
		    strAlignmentType << strTypeSpaces << '|';
		    strAlignmentRow2 << strMinSpaces << vecAlignmentRow2[i];
		}
	    }
	    else if (vecAlignmentType[i] == 'D')
	    {
		int iHelixIndex1Length = vecAlignmentRow1[i].length();
		std::string strMinSpaces (4,' ');   
		std::string strLongerSpaces ((4+iHelixIndex1Length-1), ' ');
		std::string strTypeSpaces ((4+iHelixIndex1Length), ' ');
		strAlignmentRow1 << strMinSpaces << vecAlignmentRow1[i];
		strAlignmentType << strTypeSpaces;
		strAlignmentRow2 << strLongerSpaces << '-';
	    }
	    else if (vecAlignmentType[i] == 'I')
	    {
		int iHelixIndex2Length = vecAlignmentRow2[i].length();
		std::string strMinSpaces (4,' ');   
		std::string strLongerSpaces ((4+iHelixIndex2Length-1), ' ');
		std::string strTypeSpaces ((4+iHelixIndex2Length), ' ');
		strAlignmentRow1 << strLongerSpaces << '-';
		strAlignmentType << strTypeSpaces;
		strAlignmentRow2 << strMinSpaces << vecAlignmentRow2[i];
	    }
	}
	std::cout << std::endl;
	std::cout << strAlignmentRow1.str() << std::endl;
	std::cout << strAlignmentType.str() << std::endl;
	std::cout << strAlignmentRow2.str() << std::endl;

	
	//###################################################
	//######### free resources #############
	//###################################################
	
	// due to the trees were constructed in own data structure, they have to be freed after using. 
	// in compare with this, every data contained in c++ container, for example vector. they don't need to be freed.
	
	// frees forestdist table memory
	/*for (int i=0; i<=iTree1NodeCount; i++)
	{
		delete [] aafForestDist[i];
	}
	delete [] aafForestDist;*/

	// frees treedist table memory
	/*for (int i=0; i<iTree1NodeCount; i++)
	{
		delete [] aafTreeDist[i];
	}
	delete [] aafTreeDist;*/

	// free trees
	delete ptTree1;
	delete ptTree2;

	m_strOutput = tOutputStream.str();
	m_fLastScore = fDistance;

	
	return fDistance;
}


}; // namespace algorithms
}; // namespace comp
}; // namespace hishape
