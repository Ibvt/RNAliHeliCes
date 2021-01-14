//*****************************************************************************
#ifndef __SHAPIROALIGNER_H
#define __SHAPIROALIGNER_H 1


/**
@file ShapiroAligner.h
****************************
@brief ShapiroAligner.h is the header file for ZS89 algorithm.

ShapiroAligner.h is the header file of the implementation of
Zhang & Shasha '89 algorithm for trees alignment.

@see hishape::comp::algorithms

*/
#include "PairAligner.hpp"
//#include "Aligner.hpp"    // SSAligner -> Alinger
//#include "Hishape.hpp"    // from SSAligner
//#include "ScoreScheme.hpp"// from SSAlinger
#include "TreeNode.hpp"
#include "HelixIndex.hpp"

namespace hishape
{
namespace comp
{
namespace algorithms
{
/**
@class CZhangShasha89Aligner
****************************
@brief Zhang & Shasha '89 algorithm for trees alignment implementation.

@c CZhangShasha89Aligner is the implementation of:\n
"Simple fast algorithms for the editing distance between trees and related
problems", K. Zhang, D. Shasha, Society for Industrial and Applied Mathematics
Journal Volume 18, Philadelphia, PA, USA, ISSN:0097-5397, Issue 6
(December 1989).
This algorithm can align 2 secondary structure seen as trees.

*/
class CShapiroAligner
: public CPairAligner
{
public:
	//! A tree node that contains a pair of nucleotide type.
	//typedef hishape::comp::tools::CTreeNode<hishape::comp::CNucleotidePair> CNucleotidePairTreeNode;
	typedef val::biocpp::tools::CTreeNode<hishape::comp::CHelixIndex> CHelixIndexTreeNode;

private:
	/*
	@var SNodePtrCompare
	********************
	Functor used to map tree nodes pointers in a map.

	@param CHelixIndexTreeNode::CTreeNodePointer ptNode1: pointer to node 1

	@param CHelixIndexTreeNode::CTreeNodePointer ptNode2: pointer to node 2

	@return bool: true if pointer 1 is before pointer 2.
	*/
	struct SNodePtrCompare
	{
		bool operator()(CHelixIndexTreeNode::CTreeNodePointer ptNode1, CHelixIndexTreeNode::CTreeNodePointer ptNode2) const
		{
			return (ptNode1 < ptNode2);
		}
	};
	
	struct IntegerCompare
	{
		bool operator()(int iNodeIndex1, int iNodeIndex2) const
		{
			return (iNodeIndex1 < iNodeIndex2);
		}
	};
	//! algorithm init status
	bool m_bInited;
	//! Returns the list of key roots
	CHelixIndexTreeNode::CTreeNodePointerList GetLRKeyRoots(CHelixIndexTreeNode::CTreeNodePointer ptNode) const;
	
	//! Returns the minimum value between two floats
	const float& GetMin(const float& fFirstValue, const float& fSecondValue);

	//! Returns the maximum value between two floats
	const float& GetMax(const float& fFirstValue, const float& fSecondValue);
	
public:
    /*    //! The first structure.
	CHishape m_tFirst;     // from SSAligner
	//! The second structure.
	CHishape m_tSecond;    // from SSAligner
*/	
	
	//! Standard constructor.
	CShapiroAligner(void);
	//CShapiroAligner(const CHishape &hishape1, CHishape &hishape2);
	CShapiroAligner(const CHelixIndexTreeNode::CTreeNodePointer &ptTree1, const CHelixIndexTreeNode::CTreeNodePointer &ptTree2, float ratio/*, int hishape_type*/);
	//! Non-virtual destructor.
	~CShapiroAligner(void);
	//! The main alignement function.
	float Align(bool bStoreDistMatrix=true, bool bComputeBacktracking=true);
	//! Cycles estimator.
	//int EstimateCycles(void) const;
	
protected:
	//! The first structure.
	//CHishape m_hiShape1;
	//! The second structure.
	//CHishape m_hiShape2;
	//! The alignment function without backtracking.
	float AlignNoBacktracking(bool bStoreScoreMatrix=true);
	//! The alignment function with backtracking.
	float AlignWithBacktracking(bool bStoreScoreMatrix=true);
	
};
}; // namespace algorithms
}; // namespace comp
}; // namespace hishape
#endif //ifndef __SHAPIROALIGNER_H
