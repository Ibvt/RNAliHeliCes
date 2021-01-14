//*****************************************************************************
#ifndef __PAIRALIGNER_H
#define __PAIRALIGNER_H 1


/**
@file PairAligner.h
*******************
@brief PairAligner.h is the class header for
@c hishape::comp::algorithms::CPairAligner class.

PairAligner.h is the class header for @c hishape::comp::algorithms::CPairAligner.
It defines how a nucleic acid secondary structure (NASS) aligner should behave.
This class can't be instanciated and should be used as an interface for
secondary structure alignment algorithms.

@see hishape::comp::algorithms
*/
#include "Aligner.hpp"
//#include "Hishape.hpp"
#include "TreeNode.hpp"
#include "HelixIndex.hpp"

#include "ScoreScheme.hpp"

namespace hishape
{
namespace comp
{
namespace algorithms
{

/**
@name Similarity Offsets
************************
The following constants contain the offsets of similarity rates stored in
class member @c m_afLastSimilarityRate.

*/
//@{
//! Structure similarity offset.
const unsigned int I_NASSALIGNER_STRUCT_SIM_OFFSET  = 0;
//! Bases similarity offset.
const unsigned int I_NASSALIGNER_BASES_SIM_OFFSET   = 1;
//! Overall similarity offset.
const unsigned int I_NASSALIGNER_OVERALL_SIM_OFFSET = 2;
//@}


typedef val::biocpp::tools::CTreeNode<hishape::comp::CHelixIndex> CHelixIndexTreeNode;
/**
@class CPairAligner   PairAligner.h "PairAligner.hpp"
*******************
@brief CPairAligner is an interface for secondary structure alignement.

@c CPairAligner is an interface that describes the behaviour of a nucleic acid
secondary structure alignement algorithm.
*/
class CPairAligner
: public val::biocpp::algorithms::CAligner
{
public:
	//! Standard constructor.
	CPairAligner(void);
	//! Copy constructor.
	CPairAligner(const CPairAligner &tAligner);
	// Constructor with a score scheme.
	CPairAligner(const CScoreScheme &tScoreScheme);
	//! Virtual destructor.
	virtual ~CPairAligner(void);
	//! Sets the first tree to work on.
	virtual void SetFirstTree(const CHelixIndexTreeNode::CTreeNodePointer &tFirst);
	//! Sets the second tree to align with the first tree.
	virtual void SetSecondTree(const CHelixIndexTreeNode::CTreeNodePointer &tSecond);
	//! Sets the score scheme to use.
	virtual void SetScoreScheme(const hishape::comp::CScoreScheme &tScoreScheme);
	//! Returns current score scheme.
	virtual const CScoreScheme &GetScoreScheme() const;
	//! Returns the last similarity rate computed by @c Align.
	virtual float GetStructureSimilarityRate() const;
	//! Returns the last similarity rate computed by @c Align.
	virtual float GetBasesSimilarityRate() const;
	//! Returns the last similarity rate computed by @c Align.
	virtual float GetSimilarityRate() const;
	//! Assignation operator.
	virtual CPairAligner& operator=(const CPairAligner &tAligner);
	
	
	//! The first structure.
	CHelixIndexTreeNode::CTreeNodePointer m_tFirst;
	//! The second structure.
	CHelixIndexTreeNode::CTreeNodePointer m_tSecond;
	
	float m_ratio;
	/*int m_hishapeType;*/
protected:

	//! NA score scheme
	hishape::comp::CScoreScheme *m_ptNAScoreScheme;
	//! Last similarity rate
	float m_afLastSimilarityRate[3];
	//! Returns sequences splitting cost
	//virtual float GetSequenceSplitCost() const;
	//! Returns sequences joining cost
	//virtual float GetSequencesJoinCost() const;
	//! Returns pair-breaking cost
	//virtual float GetPairBreakingCost() const;
	//! Returns pair-breaking cost
	//virtual float GetPairBreakingCost(const hishape::comp::CHelixIndex &tPair) const;
	//! Returns pair-creation cost
	//virtual float GetPairCreationCost() const;
	//! Returns pair-creation cost
	//virtual float GetPairCreationCost(const hishape::comp::CHelixIndex &tPair) const;
	

};

}; // namespace algorithms
}; // namespace comp
}; // namespace hishape
#endif //ifndef __NASSALIGNER_H
