//*****************************************************************************
#include "PairAligner.hpp"
#include "ScoreScheme.hpp"
#include "DebugTools.hpp"


namespace hishape
{
namespace comp
{
namespace algorithms
{

//+FIXME: faire un meilleur design pour l'utilisation des score schemes.
/**
Standard constructor
********************
A standard constructor that sets the initial score scheme.
*/
CPairAligner::CPairAligner(void)
{
	TRACEI("CPairAligner::CPairAligner(void)");
	m_ptScoreScheme = m_ptNAScoreScheme = new CScoreScheme();
	//m_ptNAScoreScheme->InitAsStandardScoreScheme();  // it doesn't need this method
	m_afLastSimilarityRate[I_NASSALIGNER_STRUCT_SIM_OFFSET] = 0;
	m_afLastSimilarityRate[I_NASSALIGNER_BASES_SIM_OFFSET] = 0;
	m_afLastSimilarityRate[I_NASSALIGNER_OVERALL_SIM_OFFSET] = 0;
}


/**
Copy constructor
****************
*/
CPairAligner::CPairAligner(const CPairAligner &tAligner)
: CAligner()
{
	TRACEI("CPairAligner::CPairAligner(const CPairAligner &tAligner)");
	m_afLastSimilarityRate[I_NASSALIGNER_STRUCT_SIM_OFFSET] = tAligner.m_afLastSimilarityRate[I_NASSALIGNER_STRUCT_SIM_OFFSET];
	m_afLastSimilarityRate[I_NASSALIGNER_BASES_SIM_OFFSET] = tAligner.m_afLastSimilarityRate[I_NASSALIGNER_BASES_SIM_OFFSET];
	m_afLastSimilarityRate[I_NASSALIGNER_OVERALL_SIM_OFFSET] = tAligner.m_afLastSimilarityRate[I_NASSALIGNER_OVERALL_SIM_OFFSET];
}


/**
Constructor with score scheme
*****************************
Constructor with a score scheme provided.
*/
CPairAligner::CPairAligner(const CScoreScheme &tScoreScheme)
: CAligner()
{
	TRACEI("CAligner::CAligner(const CScoreScheme &tScoreScheme)");
	m_ptScoreScheme = m_ptNAScoreScheme = new CScoreScheme(tScoreScheme);
	m_afLastSimilarityRate[I_NASSALIGNER_STRUCT_SIM_OFFSET] = 0;
	m_afLastSimilarityRate[I_NASSALIGNER_BASES_SIM_OFFSET] = 0;
	m_afLastSimilarityRate[I_NASSALIGNER_OVERALL_SIM_OFFSET] = 0;
}


/**
Destructor
**********
A simple virtual destructor.
*/
CPairAligner::~CPairAligner(void)
{
	TRACEI("CPairAligner::~CPairAligner(void)");
}


/**
@fn void CPairAligner::SetFirstTree(const CHelixIndexTreeNode::CTreeNodePointer &tFirst)
****************************************************************************
Selects the first tree as the reference tree.

@param const CHelixIndexTreeNode::CTreeNodePointer &tFirst:\n
 the first tree.
*/
void CPairAligner::SetFirstTree(const CHelixIndexTreeNode::CTreeNodePointer &tFirst)
{
	TRACEI("void CPairAligner::SetFirstStructure(const CNASecondaryStructure &tFirst)");
	//int i;
	//for (i=0; i < tFirst.GetLength(); i++) {
	//    std::cout << "PairAligner.cpp::tFirst(" << i << ")=" << tFirst.GetPairedBaseIndex(i) << std::endl;
	//}
	// operator= was not implemented
	// TODO: improve it with deep copy
	m_tFirst = tFirst;
}


/**
@fn void CPairAligner::SetSecondTree(const CHelixIndexTreeNode::CTreeNodePointer &tSecond)
******************************************************************************
Selects the second tree to align with the first one.

@param const CHelixIndexTreeNode::CTreeNodePointer &tSecond:\n
 the second tree.                     
*/
void CPairAligner::SetSecondTree(const CHelixIndexTreeNode::CTreeNodePointer &tSecond)
{
        // TODO: improve it with deep copy
	m_tSecond = tSecond;
}


/**
@fn void CAligner::SetScoreScheme(const CScoreScheme &tScoreScheme)
*******************************************************************
Sets the score scheme to use for alignements. This methode is the only one
always called each time the score scheme needs to be changed.

@param const CScoreScheme &tScoreScheme:
 a new score scheme to use.
*/
/*
void CPairAligner::SetScoreScheme(const hishape::comp::CScoreScheme &tScoreScheme)
{
	TRACEI("void CPairAligner::SetScoreScheme(const hishape::comp::CScoreScheme &tScoreScheme)");
	delete m_ptScoreScheme;
	m_ptScoreScheme = m_ptNAScoreScheme = new CScoreScheme();
	m_ptNAScoreScheme->InitAsStandardScoreScheme();
	(*m_ptScoreScheme) = tScoreScheme;
}*/


/**
@fn void CAligner::SetScoreScheme(const CScoreScheme &tScoreScheme)
*******************************************************************
Sets the score scheme to use for alignements. This methode is the only one
always called each time the score scheme needs to be changed.

@param const CScoreScheme &tScoreScheme:
 a new score scheme to use.
*/
void CPairAligner::SetScoreScheme(const hishape::comp::CScoreScheme &tScoreScheme)
{
	TRACEI("void CPairAligner::SetScoreScheme(const CScoreScheme &tScoreScheme)");
	delete m_ptScoreScheme;
	m_ptScoreScheme = m_ptNAScoreScheme = new CScoreScheme(tScoreScheme);
}


/**
@fn const CScoreScheme &tScoreScheme GetScoreScheme() const
*****************************************************************
Returns the score scheme used by current alignment algorithm.

@return const CScoreScheme &:
 a reference to current score scheme.
*/
const CScoreScheme &CPairAligner::GetScoreScheme() const
{
	TRACEI("const CScoreScheme &CPairAligner::GetScoreScheme() const");
	return *m_ptNAScoreScheme;
}


/**
@fn float CPairAligner::GetSequenceSplitCost() const
****************************************************
Returns the operation cost of splitting a sequence into 2 sequences.

@return float:
 the cost of splitting a sequence.
*/
/*
float CPairAligner::GetSequenceSplitCost() const
{
	return m_ptNAScoreScheme->GetSequenceSplitCost();
}*/


/**
@fn float CPairAligner::GetSequencesJoinCost() const
****************************************************
Returns the operation cost of joining 2 sequences into 1.

@return float:
 the cost of joining 2 sequences into 1.
*/
/*
float CPairAligner::GetSequencesJoinCost() const
{
	return m_ptNAScoreScheme->GetSequencesJoinCost();
}*/


/**
@fn float CPairAligner::GetPairBreakingCost() const
*************************************************
Returns the operation cost of breaking a pair of bases into two independent
bases.

@return float:
 the cost of breaking a pair.
*/
/*
float CPairAligner::GetPairBreakingCost() const
{
	return m_ptNAScoreScheme->GetPairBreakingCost();
}
*/


/**
@fn float CPairAligner::GetPairBreakingCost(const hishape::comp::CNucleotidePair &tPair) const
********************************************************************************************
Returns the operation cost of breaking given pair of bases into two independent
bases.

@param const hishape::comp::CNucleotidePair &tPair:
 Pair to break.

@return float:
 the cost of breaking given pair.
*/
/*
float CPairAligner::GetPairBreakingCost(const hishape::comp::CHelixIndex &tPair) const
{
	return m_ptNAScoreScheme->GetPairBreakingCost(tPair);
}
*/


/**
@fn float CPairAligner::GetPairCreationCost() const
***************************************************
Returns the operation cost of creating a pair of bases from two independent
bases.

@return float:
 the cost of creating a pair.
*/
/*
float CPairAligner::GetPairCreationCost() const
{
	return m_ptNAScoreScheme->GetPairCreationCost();
}
*/


/**
@fn float CPairAligner::GetPairCreationCost(const hishape::comp::CNucleotidePair &tPair) const
********************************************************************************************
Returns the operation cost of creating given pair of bases from two independent
bases.

@param const hishape::comp::CNucleotidePair &tPair:
 the pair to create.

@return float:
 the cost of creating given pair.
*/
/*
float CPairAligner::GetPairCreationCost(const hishape::comp::CHelixIndex &tPair) const
{
	return m_ptNAScoreScheme->GetPairCreationCost(tPair);
}
*/


/**
@fn float CPairAligner::GetStructureSimilarityRate()
************************************************
Returns the last structure similarity rate computed by @c Align.

@return float:
 last alignment structure similarity rate computed by the last @c Align call.
*/
float CPairAligner::GetStructureSimilarityRate() const
{
	TRACEI("float CPairAligner::GetStructureSimilarityRate() const");
	return m_afLastSimilarityRate[I_NASSALIGNER_STRUCT_SIM_OFFSET];
}


/**
@fn float CPairAligner::GetBasesSimilarityRate()
***********************************************
Returns the last bases similarity rate computed by @c Align.

@return float:
 last alignment bases similarity rate computed by the last @c Align call.
*/
float CPairAligner::GetBasesSimilarityRate() const
{
	TRACEI("float CPairAligner::GetBasesSimilarityRate() const");
	return m_afLastSimilarityRate[I_NASSALIGNER_BASES_SIM_OFFSET];
}


/**
@fn float CPairAligner::GetSimilarityRate()
***************************************
Returns the last similarity rate computed by @c Align.

@return float:
 last alignment similarity rate computed by the last @c Align call.
*/
float CPairAligner::GetSimilarityRate() const
{
	TRACEI("float CPairAligner::GetSimilarityRate() const");
	return m_afLastSimilarityRate[I_NASSALIGNER_OVERALL_SIM_OFFSET];
}


/**
@fn CPairAligner& CPairAligner::operator=(const CPairAligner &tAligner)
***********************************************************************
*/
CPairAligner& CPairAligner::operator=(const CPairAligner &tAligner)
{
	TRACEI("CPairAligner& CPairAligner::operator=(const CPairAligner &tAligner)");
	if (this != &tAligner)
	{
		CAligner::operator=(tAligner);
		m_afLastSimilarityRate[I_NASSALIGNER_STRUCT_SIM_OFFSET] = tAligner.m_afLastSimilarityRate[I_NASSALIGNER_STRUCT_SIM_OFFSET];
		m_afLastSimilarityRate[I_NASSALIGNER_BASES_SIM_OFFSET] = tAligner.m_afLastSimilarityRate[I_NASSALIGNER_BASES_SIM_OFFSET];
		m_afLastSimilarityRate[I_NASSALIGNER_OVERALL_SIM_OFFSET] = tAligner.m_afLastSimilarityRate[I_NASSALIGNER_OVERALL_SIM_OFFSET];
	}
	return *this;
}


}; // namespace algorithms
}; // namespace comp
}; // namespace hishape
