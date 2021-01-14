//*****************************************************************************
#include "NASecondaryStructure.hpp"
#include <algorithm>
#include <list>
#include <sstream>
#include <stdexcept>
#include <stack>

#include "DebugTools.hpp"


namespace val
{
namespace biocpp
{
/**
Default constructor
*******************
Returns a new instance of @c CNASecondaryStructure.

@return CNASecondaryStructure:\n
 a new instance of @c CNASecondaryStructure.
*/
CNASecondaryStructure::CNASecondaryStructure(void)
: CNAPrimaryStructure()
, m_bIsStructureParsed(false)
, m_iExternalBasesCount(0)
, m_iStackedPairBasesCount(0)
, m_iBulgeBasesCount(0)
, m_iInternalLoopBasesCount(0)
, m_iHairpinLoopBasesCount(0)
, m_iMultibrachedLoopBasesCount(0)
, m_iPseudoknotBasesCount(0)
, m_iKissingHairpinBasesCount(0)
, m_iKissingLoopBasesCount(0)
, m_iTripleHelixBasesCount(0)
, m_iStemSplitPosition(0)
, m_tNucleotidesPairs(1)
, m_tNucleotidesTertiaryInteractions(1)
, m_tNucleotidesConfig()
, m_tNucleotidesConfigTertiaryStructure()
, m_iStemsCount()
, m_iStemLoopsCount()
, m_tStemsAndStemLoopsNature()
, m_tStemsAndStemLoopsStartingBases()
{
	TRACEI("CNASecondaryStructure::CNASecondaryStructure(void)");
}


/**
copy constructor
****************
Returns a new instance of @c CNASecondaryStructure which is a copy of the
given instance.

@param const CNASecondaryStructure& tSecondaryStructure:\n
 a secondary structure to duplicate.

@return CNASecondaryStructure:\n
 a new instance of @c CNASecondaryStructure.
*/
CNASecondaryStructure::CNASecondaryStructure(const CNASecondaryStructure& tSecondaryStructure)
: CNAPrimaryStructure(tSecondaryStructure)
, m_bIsStructureParsed(tSecondaryStructure.m_bIsStructureParsed)
, m_iExternalBasesCount(tSecondaryStructure.m_iExternalBasesCount)
, m_iStackedPairBasesCount(tSecondaryStructure.m_iStackedPairBasesCount)
, m_iBulgeBasesCount(tSecondaryStructure.m_iBulgeBasesCount)
, m_iInternalLoopBasesCount(tSecondaryStructure.m_iInternalLoopBasesCount)
, m_iHairpinLoopBasesCount(tSecondaryStructure.m_iHairpinLoopBasesCount)
, m_iMultibrachedLoopBasesCount(tSecondaryStructure.m_iMultibrachedLoopBasesCount)
, m_iPseudoknotBasesCount(tSecondaryStructure.m_iPseudoknotBasesCount)
, m_iKissingHairpinBasesCount(tSecondaryStructure.m_iKissingHairpinBasesCount)
, m_iKissingLoopBasesCount(tSecondaryStructure.m_iKissingLoopBasesCount)
, m_iTripleHelixBasesCount(tSecondaryStructure.m_iTripleHelixBasesCount)
, m_iStemSplitPosition(tSecondaryStructure.m_iStemSplitPosition)
, m_tNucleotidesPairs(tSecondaryStructure.m_tNucleotidesPairs)
, m_tNucleotidesTertiaryInteractions(tSecondaryStructure.m_tNucleotidesTertiaryInteractions)
, m_tNucleotidesConfig(tSecondaryStructure.m_tNucleotidesConfig)
, m_tNucleotidesConfigTertiaryStructure(tSecondaryStructure.m_tNucleotidesConfigTertiaryStructure)
, m_iStemsCount(tSecondaryStructure.m_iStemsCount)
, m_iStemLoopsCount(tSecondaryStructure.m_iStemLoopsCount)
, m_tStemsAndStemLoopsNature(tSecondaryStructure.m_tStemsAndStemLoopsNature)
, m_tStemsAndStemLoopsStartingBases(tSecondaryStructure.m_tStemsAndStemLoopsStartingBases)
{
	TRACEI("CNASecondaryStructure::CNASecondaryStructure(const CNASecondaryStructure& tSecondaryStructure)");
}


/**
Destructor
**********
A standard virtual destructor.
*/
CNASecondaryStructure::~CNASecondaryStructure(void)
{
	TRACEI("CNASecondaryStructure::~CNASecondaryStructure(void)");
}


/*
1-based
*/
CNASecondaryStructure::CDequeIntIterator CNASecondaryStructure::GetPairsDataIteratorAt(int iAtPosition)
{
	TRACEI("CDequeIntIterator CNASecondaryStructure::GetPairsDataIteratorAt(int iAtPosition)");
	// keep in mind that iAtPosition is not a 0-based index
	// ie. first position is 1
	if ((iAtPosition <= 0) || ((signed) m_tNucleotidesPairs.size() < iAtPosition))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::GetPairsDataIteratorAt(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	return (m_tNucleotidesPairs.begin() + iAtPosition);
}


/*
1-based
*/
CNASecondaryStructure::CDequeIntIterator CNASecondaryStructure::GetTertiaryInteractionsDataIteratorAt(int iAtPosition)
{
	TRACEI("CDequeIntIterator CNASecondaryStructure::GetPairsDataIteratorAt(int iAtPosition)");
	// keep in mind that iAtPosition is not a 0-based index
	// ie. first position is 1
	if ((iAtPosition <= 0) || ((signed) m_tNucleotidesTertiaryInteractions.size() < iAtPosition))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::GetTertiaryInteractionsDataIteratorAt(" << iAtPosition << ": " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	return (m_tNucleotidesTertiaryInteractions.begin() + iAtPosition);
}


/**
@fn CNASecondaryStructure& CNASecondaryStructure::operator=(const CNASecondaryStructure &tSecondaryStructure)
*************************************************************************************************************
Assignation operator. Copies the given secondary structure (including links)
into this one.

@param const CNASecondaryStructure &tSecondaryStructure:\n
 the secondary structure to assign to the current one.

@return CNASecondaryStructure&:\n
 returns (*this) object.
*/
CNASecondaryStructure& CNASecondaryStructure::operator=(const CNASecondaryStructure &tSecondaryStructure)
{
	TRACEI("CNASecondaryStructure& CNASecondaryStructure::operator=(const CNASecondaryStructure &tSecondaryStructure)");

	// check if it's a self-assignation
	if (this == &tSecondaryStructure)
	{
		return *this;
	}
	// call inherited = operator
	CNAPrimaryStructure::operator=(tSecondaryStructure);
	// copy fields
	m_bIsStructureParsed = tSecondaryStructure.m_bIsStructureParsed;
	m_iExternalBasesCount = tSecondaryStructure.m_iExternalBasesCount;
	m_iStackedPairBasesCount = tSecondaryStructure.m_iStackedPairBasesCount;
	m_iBulgeBasesCount = tSecondaryStructure.m_iBulgeBasesCount;
	m_iInternalLoopBasesCount = tSecondaryStructure.m_iInternalLoopBasesCount;
	m_iHairpinLoopBasesCount = tSecondaryStructure.m_iHairpinLoopBasesCount;
	m_iMultibrachedLoopBasesCount = tSecondaryStructure.m_iMultibrachedLoopBasesCount;
	m_iPseudoknotBasesCount = tSecondaryStructure.m_iPseudoknotBasesCount;
	m_iKissingHairpinBasesCount = tSecondaryStructure.m_iKissingHairpinBasesCount;
	m_iKissingLoopBasesCount = tSecondaryStructure.m_iKissingLoopBasesCount;
	m_iTripleHelixBasesCount = tSecondaryStructure.m_iTripleHelixBasesCount;
	m_iStemSplitPosition = tSecondaryStructure.m_iStemSplitPosition;
	m_tNucleotidesPairs = tSecondaryStructure.m_tNucleotidesPairs;
	m_tNucleotidesTertiaryInteractions = tSecondaryStructure.m_tNucleotidesTertiaryInteractions;
	m_tNucleotidesConfig = tSecondaryStructure.m_tNucleotidesConfig;
	m_tNucleotidesConfigTertiaryStructure = tSecondaryStructure.m_tNucleotidesConfigTertiaryStructure;
	m_iStemsCount = tSecondaryStructure.m_iStemsCount;
	m_iStemLoopsCount = tSecondaryStructure.m_iStemLoopsCount;
	m_tStemsAndStemLoopsNature = tSecondaryStructure.m_tStemsAndStemLoopsNature;
	m_tStemsAndStemLoopsStartingBases = tSecondaryStructure.m_tStemsAndStemLoopsStartingBases;
	return *this;
}


/*
Inherited assignation operator
*/
CNAPrimaryStructure& CNASecondaryStructure::operator=(const val::biocpp::CNAPrimaryStructure &tSequence)
{
	TRACEI("CNAPrimaryStructure& CNASecondaryStructure::operator=(const val::biocpp::CNAPrimaryStructure &tSequence)");
	// check if it's a self-assignation
	if (this == &tSequence)
	{
		return *this;
	}
	// call inherited = operator
	CNAPrimaryStructure::operator=(tSequence);
	// copy fields
	m_bIsStructureParsed = false;
	m_iExternalBasesCount = 0;
	m_iStackedPairBasesCount = 0;
	m_iBulgeBasesCount = 0;
	m_iInternalLoopBasesCount = 0;
	m_iHairpinLoopBasesCount = 0;
	m_iMultibrachedLoopBasesCount = 0;
	m_iPseudoknotBasesCount = 0;
	m_iKissingHairpinBasesCount = 0;
	m_iKissingLoopBasesCount = 0;
	m_iTripleHelixBasesCount = 0;
	m_iStemSplitPosition = 0;
	int iPairsDataSize = GetLength() + 1;
	m_tNucleotidesPairs.resize(iPairsDataSize);
	m_tNucleotidesTertiaryInteractions.resize(iPairsDataSize);
	for (int i=0; i < iPairsDataSize; i++)
	{
		m_tNucleotidesPairs[i] = 0;
		m_tNucleotidesTertiaryInteractions[i] = 0;
	}
	m_tNucleotidesConfig.clear();
	m_tNucleotidesConfigTertiaryStructure.clear();
	m_iStemsCount = 0;
	m_iStemLoopsCount = 0;
	m_tStemsAndStemLoopsNature.clear();
	m_tStemsAndStemLoopsStartingBases.clear();
	return *this;
}


/*
Concatenates given sequence to the structure
*/
CNAPrimaryStructure& CNASecondaryStructure::operator+=(const val::biocpp::CNAPrimaryStructure &tRightSequence)
{
	TRACEI("CNAPrimaryStructure& CNASecondaryStructure::operator+=(const val::biocpp::CNAPrimaryStructure &tRightSequence)");

	// call inherited = operator
	CNAPrimaryStructure::operator+=(tRightSequence);
	// copy fields
	m_bIsStructureParsed = false;
	m_iExternalBasesCount = 0;
	m_iStackedPairBasesCount = 0;
	m_iBulgeBasesCount = 0;
	m_iInternalLoopBasesCount = 0;
	m_iHairpinLoopBasesCount = 0;
	m_iMultibrachedLoopBasesCount = 0;
	m_iPseudoknotBasesCount = 0;
	m_iKissingHairpinBasesCount = 0;
	m_iKissingLoopBasesCount = 0;
	m_iTripleHelixBasesCount = 0;
	m_iStemSplitPosition = 0;

	int iPairsDataSize = GetLength() + 1;
	m_tNucleotidesPairs.resize(iPairsDataSize, 0);
	m_tNucleotidesTertiaryInteractions.resize(iPairsDataSize, 0);

	m_tNucleotidesConfig.clear();
	m_tNucleotidesConfigTertiaryStructure.clear();
	m_iStemsCount = 0;
	m_iStemLoopsCount = 0;
	m_tStemsAndStemLoopsNature.clear();
	m_tStemsAndStemLoopsStartingBases.clear();
	return *this;
}


/*
Adds a nucleotide on the 3' end of the sequence.
*/
void CNASecondaryStructure::AddOn3Prime(const nucleic_acids::SNucleotideCode &tNucleotide)
{
	TRACEI("void CNASecondaryStructure::AddOn3Prime(const nucleic_acids::SNucleotideCode &tNucleotide)");
	// call parent method to update sequence
	CNAPrimaryStructure::AddOn3Prime(tNucleotide);

	// update structure info
	// new base has no interaction set
	m_tNucleotidesPairs.push_back(0);
	m_tNucleotidesTertiaryInteractions.push_back(0);
	// add base as external
	++m_iExternalBasesCount;
	m_tNucleotidesConfig.push_back(val::biocpp::na_secondary_structure::I_EXTERNAL_BASE);
	m_tNucleotidesConfigTertiaryStructure.push_back(val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG);	
}


/*
Adds a nucleotide on the 5' end of the sequence.
*/
void CNASecondaryStructure::AddOn5Prime(const nucleic_acids::SNucleotideCode &tNucleotide)
{
	TRACEI("void CNASecondaryStructure::AddOn5Prime(const nucleic_acids::SNucleotideCode &tNucleotide)");
	// call parent method to update sequence
	CNAPrimaryStructure::AddOn5Prime(tNucleotide);

	// update structure info
	// new base has no interaction set
	// note: 'front' element is not a part of the structure
	//       but inserting a 0 at front has the same behaviour
	//       as inserting it right after 'front'.
	m_tNucleotidesPairs.push_front(0);
	m_tNucleotidesTertiaryInteractions.push_front(0);

	// add as external base
	++m_iExternalBasesCount;
	m_tNucleotidesConfig.push_front(val::biocpp::na_secondary_structure::I_EXTERNAL_BASE);
	m_tNucleotidesConfigTertiaryStructure.push_front(val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG);

	// update other bases index and references
	//  update pairs index
	for (unsigned int i=2; i < m_tNucleotidesPairs.size(); i++)
	{
		// check if paired
		if (0 != m_tNucleotidesPairs[i])
		{
			m_tNucleotidesPairs[i] = m_tNucleotidesPairs[i] + 1;
		}
		// check if interaction
		if (0 != m_tNucleotidesTertiaryInteractions[i])
		{
			m_tNucleotidesTertiaryInteractions[i] = m_tNucleotidesTertiaryInteractions[i] + 1;
		}
	}
	// update stem/stem-loop index
	for (int i=0; i < (signed)m_tStemsAndStemLoopsStartingBases.size(); i++)
	{
		m_tStemsAndStemLoopsStartingBases[i].first = m_tStemsAndStemLoopsStartingBases[i].first + 1;
		m_tStemsAndStemLoopsStartingBases[i].second = m_tStemsAndStemLoopsStartingBases[i].second + 1;
	}
}


/*
Appends a nucleotide on the 3' end of the sequence.
*/
void CNASecondaryStructure::Append(const nucleic_acids::SNucleotideCode &tNucleotide)
{
	TRACEI("void CNASecondaryStructure::Append(const nucleic_acids::SNucleotideCode &tNucleotide)");
	// call parent method to update sequence
	CNAPrimaryStructure::Append(tNucleotide);

	// update structure info
	// new base has no interaction set
	m_tNucleotidesPairs.push_back(0);
	m_tNucleotidesTertiaryInteractions.push_back(0);
	// add base as external
	++m_iExternalBasesCount;
	m_tNucleotidesConfig.push_back(val::biocpp::na_secondary_structure::I_EXTERNAL_BASE);
	m_tNucleotidesConfigTertiaryStructure.push_back(val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG);	
}


/*
Appends a sequence to this sequence on the 3' end.
*/
int CNASecondaryStructure::AppendSequence(const std::string &strSequence, int iFromOffset, int iAppendLength)
{
	TRACEI("int CNASecondaryStructure::AppendSequence(const std::string &strSequence, int iFromOffset, int iAppendLength)");
	// call parent method to update sequence
	int iRet = CNAPrimaryStructure::AppendSequence(strSequence, iFromOffset, iAppendLength);

	// update structure info
	// new bases have no interaction set
	m_tNucleotidesPairs.insert(m_tNucleotidesPairs.end(), iRet, 0);
	m_tNucleotidesTertiaryInteractions.insert(m_tNucleotidesTertiaryInteractions.end(), iRet, 0);
	// add bases as externals
	m_iExternalBasesCount += iRet;
	m_tNucleotidesConfig.insert(m_tNucleotidesConfig.end(), iRet, val::biocpp::na_secondary_structure::I_EXTERNAL_BASE);
	m_tNucleotidesConfigTertiaryStructure.insert(m_tNucleotidesConfigTertiaryStructure.end(), iRet, val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG);

	return iRet;
}


/*
Appends a sequence to this sequence on the 3' end.
*/
int CNASecondaryStructure::AppendSequence(const CNAPrimaryStructure &tSequence, int iSequenceOffset, int iAppendLength)
{
	TRACEI("int CNASecondaryStructure::AppendSequence(const CNAPrimaryStructure &tSequence, int iSequenceOffset, int iAppendLength)");
	// call parent method to update sequence
	int iRet = CNAPrimaryStructure::AppendSequence(tSequence, iSequenceOffset, iAppendLength);

	// update structure info
	// new bases have no interactions set
	m_tNucleotidesPairs.insert(m_tNucleotidesPairs.end(), iRet, 0);
	m_tNucleotidesTertiaryInteractions.insert(m_tNucleotidesTertiaryInteractions.end(), iRet, 0);
	// add bases as externals
	m_iExternalBasesCount += iRet;
	m_tNucleotidesConfig.insert(m_tNucleotidesConfig.end(), iRet, val::biocpp::na_secondary_structure::I_EXTERNAL_BASE);
	m_tNucleotidesConfigTertiaryStructure.insert(m_tNucleotidesConfigTertiaryStructure.end(), iRet, val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG);

	return iRet;
}


/*
Inserts iCopyCount copies of tBaseCode at the given position.
*/
void CNASecondaryStructure::InsertAt(const nucleic_acids::SNucleotideCode &tNucleotide, int iAtPosition, int iCopyCount)
{
	TRACEI("void CNASecondaryStructure::InsertAt(const nucleic_acids::SNucleotideCode &tNucleotide, int iAtPosition, int iCopyCount)");
	CNAPrimaryStructure::InsertAt(tNucleotide, iAtPosition, iCopyCount);
	m_bIsStructureParsed = false;

	m_tNucleotidesPairs.insert(GetPairsDataIteratorAt(iAtPosition), iCopyCount, 0);
	m_tNucleotidesTertiaryInteractions.insert(GetTertiaryInteractionsDataIteratorAt(iAtPosition), iCopyCount, 0);
	for (unsigned int i=1; i < m_tNucleotidesPairs.size(); i++)
	{
		// check if paired
		if ((0 != m_tNucleotidesPairs[i]) && (iAtPosition <= m_tNucleotidesPairs[i]))
		{
			m_tNucleotidesPairs[i] = m_tNucleotidesPairs[i] + iCopyCount;
		}
		// check if interaction
		if ((0 != m_tNucleotidesTertiaryInteractions[i]) && (iAtPosition <= m_tNucleotidesTertiaryInteractions[i]))
		{
			m_tNucleotidesTertiaryInteractions[i] = m_tNucleotidesTertiaryInteractions[i] + iCopyCount;
		}
	}
}


/*
Inserts a sequence of nucleotides at the given position.
*/
int CNASecondaryStructure::InsertSequenceAt(const std::string &strSequence, int iAtPosition, int iFromOffset, int iInsertLength, int iCopyCount)
{
	TRACEI("int CNASecondaryStructure::InsertSequenceAt(const std::string &strSequence, int iAtPosition, int iFromOffset, int iInsertLength, int iCopyCount)");
	int iRet = CNAPrimaryStructure::InsertSequenceAt(strSequence, iAtPosition, iFromOffset, iInsertLength, iCopyCount);
	m_bIsStructureParsed = false;

	m_tNucleotidesPairs.insert(GetPairsDataIteratorAt(iAtPosition), iRet, 0);
	m_tNucleotidesTertiaryInteractions.insert(GetTertiaryInteractionsDataIteratorAt(iAtPosition), iRet, 0);
	for (unsigned int i=1; i < m_tNucleotidesPairs.size(); i++)
	{
		// check if paired
		if ((0 != m_tNucleotidesPairs[i]) && (iAtPosition <= m_tNucleotidesPairs[i]))
		{
			m_tNucleotidesPairs[i] = m_tNucleotidesPairs[i] + iRet;
		}
		// check if interaction
		if ((0 != m_tNucleotidesTertiaryInteractions[i]) && (iAtPosition <= m_tNucleotidesTertiaryInteractions[i]))
		{
			m_tNucleotidesTertiaryInteractions[i] = m_tNucleotidesTertiaryInteractions[i] + iRet;
		}
	}
	return iRet;
}


/*
Inserts a sequence of nucleotides at the given position.
*/
void CNASecondaryStructure::InsertSequenceAt(const CNAPrimaryStructure &tSequence, int iAtPosition, int iFromOffset, int iInsertLength, int iCopyCount)
{
	TRACEI("void CNASecondaryStructure::InsertSequenceAt(const CNAPrimaryStructure &tSequence, int iAtPosition, int iFromOffset, int iInsertLength, int iCopyCount)");
	int iOldLength = GetLength();
	CNAPrimaryStructure::InsertSequenceAt(tSequence, iAtPosition, iFromOffset, iInsertLength, iCopyCount);
	m_bIsStructureParsed = false;

	int iRealInsertLength = GetLength() - iOldLength;

	m_tNucleotidesPairs.insert(GetPairsDataIteratorAt(iAtPosition), iRealInsertLength, 0);
	m_tNucleotidesTertiaryInteractions.insert(GetTertiaryInteractionsDataIteratorAt(iAtPosition), iRealInsertLength, 0);
	for (unsigned int i=1; i < m_tNucleotidesPairs.size(); i++)
	{
		// check if paired
		if ((0 != m_tNucleotidesPairs[i]) && (iAtPosition <= m_tNucleotidesPairs[i]))
		{
			m_tNucleotidesPairs[i] = m_tNucleotidesPairs[i] + iRealInsertLength;
		}
		// check if interaction
		if ((0 != m_tNucleotidesTertiaryInteractions[i]) && (iAtPosition <= m_tNucleotidesTertiaryInteractions[i]))
		{
			m_tNucleotidesTertiaryInteractions[i] = m_tNucleotidesTertiaryInteractions[i] + iRealInsertLength;
		}
	}
}


/*
Removes iNucleotideCount nucleotides on the 3' end.
*/
void CNASecondaryStructure::RemoveOn3Prime(int iNucleotideCount)
{
	TRACEI("void CNASecondaryStructure::RemoveOn3Prime(int iNucleotideCount)");
	CNAPrimaryStructure::RemoveOn3Prime(iNucleotideCount);
	m_bIsStructureParsed = false;
	for (int i = 0; iNucleotideCount > i; i++)
	{
		// remove tertiary interactions if needed
		if (m_tNucleotidesTertiaryInteractions.back() != 0)
		{
			if (m_tNucleotidesPairs.back() != 0)
			{
				// base involved in both secondary and tertirary interactions
				m_tNucleotidesTertiaryInteractions[m_tNucleotidesPairs.back()] = 0;
			}
			else
			{
				// base is only involved in tertiary interactions
				// remove tertiary interactions
				m_tNucleotidesTertiaryInteractions[m_tNucleotidesPairs[m_tNucleotidesTertiaryInteractions.back()]] = 0;
			}
			m_tNucleotidesTertiaryInteractions[m_tNucleotidesTertiaryInteractions.back()] = 0;
		}
		m_tNucleotidesTertiaryInteractions.pop_back();	
		// remove secondary interaction
		if (m_tNucleotidesPairs.back() != 0)
		{
			m_tNucleotidesPairs[m_tNucleotidesPairs.back()] = 0;
		}
		m_tNucleotidesPairs.pop_back();
	}
}


/*
Removes iNucleotideCount nucleotides on the 5' end.
*/
void CNASecondaryStructure::RemoveOn5Prime(int iNucleotideCount)
{
	TRACEI("void CNASecondaryStructure::RemoveOn5Prime(int iNucleotideCount)");
	CNAPrimaryStructure::RemoveOn5Prime(iNucleotideCount);
	m_bIsStructureParsed = false;
	for (int i = 1; iNucleotideCount >= i; i++)
	{
		// remove tertiary interactions if needed
		if (m_tNucleotidesTertiaryInteractions[i] != 0)
		{
			if (m_tNucleotidesPairs[i] != 0)
			{
				// base involved in both secondary and tertirary interactions
				m_tNucleotidesTertiaryInteractions[m_tNucleotidesPairs[i]] = 0;
			}
			else
			{
				// base is only involved in tertiary interactions
				// remove tertiary interactions
				m_tNucleotidesTertiaryInteractions[m_tNucleotidesPairs[m_tNucleotidesTertiaryInteractions[i]]] = 0;
			}
			m_tNucleotidesTertiaryInteractions[m_tNucleotidesTertiaryInteractions[i]] = 0;
		}
		// remove pairs
		if (m_tNucleotidesPairs[i] != 0)
		{
			m_tNucleotidesPairs[m_tNucleotidesPairs[i]] = 0;
		}
	}

	// update index
	for (unsigned int i=iNucleotideCount+1; i < m_tNucleotidesPairs.size(); i++)
	{
		// check if paired
		if (0 != m_tNucleotidesPairs[i])
		{
			m_tNucleotidesPairs[i] = m_tNucleotidesPairs[i] - iNucleotideCount;
		}
		// check if interaction
		if (0 != m_tNucleotidesTertiaryInteractions[i])
		{
			m_tNucleotidesTertiaryInteractions[i] = m_tNucleotidesTertiaryInteractions[i] - iNucleotideCount;
		}
	}

	m_tNucleotidesTertiaryInteractions.erase(GetTertiaryInteractionsDataIteratorAt(1), GetTertiaryInteractionsDataIteratorAt(iNucleotideCount) + 1);
	m_tNucleotidesPairs.erase(GetPairsDataIteratorAt(1), GetPairsDataIteratorAt(iNucleotideCount) + 1);
}


/*
Removes iNucleotideCount nucleotides at the given position.
*/
void CNASecondaryStructure::RemoveAt(int iAtPosition, int iNucleotideCount)
{
	TRACEI("void CNASecondaryStructure::RemoveAt(int iAtPosition, int iNucleotideCount)");
	CNAPrimaryStructure::RemoveAt(iAtPosition, iNucleotideCount);
	m_bIsStructureParsed = false;

	for (int i=1; (signed)m_tNucleotidesPairs.size() > i; i++)
	{
		if ((i >= iAtPosition) && (i < iAtPosition + iNucleotideCount))
		{
			// remove pair of deleted bases
			if (0 != m_tNucleotidesPairs[i])
			{
				m_tNucleotidesPairs[m_tNucleotidesPairs[i]] = 0;
			}
			if (0 != m_tNucleotidesTertiaryInteractions[i])
			{
				m_tNucleotidesTertiaryInteractions[m_tNucleotidesTertiaryInteractions[i]] = 0;
			}
		}
		else
		{
			// check if paired and adjust pair
			if ((0 != m_tNucleotidesPairs[i]) && (iAtPosition <= m_tNucleotidesPairs[i]))
			{
				m_tNucleotidesPairs[i] = m_tNucleotidesPairs[i] - iNucleotideCount;
			}
			// check if interaction
			if ((0 != m_tNucleotidesTertiaryInteractions[i]) && (iAtPosition <= m_tNucleotidesTertiaryInteractions[i]))
			{
				m_tNucleotidesTertiaryInteractions[i] = m_tNucleotidesTertiaryInteractions[i] - iNucleotideCount;
			}
		}
	}

	CDequeIntIterator tIterBegin = GetPairsDataIteratorAt(iAtPosition);
	CDequeIntIterator tIterEnd = tIterBegin + iNucleotideCount;
	m_tNucleotidesPairs.erase(tIterBegin, tIterEnd);
	tIterBegin = GetTertiaryInteractionsDataIteratorAt(iAtPosition);
	tIterEnd = tIterBegin + iNucleotideCount;
	m_tNucleotidesTertiaryInteractions.erase(tIterBegin, tIterEnd);
}


/*
Removes iNucleotideCount nucleotides before the given nucleotide.
*/
void CNASecondaryStructure::RemoveBefore(int iAtPosition, int iNucleotideCount)
{
	TRACEI("void CNASecondaryStructure::RemoveBefore(int iAtPosition, int iNucleotideCount)");
	CNAPrimaryStructure::RemoveBefore(iAtPosition, iNucleotideCount);
	m_bIsStructureParsed = false;

	int iStartPosition = iAtPosition - iNucleotideCount;
	for (int i=0; i < (signed)m_tNucleotidesPairs.size(); i++)
	{
		if ((i >= iStartPosition) && (i < iAtPosition))
		{
			// remove pair
			if (0 != m_tNucleotidesPairs[i])
			{
				m_tNucleotidesPairs[m_tNucleotidesPairs[i]] = 0;
			}
			if (0 != m_tNucleotidesTertiaryInteractions[i])
			{
				m_tNucleotidesTertiaryInteractions[m_tNucleotidesTertiaryInteractions[i]] = 0;
			}
		}
		else
		{
			// check if paired and adjust pair
			if ((0 != m_tNucleotidesPairs[i]) && (iStartPosition <= m_tNucleotidesPairs[i]))
			{
				m_tNucleotidesPairs[i] = m_tNucleotidesPairs[i] - iNucleotideCount;
			}
			// check if interaction
			if ((0 != m_tNucleotidesTertiaryInteractions[i]) && (iStartPosition <= m_tNucleotidesTertiaryInteractions[i]))
			{
				m_tNucleotidesTertiaryInteractions[i] = m_tNucleotidesTertiaryInteractions[i] - iNucleotideCount;
			}
		}
	}

	CDequeIntIterator tIterBegin = GetPairsDataIteratorAt(iStartPosition);
	CDequeIntIterator tIterEnd = tIterBegin + iNucleotideCount;
	m_tNucleotidesPairs.erase(tIterBegin, tIterEnd);
	tIterBegin = GetTertiaryInteractionsDataIteratorAt(iStartPosition);
	tIterEnd = tIterBegin + iNucleotideCount;
	m_tNucleotidesTertiaryInteractions.erase(tIterBegin, tIterEnd);
}


/*
Removes iNucleotideCount nucleotides after the given nucleotide.
*/
void CNASecondaryStructure::RemoveAfter(int iAtPosition, int iNucleotideCount)
{
	TRACEI("void CNASecondaryStructure::RemoveAfter(int iAtPosition, int iNucleotideCount)");
	CNAPrimaryStructure::RemoveAfter(iAtPosition, iNucleotideCount);
	m_bIsStructureParsed = false;

	int iStartPosition = iAtPosition + 1;
	for (int i=0; i < (signed)m_tNucleotidesPairs.size(); i++)
	{
		if ((i >= iStartPosition) && (i < iStartPosition + iNucleotideCount))
		{
			// remove pair
			if (0 != m_tNucleotidesPairs[i])
			{
				m_tNucleotidesPairs[m_tNucleotidesPairs[i]] = 0;
			}
			if (0 != m_tNucleotidesTertiaryInteractions[i])
			{
				m_tNucleotidesTertiaryInteractions[m_tNucleotidesTertiaryInteractions[i]] = 0;
			}
		}
		else
		{
			// check if paired and adjust pair
			if ((0 != m_tNucleotidesPairs[i]) && (iStartPosition <= m_tNucleotidesPairs[i]))
			{
				m_tNucleotidesPairs[i] = m_tNucleotidesPairs[i] - iNucleotideCount;
			}
			// check if interaction
			if ((0 != m_tNucleotidesTertiaryInteractions[i]) && (iStartPosition <= m_tNucleotidesTertiaryInteractions[i]))
			{
				m_tNucleotidesTertiaryInteractions[i] = m_tNucleotidesTertiaryInteractions[i] - iNucleotideCount;
			}
		}
	}

	CDequeIntIterator tIterBegin = GetPairsDataIteratorAt(iStartPosition);
	CDequeIntIterator tIterEnd = tIterBegin + iNucleotideCount;
	m_tNucleotidesPairs.erase(tIterBegin, tIterEnd);
	tIterBegin = GetTertiaryInteractionsDataIteratorAt(iStartPosition);
	tIterEnd = tIterBegin + iNucleotideCount;
	m_tNucleotidesTertiaryInteractions.erase(tIterBegin, tIterEnd);
}


/**
@fn void CNASecondaryStructure::Pair(int iFirstPosition, int iSecondPosition)
*******************************************************************************************
Pairs together the two nucleotide specified by their positions in the sequence.
Position are 1-based, which means the first position in the sequence is 1.
If any of the given positions were already paired, the previous pair will be removed.

@param int iFirstPosition:\n
 position of the first nucleotide to pair.

@param int iSecondPosition:\n
 position of the other nucleotide to complete the pair.

@throw std::out_of_range:\n
 if the given position is outside the sequence.

*/
void CNASecondaryStructure::Pair(int iFirstPosition, int iSecondPosition)
{
	TRACEI("void CNASecondaryStructure::Pair(int iFirstPosition, int iSecondPosition)");
	
	int iLength = (signed)m_tNucleotidesPairs.size();
	// check given index
	if ((iFirstPosition <= 0) || (iSecondPosition <= 0) || (iFirstPosition > iLength) || (iSecondPosition > iLength))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::Pair(" << iFirstPosition << ", " << iSecondPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}

	m_bIsStructureParsed = false;
	// check if bases were already paired
	if (0 != m_tNucleotidesPairs[iFirstPosition])
	{
		// break other pair
		m_tNucleotidesPairs[m_tNucleotidesPairs[iFirstPosition]] = 0;
	}
	if (0 != m_tNucleotidesPairs[iSecondPosition])
	{
		// break other pair
		m_tNucleotidesPairs[m_tNucleotidesPairs[iSecondPosition]] = 0;
	}
	m_tNucleotidesPairs[iFirstPosition] = iSecondPosition;
	m_tNucleotidesPairs[iSecondPosition] = iFirstPosition;
}


/**
@fn void CNASecondaryStructure::PairInterval(int iFirstStartPosition, int iSecondStartPosition, int iPairsCount, bool bReverseDirectionForSecond)
*****************************************************************************************************************************************************************************
Pairs together the two intervals of nucleotides specified by the positions of
their first nucleotides. Position are 1-based, which means the first position
in the sequence is 1. If @c bReverseDirectionForSecond is set to true, the
second position can be considered as the last position of the nucleotide of
the second interval.

@param int iFirstStartPosition:\n
 position of the first nucleotide of the interval to pair.

@param int iSecondStartPosition:\n
 position of the first nucleotide of the other interval to complete pairing.

@param int iPairsCount:\n
 number of consecutive nucleotides to pair.

@param bool bIncreaseSecondPosition:\n
 if true, the first and the second positions are increased at the same time.\n
 By default, the first position is increased while the second position is
 decreased.
 Default value: false

@throw std::out_of_range:\n
 if the given positions lead outside the sequence.

*/
void CNASecondaryStructure::PairInterval(int iFirstStartPosition, int iSecondStartPosition, int iPairsCount, bool bIncreaseSecondPosition)
{
	TRACEI("void CNASecondaryStructure::PairInterval(int iFirstStartPosition, int iSecondStartPosition, int iPairsCount, bool bIncreaseSecondPosition)");
	int iLength = (signed)m_tNucleotidesPairs.size();
	// check given index
	if ((iFirstStartPosition <= 0)
	 || (iFirstStartPosition + iPairsCount > iLength)
	 || (iSecondStartPosition <= 0)
	 || ((true == bIncreaseSecondPosition) && (iSecondStartPosition + iPairsCount > iLength))
	 || ((false == bIncreaseSecondPosition) && ((iSecondStartPosition > iLength) || (iSecondStartPosition - iPairsCount <= 0))))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::PairInterval(" << iFirstStartPosition << ", " << iSecondStartPosition << ", " << iPairsCount << ", " << bIncreaseSecondPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}

	m_bIsStructureParsed = false;
	int iBase1 = iFirstStartPosition;
	int iBase2 = iSecondStartPosition ;
	if (true == bIncreaseSecondPosition)
	{
		for (int i=0; i < iPairsCount; i++)
		{
			// check if bases were already paired
			if (0 != m_tNucleotidesPairs[iBase1])
			{
				// break other pair
				m_tNucleotidesPairs[m_tNucleotidesPairs[iBase1]] = 0;
			}
			if (0 != m_tNucleotidesPairs[iBase2])
			{
				// break other pair
				m_tNucleotidesPairs[m_tNucleotidesPairs[iBase2]] = 0;
			}
			m_tNucleotidesPairs[iBase1] = iBase2;
			m_tNucleotidesPairs[iBase2] = iBase1;
			++iBase1;
			++iBase2;
		}
	}
	else
	{
		for (int i=0; i < iPairsCount; i++)
		{
			// check if bases were already paired
			if (0 != m_tNucleotidesPairs[iBase1])
			{
				// break other pair
				m_tNucleotidesPairs[m_tNucleotidesPairs[iBase1]] = 0;
			}
			if (0 != m_tNucleotidesPairs[iBase2])
			{
				// break other pair
				m_tNucleotidesPairs[m_tNucleotidesPairs[iBase2]] = 0;
			}
			m_tNucleotidesPairs[iBase1] = iBase2;
			m_tNucleotidesPairs[iBase2] = iBase1;
			++iBase1;
			--iBase2;
		}
	}
}


/**
@fn void CNASecondaryStructure::BreakPair(int iAtPosition)
*****************************************************************

*/
void CNASecondaryStructure::BreakPair(int iAtPosition)
{
	TRACEI("void CNASecondaryStructure::BreakPair(int iAtPosition)");
	int iLength = (signed)m_tNucleotidesPairs.size();
	// check given index
	if ((iAtPosition <= 0) || (iAtPosition > iLength))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::BreakPair(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}

	m_bIsStructureParsed = false;
	// check if base was already paired
	if (0 != m_tNucleotidesPairs[iAtPosition])
	{
		// break pair
		m_tNucleotidesPairs[m_tNucleotidesPairs[iAtPosition]] = 0;
		m_tNucleotidesPairs[iAtPosition] = 0;
	}
}


/**
@fn void CNASecondaryStructure::BreakPairInterval(int iStartPosition, int iPairsCount, bool bIncreasePosition)
********************************************************************************************************************

*/
void CNASecondaryStructure::BreakPairInterval(int iStartPosition, int iPairsCount, bool bIncreasePosition)
{
	TRACEI("void CNASecondaryStructure::BreakPairInterval(int iStartPosition, int iPairsCount, bool bIncreasePosition)");
	int iLength = (signed)m_tNucleotidesPairs.size();
	// check given index
	if ((iStartPosition <= 0)
	 || ((true == bIncreasePosition) && (iStartPosition + iPairsCount > iLength))
	 || ((false == bIncreasePosition) && (iStartPosition - iPairsCount <= 0)))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::BreakPairInterval(" << iStartPosition << ", " << iPairsCount << ", " << bIncreasePosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}

	m_bIsStructureParsed = false;
	int iBase1 = iStartPosition;
	if (true == bIncreasePosition)
	{
		for (int i=0; i < iPairsCount; i++)
		{
			// check if bases were already paired
			if (0 != m_tNucleotidesPairs[iBase1])
			{
				// break pair
				m_tNucleotidesPairs[m_tNucleotidesPairs[iBase1]] = 0;
			}
			m_tNucleotidesPairs[iBase1] = 0;
			++iBase1;
		}
	}
	else
	{
		for (int i=0; i < iPairsCount; i++)
		{
			// check if bases were already paired
			if (0 != m_tNucleotidesPairs[iBase1])
			{
				// break pair
				m_tNucleotidesPairs[m_tNucleotidesPairs[iBase1]] = 0;
			}
			m_tNucleotidesPairs[iBase1] = 0;
			iBase1--;
		}
	}
}


/**
@fn std::string CNASecondaryStructure::ToString(int iDisplayMode) const
***********************************************************************
Returns a string describing the secondary structure.
@li Mode I_ONE_LINE_NO_DOTBRACKETS returns the simple sequence.
@li Mode I_ONE_LINE_WITH_DOTBRACKETS returns a sequence mixed with brackets.
For example if there is only a G-C pair in the following structure, the
returned string will look like this: "AAAG(AAA)CAAA".
@li Mode I_TWO_LINES returns the sequence annoted with brackets on the second
line:
@code
AAAGAAACAAA
...(...)...
@endcode
@li Mode I_TEXT_STEMLOOPS returns a kind of graphical representation of stems
and stem-loops of the structure. The structure is split in stems and stem-loops
elements and each element includes the free bases on its 5' ends. The last
element also includes free bases on its last 3' end.
Elements are displayed one after the other and include 5' and 3' links
representation and also starting and ending bases of each segment.

@param int iDisplayMode:\n
 see @link #I_ONE_LINE_NO_BRACKETS String format constants @endlink for
 recognized values.\n
 Default value: I_ONE_LINE_NO_BRACKETS

@return std::string:\n
 a string of the sequence in the specified format.

@see @link #I_ONE_LINE_NO_BRACKETS String format constants @endlink
@see I_ONE_LINE_NO_BRACKETS
@see I_ONE_LINE_WITH_BRACKETS
@see I_TWO_LINES
@see I_TEXT_STEMLOOPS
@see I_TWO_LINES_STEMLOOPS
*/
std::string CNASecondaryStructure::ToString(int iDisplayMode) const
{
	TRACEI("std::string CNASecondaryStructure::ToString(int iDisplayMode) const");
	std::string strReturnString;

	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::ToString: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
/*
//+debug...
	for (unsigned int j=1; j<m_tNucleotidesPairs.size();j++)
	{
		std::cout << j << ": pair=" << m_tNucleotidesPairs[j] << ", ti=" << m_tNucleotidesTertiaryInteractions[j] << ", cfg=" << m_tNucleotidesConfig[j] << ", tcfg=" << m_tNucleotidesConfigTertiaryStructure[j] << std::endl;
	}
//...+debug
*/

	switch (iDisplayMode)
	{
	case na_secondary_structure::I_ONE_LINE_WITH_BRACKETS:
		{
			for (int i = 1; i < (signed)m_tNucleotidesPairs.size(); i++)
			{
				if (i == m_iStemSplitPosition)
				{
						strReturnString += val::biocpp::na_secondary_structure::C_STEM_SEPARATOR;
				}
				if (0 != m_tNucleotidesPairs[i])
				{
					// adds the brackets
					if (i < m_tNucleotidesPairs[i])
					{
						strReturnString += val::biocpp::SNucleotideCodeToChar(GetAt(i));
						strReturnString += '(';
					}
					else
					{
						strReturnString += ')';
						strReturnString += val::biocpp::SNucleotideCodeToChar(GetAt(i));
					}
				}
				else
				{
					// no H-bond, just add the base code
					strReturnString += val::biocpp::SNucleotideCodeToChar(GetAt(i));
				}
			}
			break;
		}
	case na_secondary_structure::I_TWO_LINES:
		{
			std::list<int> tFirstBaseList;
			std::string strFirstLine;
			std::string strSecondLine;
			unsigned char cCurrentChar = '.';
			bool bSamePseudoknot = false;
			for (int i = 1; i < (signed)m_tNucleotidesPairs.size(); i++)
			{
				if (i == m_iStemSplitPosition)
				{
					strFirstLine += val::biocpp::na_secondary_structure::C_STEM_SEPARATOR;
					strSecondLine += val::biocpp::na_secondary_structure::C_STEM_SEPARATOR;
				}
				strFirstLine += val::biocpp::SNucleotideCodeToChar(GetAt(i));
				if (0 != m_tNucleotidesPairs[i])
				{
					// is paired, check which bracket to use
					if (i < m_tNucleotidesPairs[i])
					{
						tFirstBaseList.push_back(i);
						strSecondLine += '(';
					}
					else if (i > m_tNucleotidesPairs[i])
					{
						if (tFirstBaseList.back() == m_tNucleotidesPairs[i])
						{
							// not a pseudoknot
							bSamePseudoknot = false;
							strSecondLine += ')';
							tFirstBaseList.pop_back();
						}
						else
						{
							// we got a pseudoknot
							if (true == bSamePseudoknot)
							{
								// still in the same pseudoknot
								strSecondLine += cCurrentChar;
							}
							else
							{
								// a new pseudoknot
								bSamePseudoknot = true;
								strSecondLine += cCurrentChar;
							}
							strSecondLine[m_tNucleotidesPairs[i] - 1] = cCurrentChar;
							// remove first base of the pair
							tFirstBaseList.remove(m_tNucleotidesPairs[i]);
						}
					}
				}
				else
				{
					// no pair, add a dot on the second line
					strSecondLine += '.';
				}
			}
			strReturnString += strFirstLine;
			strReturnString += '\n';
			strReturnString += strSecondLine;
			strReturnString += '\n';
			break;
		}
	case na_secondary_structure::I_TEXT_STEMLOOPS:
		{
			// check if we got stems in the structure
			if (0 == m_tStemsAndStemLoopsStartingBases.size())
			{
				// no stem
				strReturnString = "5'-" + CNAPrimaryStructure::ToString() + "-3'\n";
			}
			else
			{
				// stems or stem-loops
				for (int iStemIndex = 0; (signed)m_tStemsAndStemLoopsStartingBases.size() > iStemIndex; iStemIndex++)
				{
					std::string strLine1;
					std::string strLine2;
					std::string strLine3;
					std::string strLine4;
					std::string strLine5;
					int iStemStart = m_tStemsAndStemLoopsStartingBases[iStemIndex].first;
					int iStemEnd = m_tStemsAndStemLoopsStartingBases[iStemIndex].second;
					// get free bases on 5' end of the stem/stem-loop
					int iBaseIndex = iStemStart;
					--iBaseIndex;
					while ((0 < iBaseIndex) && (0 == m_tNucleotidesPairs[iBaseIndex]))
					{
						strLine1 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex--));
					}
					int iCutStart = iBaseIndex + 1;
					// get free bases on 3' end of the stem/stem-loop
					iBaseIndex = m_tNucleotidesPairs[iStemStart];
					++iBaseIndex;
					while (((signed)m_tNucleotidesPairs.size() > iBaseIndex) && (0 == m_tNucleotidesPairs[iBaseIndex]))
					{
						strLine5 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++));
					}
					int iCutEnd = iBaseIndex - 1;
					// check if we reached the 3' end
					if ((signed)m_tNucleotidesPairs.size() != iBaseIndex)
					{
						// not reached: drop bases (they'll be displayed on next stem/stem-loop)
						strLine5 = "";
						iCutEnd = m_tNucleotidesPairs[iStemStart];
					}
					// adjust string length
					while (strLine1.size() < strLine5.size())
					{
						strLine1 += "-";
					}
					while (strLine1.size() > strLine5.size())
					{
						strLine5 += "-";
					}
					// reverse the strings to have base in correct order
					std::reverse(strLine1.begin(), strLine1.end());
					std::reverse(strLine5.begin(), strLine5.end());
					// adjust other strings
					strLine2.insert(strLine2.begin(), strLine1.size(), ' ');
					strLine4 = strLine3 = strLine2;
					// start working on the stem/stem-loop
					iBaseIndex = iStemStart;
					int iBaseIndex2 = m_tNucleotidesPairs[iBaseIndex];
					if (val::biocpp::na_secondary_structure::I_STEMLOOP_TYPE == m_tStemsAndStemLoopsNature[iStemIndex])
					{
						// stem-loop
						// work on the stem part
						while ((iBaseIndex2 > iBaseIndex) && (((iBaseIndex2 - iBaseIndex) >= 3) || (val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP != m_tNucleotidesConfig[iBaseIndex])))
						{
							// check for pairs
							if (m_tNucleotidesPairs[iBaseIndex] == iBaseIndex2)
							{
								// pair type
								if (val::biocpp::SNucleotideCodeToIntCode(val::biocpp::GetComplementSNucleotideCode(GetAt(iBaseIndex))) == val::biocpp::SNucleotideCodeToIntCode(GetAt(iBaseIndex2)))
								{
									// canonical
									strLine1 += " ";
									strLine2 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++));
									strLine3 += "|";
									strLine4 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex2--));
									strLine5 += " ";
								}
								else if ((val::biocpp::nucleic_acids::Guanine == GetAt(iBaseIndex))
								          && (val::biocpp::nucleic_acids::Uridine == GetAt(iBaseIndex2))
								         || (val::biocpp::nucleic_acids::Uridine == GetAt(iBaseIndex))
								          && (val::biocpp::nucleic_acids::Guanine == GetAt(iBaseIndex2)))
								{
									// Wobble
									strLine1 += " ";
									strLine2 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++));
									strLine3 += "o";
									strLine4 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex2--));
									strLine5 += " ";
								}
								else
								{
									// non-canonical
									strLine1 += " ";
									strLine2 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++));
									strLine3 += ":";
									strLine4 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex2--));
									strLine5 += " ";
								}
							}
							else if ((0 == m_tNucleotidesPairs[iBaseIndex]) && (0 == m_tNucleotidesPairs[iBaseIndex2]))
							{
								// no pair
								strLine1 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++));
								strLine2 += " ";
								strLine3 += " ";
								strLine4 += " ";
								strLine5 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex2--));
							}
							else if (0 == m_tNucleotidesPairs[iBaseIndex2])
							{
								// more free bases on the 3' side
								strLine1 += "-";
								strLine2 += " ";
								strLine3 += " ";
								strLine4 += " ";
								strLine5 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex2--));
							}
							else if (0 == m_tNucleotidesPairs[iBaseIndex])
							{
								// more free bases on the 5' side
								strLine1 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++));
								strLine2 += " ";
								strLine3 += " ";
								strLine4 += " ";
								strLine5 += "-";
							}
							else
							{
								// ERROR: both 3' and 5' side bases are paired but with 2 other different basses!
								//+FIXME: throw an exception
							}
						}
						// work on the hairpin-loop terminal part
						if ((iBaseIndex2 - iBaseIndex) == 2)
						{
							if (val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP == m_tNucleotidesConfig[iBaseIndex - 1])
							{
								strLine1 += " ";
								strLine2 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++));
								strLine3 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++));
								strLine4 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex2--));
								strLine5 += " ";
							}
							else
							{
								strLine1 += "  ";
								strLine2 = strLine2 + val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++)) + " ";
								strLine3 = strLine3 + " " + val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++));
								strLine4 = strLine4 + val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex2--)) + " ";
								strLine5 += "  ";
							}
						}
						else if ((iBaseIndex2 - iBaseIndex) == 1)
						{
							if (val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP == m_tNucleotidesConfig[iBaseIndex - 1])
							{
								strLine1 += " ";
								strLine2 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++));
								strLine3 += " ";
								strLine4 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex2--));
								strLine5 += " ";
							}
							else
							{
								strLine1 += " ";
								strLine2 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++));
								strLine3 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex2--));
								strLine4 += " ";
								strLine5 += " ";
							}
						}
						else if (iBaseIndex2 == iBaseIndex)
						{
							strLine1 += " ";
							strLine2 += " ";
							strLine3 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++));
							strLine4 += " ";
							strLine5 += " ";
						}
						if (' ' == strLine1[0])
						{
							strReturnString += "   " + strLine1 + "\n";
							strReturnString += "5'-" + strLine2 + "\n";
							strReturnString += "   " + strLine3 + "\n";
							strReturnString += "3'-" + strLine4 + "\n";
							strReturnString += "   " + strLine5 + "\n";
						}
						else
						{
							strReturnString += "5'-" + strLine1 + "\n";
							strReturnString += "   " + strLine2 + "\n";
							strReturnString += "   " + strLine3 + "\n";
							strReturnString += "   " + strLine4 + "\n";
							strReturnString += "3'-" + strLine5 + "\n";
						}
						std::ostringstream stmStemData;
						stmStemData << "Stem-loop: " << iCutStart << "-" << iCutEnd << std::endl << std::endl;
						strReturnString += stmStemData.str();
					}
					else
					{
						// stem
						while ((iBaseIndex <= m_tNucleotidesPairs[iStemEnd]) || (iBaseIndex2 >= iStemEnd))
						{
							// check for pairs
							if (m_tNucleotidesPairs[iBaseIndex] == iBaseIndex2)
							{
								// pair type
								if (val::biocpp::SNucleotideCodeToIntCode(val::biocpp::GetComplementSNucleotideCode(GetAt(iBaseIndex))) == val::biocpp::SNucleotideCodeToIntCode(GetAt(iBaseIndex2)))
								{
									// canonical
									strLine1 += " ";
									strLine2 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++));
									strLine3 += "|";
									strLine4 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex2--));
									strLine5 += " ";
								}
								else if ((val::biocpp::nucleic_acids::Guanine == GetAt(iBaseIndex))
								          && (val::biocpp::nucleic_acids::Uridine == GetAt(iBaseIndex2))
								         || (val::biocpp::nucleic_acids::Uridine == GetAt(iBaseIndex))
								          && (val::biocpp::nucleic_acids::Guanine == GetAt(iBaseIndex2)))
								{
									// Wobble
									strLine1 += " ";
									strLine2 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++));
									strLine3 += "o";
									strLine4 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex2--));
									strLine5 += " ";
								}
								else
								{
									// non-canonical
									strLine1 += " ";
									strLine2 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++));
									strLine3 += ":";
									strLine4 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex2--));
									strLine5 += " ";
								}
							}
							else if ((0 == m_tNucleotidesPairs[iBaseIndex]) && (0 == m_tNucleotidesPairs[iBaseIndex2]))
							{
								// no pair
								strLine1 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++));
								strLine2 += " ";
								strLine3 += " ";
								strLine4 += " ";
								strLine5 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex2--));
							}
							else if (0 == m_tNucleotidesPairs[iBaseIndex2])
							{
								// more free bases on the 3' side
								strLine1 += "-";
								strLine2 += " ";
								strLine3 += " ";
								strLine4 += " ";
								strLine5 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex2--));
							}
							else if (0 == m_tNucleotidesPairs[iBaseIndex])
							{
								// more free bases on the 5' side
								strLine1 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++));
								strLine2 += " ";
								strLine3 += " ";
								strLine4 += " ";
								strLine5 += "-";
							}
							else
							{
								// ERROR: both 3' and 5' side bases are paired but with 2 other different basses!
								//+FIXME: throw an exception
							}
						}
						// add free bases
						if (0 == m_iStemSplitPosition)
						{
							while (0 == m_tNucleotidesPairs[iBaseIndex2])
							{
								strLine1 += "-";
								strLine2 += " ";
								strLine3 += " ";
								strLine4 += " ";
								strLine5 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex2--));
							}
						}
						else
						{
							while (m_iStemSplitPosition > iBaseIndex)
							{
								strLine1 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex++));
							}
							while (m_iStemSplitPosition <= iBaseIndex2)
							{
								strLine5 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex2--));
							}
							// adjust string length
							while (strLine1.size() < strLine5.size())
							{
								strLine1 += "-";
							}
							while (strLine1.size() > strLine5.size())
							{
								strLine5 += "-";
							}
							// adjust other strings
							while (strLine1.size() > strLine2.size())
							{
								strLine2 += " ";
								strLine3 += " ";
								strLine4 += " ";
							}
						}
						if (' ' == strLine1[0])
						{
							strLine1 = "   " + strLine1;
							strLine2 = "5'-" + strLine2;
							strLine3 = "   " + strLine3;
							strLine4 = "3'-" + strLine4;
							strLine5 = "   " + strLine5;
						}
						else
						{
							strLine1 = "5'-" + strLine1;
							strLine2 = "   " + strLine2;
							strLine3 = "   " + strLine3;
							strLine4 = "   " + strLine4;
							strLine5 = "3'-" + strLine5;
						}
						if (' ' == strLine1[strLine1.size() - 1])
						{
							strLine1 += "   ";
							strLine2 += "-3'";
							strLine3 += "   ";
							strLine4 += "-5'";
							strLine5 += "   ";
						}
						else
						{
							strLine1 += "-3'";
							strLine2 += "   ";
							strLine3 += "   ";
							strLine4 += "   ";
							strLine5 += "-5'";
						}
						strReturnString += strLine1 + "\n";
						strReturnString += strLine2 + "\n";
						strReturnString += strLine3 + "\n";
						strReturnString += strLine4 + "\n";
						strReturnString += strLine5 + "\n";
						std::ostringstream stmStemData;
						stmStemData << "Stem: " << iCutStart << "-" << iBaseIndex - 1 /*m_tNucleotidesPairs[iStemEnd]*/ << ", " << iBaseIndex2 + 1 << "-" << iCutEnd << std::endl << std::endl;
						strReturnString += stmStemData.str();
					}
				}
			}
			break;
		}
	case na_secondary_structure::I_TWO_LINES_STEMLOOPS:
		{
			// check if we got stems in the structure
			if (0 == m_tStemsAndStemLoopsStartingBases.size())
			{
				// no stem
				strReturnString = CNAPrimaryStructure::ToString();
				strReturnString = strReturnString + "\n" + std::string(strReturnString.length(), '.') + "\n";
			}
			else
			{
				// stems or stem-loops
				for (int iStemIndex = 0; (signed)m_tStemsAndStemLoopsStartingBases.size() > iStemIndex; iStemIndex++)
				{
					int iStemStart = m_tStemsAndStemLoopsStartingBases[iStemIndex].first;
					int iStemEnd = m_tStemsAndStemLoopsStartingBases[iStemIndex].second;
					std::string strLine1;
					std::string strLine2;
					int iPairBasesCount = 0;
					int iBulgeBasesCount = 0;
					int iInternalLoopBasesCount = 0;
					int iHairpinLoopBasesCount = 0;
					int iExternalBasesCount = 0;

					if (val::biocpp::na_secondary_structure::I_STEMLOOP_TYPE == m_tStemsAndStemLoopsNature[iStemIndex])
					{
						// stem-loop
						// start working on the stem/stem-loop
						//get previous free bases
						int iBaseIndex = iStemStart - 1;
						while ((0 < iBaseIndex) && (0 == m_tNucleotidesPairs[iBaseIndex]))
						{
							--iBaseIndex;
						}
						iStemStart = iBaseIndex + 1;
	
						//get following free bases
						iBaseIndex = iStemEnd + 1;
						while (((signed)m_tNucleotidesPairs.size() > iBaseIndex) && (0 == m_tNucleotidesPairs[iBaseIndex]))
						{
							++iBaseIndex;
						}
						iStemEnd = iBaseIndex - 1;
	
						for (iBaseIndex = iStemStart; iBaseIndex  <= iStemEnd; ++iBaseIndex)
						{
							strLine1 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex));
							// check for pairs
							if (0 < m_tNucleotidesPairs[iBaseIndex])
							{
								//we got a pair
								if (iBaseIndex < m_tNucleotidesPairs[iBaseIndex])
								{
									strLine2 += "(";
								}
								else
								{
									strLine2 += ")";
								}
							}
							else
							{
								// no pair
								strLine2 += ".";
							}
							// get config info
							switch (m_tNucleotidesConfig[iBaseIndex])
							{
								case val::biocpp::na_secondary_structure::I_STACKED_PAIR:
								{
									++iPairBasesCount;
									break;
								}
								case val::biocpp::na_secondary_structure::I_BULGE:
								{
									++iBulgeBasesCount;
									break;
								}
								case val::biocpp::na_secondary_structure::I_INTERNAL_LOOP:
								{
									++iInternalLoopBasesCount;
									break;
								}
								case val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP:
								{
									++iHairpinLoopBasesCount;
									break;
								}
								case val::biocpp::na_secondary_structure::I_EXTERNAL_BASE:
								case val::biocpp::na_secondary_structure::I_MULTIBRANCHED_LOOP:
								{
									++iExternalBasesCount;
									break;
								}
								case val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG:
								default:
								{
									break;
								}
							}
						}
						std::ostringstream stmStemData;
						stmStemData << strLine1 << std::endl << strLine2 << std::endl << "Stem-loop: " << m_tStemsAndStemLoopsStartingBases[iStemIndex].first << "-" << m_tStemsAndStemLoopsStartingBases[iStemIndex].second << std::endl << "Pair: " << iPairBasesCount << "; Bulge: " << iBulgeBasesCount << "; Internal: " << iInternalLoopBasesCount << "; Hairpin: " << iHairpinLoopBasesCount << "; External: " << iExternalBasesCount << std::endl << std::endl;
						strReturnString += stmStemData.str();
					}
					else
					{
						// stem

						// start working on the stem/stem-loop
						iStemEnd = m_tNucleotidesPairs[iStemEnd];
						int iStemStart2 = m_tStemsAndStemLoopsStartingBases[iStemIndex].second;
						int iStemEnd2 = m_tNucleotidesPairs[iStemStart];
						// we got: ..._iStemStart-iStemEnd_..._iStemStart2-iStemEnd2_...

						//get previous free bases (first segment)
						int iBaseIndex = iStemStart - 1;
						while ((0 < iBaseIndex) && (0 == m_tNucleotidesPairs[iBaseIndex]))
						{
							--iBaseIndex;
						}
						iStemStart = iBaseIndex + 1;
	
						//get following free bases (first segment)
						iBaseIndex = iStemEnd + 1;
						while (((signed)m_tNucleotidesPairs.size() > iBaseIndex) && (0 == m_tNucleotidesPairs[iBaseIndex]))
						{
							++iBaseIndex;
						}
						iStemEnd = iBaseIndex - 1;
	
						iBaseIndex = iStemStart2 - 1;
						//get previous free bases (second segment)
						while ((0 < iBaseIndex) && (0 == m_tNucleotidesPairs[iBaseIndex]))
						{
							--iBaseIndex;
						}
						iStemStart2 = iBaseIndex + 1;
	
						//get following free bases (second segment)
						iBaseIndex = iStemEnd2 + 1;
						while (((signed)m_tNucleotidesPairs.size() > iBaseIndex) && (0 == m_tNucleotidesPairs[iBaseIndex]))
						{
							++iBaseIndex;
						}
						iStemEnd2 = iBaseIndex - 1;
	
						for (iBaseIndex = iStemStart; iBaseIndex  <= iStemEnd; ++iBaseIndex)
						{
							strLine1 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex));
							// check for pairs
							if (0 < m_tNucleotidesPairs[iBaseIndex])
							{
								// we got a pair
								strLine2 += "(";
							}
							else
							{
								// no pair
								strLine2 += ".";
							}
							// get config info
							switch (m_tNucleotidesConfig[iBaseIndex])
							{
								case val::biocpp::na_secondary_structure::I_STACKED_PAIR:
								{
									++iPairBasesCount;
									break;
								}
								case val::biocpp::na_secondary_structure::I_BULGE:
								{
									++iBulgeBasesCount;
									break;
								}
								case val::biocpp::na_secondary_structure::I_INTERNAL_LOOP:
								{
									++iInternalLoopBasesCount;
									break;
								}
								case val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP:
								{
									++iHairpinLoopBasesCount;
									break;
								}
								case val::biocpp::na_secondary_structure::I_EXTERNAL_BASE:
								case val::biocpp::na_secondary_structure::I_MULTIBRANCHED_LOOP:
								{
									++iExternalBasesCount;
									break;
								}
								case val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG:
								default:
								{
									break;
								}
							}
						}
						strLine1 += val::biocpp::na_secondary_structure::C_STEM_SEPARATOR;
						strLine2 += val::biocpp::na_secondary_structure::C_STEM_SEPARATOR;
						for (iBaseIndex = iStemStart2; iBaseIndex  <= iStemEnd2; ++iBaseIndex)
						{
							strLine1 += val::biocpp::SNucleotideCodeToChar(GetAt(iBaseIndex));
							// check for pairs
							if (0 < m_tNucleotidesPairs[iBaseIndex])
							{
								//we got a pair
								strLine2 += ")";
							}
							else
							{
								// no pair
								strLine2 += ".";
							}
							// get config info
							switch (m_tNucleotidesConfig[iBaseIndex])
							{
								case val::biocpp::na_secondary_structure::I_STACKED_PAIR:
								{
									++iPairBasesCount;
									break;
								}
								case val::biocpp::na_secondary_structure::I_BULGE:
								{
									++iBulgeBasesCount;
									break;
								}
								case val::biocpp::na_secondary_structure::I_INTERNAL_LOOP:
								{
									++iInternalLoopBasesCount;
									break;
								}
								case val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP:
								{
									++iHairpinLoopBasesCount;
									break;
								}
								case val::biocpp::na_secondary_structure::I_EXTERNAL_BASE:
								case val::biocpp::na_secondary_structure::I_MULTIBRANCHED_LOOP:
								{
									++iExternalBasesCount;
									break;
								}
								case val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG:
								default:
								{
									break;
								}
							}
						}
						std::ostringstream stmStemData;
						stmStemData << strLine1 << std::endl << strLine2 << std::endl << "Stem: " << m_tStemsAndStemLoopsStartingBases[iStemIndex].first << "-" << m_tNucleotidesPairs[m_tStemsAndStemLoopsStartingBases[iStemIndex].second] << ", " << m_tStemsAndStemLoopsStartingBases[iStemIndex].second << "-" << m_tNucleotidesPairs[m_tStemsAndStemLoopsStartingBases[iStemIndex].first] << std::endl << "Pair: " << iPairBasesCount << "; Bulge: " << iBulgeBasesCount << "; Internal: " << iInternalLoopBasesCount << "; Hairpin: " << iHairpinLoopBasesCount << "; External: " << iExternalBasesCount << std::endl << std::endl;
						strReturnString += stmStemData.str();
					}
				}
			}
			break;
		}
	default:
		{
			// no secondary structure info
			strReturnString = CNAPrimaryStructure::ToString();
			break;
		}
	}

	return strReturnString;
}


/**
Returns 0 if there are no pair
*/
int CNASecondaryStructure::GetPairedBaseIndex(int iAtPosition) const
{
	TRACEI("int CNASecondaryStructure::GetPairedBaseIndex(int iAtPosition) const");
	int iLength = (signed)m_tNucleotidesPairs.size();
	// check given index
	if ((iAtPosition <= 0) || (iAtPosition > iLength))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::GetPairedBaseIndex(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	return m_tNucleotidesPairs[iAtPosition];
}


/**
Initializes nodes information (call it before using nodes info methods)
*/
void CNASecondaryStructure::ParseStructure()
{
	TRACEI("void CNASecondaryStructure::ParseStructure()");
	// check if already parsed
	if (true == m_bIsStructureParsed)
	{return;}

	m_iExternalBasesCount = 0;
	m_iStackedPairBasesCount = 0;
	m_iBulgeBasesCount = 0;
	m_iInternalLoopBasesCount = 0;
	m_iHairpinLoopBasesCount = 0;
	m_iMultibrachedLoopBasesCount = 0;
	m_iPseudoknotBasesCount = 0;
	m_iKissingHairpinBasesCount = 0;
	m_iKissingLoopBasesCount = 0;
	m_iTripleHelixBasesCount = 0;
	m_iStemSplitPosition = 0;
	
//	m_tNucleotidesConfig.clear();
	m_tNucleotidesConfig.resize(m_tNucleotidesPairs.size());
//	m_tNucleotidesConfigTertiaryStructure.clear();
	m_tNucleotidesConfigTertiaryStructure.resize(m_tNucleotidesPairs.size());
	for (unsigned int i=0; i < m_tNucleotidesPairs.size(); i++)
	{
		m_tNucleotidesConfig[i] = m_tNucleotidesConfigTertiaryStructure[i] = val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG;
	}

	m_iStemsCount = 0;
	m_iStemLoopsCount = 0;
	m_tStemsAndStemLoopsNature.clear();
	m_tStemsAndStemLoopsStartingBases.clear();

	// index1 = base_5'
	int iIndex1 = 1;
	std::deque<int> tStemsToParseStack;
	// tant_que index1 <= base_3'
	while (iIndex1 < (signed)m_tNucleotidesPairs.size())
	{
		// tant_que index1 non appariée et que index1 <= base_3'
		while ((iIndex1 < (signed)m_tNucleotidesPairs.size()) && (0 == m_tNucleotidesPairs[iIndex1]))
		{
			// marquer index1 comme I_EXTERNAL_BASE
			m_tNucleotidesConfig[iIndex1] = val::biocpp::na_secondary_structure::I_EXTERNAL_BASE;
			++m_iExternalBasesCount;
			// index1++
			++iIndex1;
		// fin_tant_que
		}
		// tant_que appariement et que index1 <= base_3'
		while ((iIndex1 < (signed)m_tNucleotidesPairs.size()) && (0 != m_tNucleotidesPairs[iIndex1]))
		{
			// si index1_paire > index1
			if (m_tNucleotidesPairs[iIndex1] > iIndex1)
			{
				// empiler index1 dans pile_paires
				tStemsToParseStack.push_back(iIndex1);
				// index1 = index1_paire + 1
				iIndex1 = m_tNucleotidesPairs[iIndex1] + 1;
			}
			// sinon // index1_paire < index1
			else
			{
				m_tNucleotidesConfig[iIndex1] = val::biocpp::na_secondary_structure::I_EXTERNAL_BASE;
				m_tNucleotidesConfig[m_tNucleotidesPairs[iIndex1]] = val::biocpp::na_secondary_structure::I_EXTERNAL_BASE;
				m_iExternalBasesCount += 2;
				// marquer index1 et index1_paire comme I_PSEUDOKNOT
				m_tNucleotidesConfigTertiaryStructure[iIndex1] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
				m_tNucleotidesConfigTertiaryStructure[m_tNucleotidesPairs[iIndex1]] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
				m_iPseudoknotBasesCount += 2;
				// index1++
				++iIndex1;
			// fin_si
			}
		// fin_tant_que
		}
	// fin_tant_que
	}
	//
	int iIndex = 0;
	// tant_que pile_paires non vide
	while (false == tStemsToParseStack.empty())
	{
		// stem_start = index1 = dépiler pile_paires
		int iStemStart = iIndex1 = tStemsToParseStack.front();
		tStemsToParseStack.pop_front();
		// index2 = index1_paire
		int iIndex2 = m_tNucleotidesPairs[iIndex1];
		// tant_que index1 < index2
		while (iIndex1 < iIndex2)
		{
			// faire
			do
			{
				// marquer index1 et index2 comme I_STACKED_PAIR
				m_tNucleotidesConfig[iIndex1] = val::biocpp::na_secondary_structure::I_STACKED_PAIR;
				m_tNucleotidesConfig[iIndex2] = val::biocpp::na_secondary_structure::I_STACKED_PAIR;
				m_iStackedPairBasesCount += 2;
				// index1++
				++iIndex1;
				// index2--
				iIndex2--;
			// fin_faire tant_que index1 apparié et index2 = index_paire et que index1 < index2
			} while ((iIndex2 == m_tNucleotidesPairs[iIndex1]) && (iIndex1 < iIndex2));
			// index1_start = index1
			int iStartIndex1 = iIndex1;
			// index2_start = index2
			int iStartIndex2 = iIndex2;
			std::deque<int> t5PrimeStack, t3PrimeStack;
			//si index1 > index2 // cas ou une tige se termine par une paire de base
			if (iIndex1 > iIndex2)
			{
//				// marquer stem_start comme un début et index1 comme fin de stem I_STEM_TYPE
//				m_tStemsAndStemLoopsNature.push_back(I_STEM_TYPE);
//				m_tStemsAndStemLoopsStartingBases.push_back(std::pair<int, int>(iStemStart, iIndex1));
				// marquer stem_start comme un début et index1 comme fin de stem I_STEM_TYPE
				m_tStemsAndStemLoopsNature.push_back(val::biocpp::na_secondary_structure::I_STEMLOOP_TYPE);
				m_tStemsAndStemLoopsStartingBases.push_back(std::pair<int, int>(iStemStart, m_tNucleotidesPairs[iStemStart]));
				++m_iStemLoopsCount;
			}
			// sinon si index1 = index2 // cas d'une seule base entre 2 bases d'une même paire
			else if (iIndex1 == iIndex2)
			{
				// si index1 appariée
				// marquer index1 comme I_HAIRPIN_LOOP
				m_tNucleotidesConfig[iIndex1] = val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP;
				++m_iHairpinLoopBasesCount;
				if (0 != m_tNucleotidesPairs[iIndex1])
				{
					//+FIXME: retirer l'appariement et le mettre dans les interactions tertiaires
					// marquer index1 comme I_PSEUDOKNOT
					m_tNucleotidesConfigTertiaryStructure[iIndex1] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
					++m_iPseudoknotBasesCount;
				// fin_si
				}
				// marquer stem_start comme un début de stem-loop I_STEMLOOP_TYPE
				m_tStemsAndStemLoopsNature.push_back(val::biocpp::na_secondary_structure::I_STEMLOOP_TYPE);
				m_tStemsAndStemLoopsStartingBases.push_back(std::pair<int, int>(iStemStart, m_tNucleotidesPairs[iStemStart]));
				++m_iStemLoopsCount;
			}
			// sinon
			else
			{
				// tant_que index1 <= index2 et que (index1 non appariée ou que index1_paire < index1_start ou que index1_paire > index2_start)
				while ((iIndex1 <= iIndex2) && ((0 == m_tNucleotidesPairs[iIndex1]) || (m_tNucleotidesPairs[iIndex1] < iStartIndex1) || (m_tNucleotidesPairs[iIndex1] > iStartIndex2)))
				{
					// empiler index1 dans pile_5'
					t5PrimeStack.push_back(iIndex1);
					// index1++
					++iIndex1;
				// fin_tant_que
				}
				// tant_que index1 < index2 et que (index2 non appariée ou que index2_paire < index1_start ou que index2_paire > index2_start)
				if (iIndex1 < iIndex2)
				{
					while ((0 == m_tNucleotidesPairs[iIndex2]) || (m_tNucleotidesPairs[iIndex2] < iStartIndex1) || (m_tNucleotidesPairs[iIndex2] > iStartIndex2))
					{
						// empiler index2 dans pile_3'
						t3PrimeStack.push_back(iIndex2);
						// index2--
						--iIndex2;
					// fin_tant_que
					}
				}
				// si index1 > index2 // aucune paire rencontrée de index1 vers index2
				if (iIndex1 > iIndex2)
				{
					// pour chaque base de la pile_5'
					while (false == t5PrimeStack.empty())
					{
						// index = depiler pile_5'
						iIndex = t5PrimeStack.back();
						t5PrimeStack.pop_back();
						//  marquer index comme I_HAIRPIN_LOOP
						m_tNucleotidesConfig[iIndex] = val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP;
						++m_iHairpinLoopBasesCount;
						// si index appariée
						if (0 != m_tNucleotidesPairs[iIndex])
						{
							// marquer index comme I_PSEUDOKNOT
							m_tNucleotidesConfigTertiaryStructure[iIndex] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
							++m_iPseudoknotBasesCount;
						// fin_si
						}
					// fin_pour
					}
					// marquer stem_start comme un début de stem-loop
					m_tStemsAndStemLoopsNature.push_back(val::biocpp::na_secondary_structure::I_STEMLOOP_TYPE);
					m_tStemsAndStemLoopsStartingBases.push_back(std::pair<int, int>(iStemStart, m_tNucleotidesPairs[iStemStart]));
					++m_iStemLoopsCount;
				}
				// sinon
				else
				{
					// si index1_paire = index2 // une seule paire rencontrée entre index1 et index2
					if (m_tNucleotidesPairs[iIndex1] == iIndex2)
					{
						//  si pile_5' et pile_3' non vides
						if ((false == t5PrimeStack.empty()) && (false == t3PrimeStack.empty()))
						{
							// dépiler et marquer chaque base de pile_5' et de pile_3' comme I_INTERNAL_LOOP
							while (false == t5PrimeStack.empty())
							{
								// index = depiler pile_5'
								iIndex = t5PrimeStack.back();
								t5PrimeStack.pop_back();
								// marquer index comme I_BULGE
								m_tNucleotidesConfig[iIndex] = val::biocpp::na_secondary_structure::I_INTERNAL_LOOP;
								++m_iInternalLoopBasesCount;
								// si index appariée
								if (0 != m_tNucleotidesPairs[iIndex])
								{
									// marquer index comme I_PSEUDOKNOT
									m_tNucleotidesConfigTertiaryStructure[iIndex] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
									++m_iPseudoknotBasesCount;
								// fin_si
								}
							}
							while (false == t3PrimeStack.empty())
							{
								// index = depiler pile_3'
								iIndex = t3PrimeStack.back();
								t3PrimeStack.pop_back();
								// marquer index comme I_BULGE
								m_tNucleotidesConfig[iIndex] = val::biocpp::na_secondary_structure::I_INTERNAL_LOOP;
								++m_iInternalLoopBasesCount;
								// si index appariée
								if (0 != m_tNucleotidesPairs[iIndex])
								{
									// marquer index comme I_PSEUDOKNOT
									m_tNucleotidesConfigTertiaryStructure[iIndex] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
									++m_iPseudoknotBasesCount;
								// fin_si
								}
							}
						}
						// sinon
						else
						{
							// pour chaque base de la pile_5'
							while (false == t5PrimeStack.empty())
							{
								// index = depiler pile_5'
								iIndex = t5PrimeStack.back();
								t5PrimeStack.pop_back();
								// marquer index comme I_BULGE
								m_tNucleotidesConfig[iIndex] = val::biocpp::na_secondary_structure::I_BULGE;
								++m_iBulgeBasesCount;
								// si index appariée
								if (0 != m_tNucleotidesPairs[iIndex])
								{
									// marquer index comme I_PSEUDOKNOT
									m_tNucleotidesConfigTertiaryStructure[iIndex] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
									++m_iPseudoknotBasesCount;
								// fin_si
								}
							// fin_pour
							}
							// pour chaque base de la pile_3'
							while (false == t3PrimeStack.empty())
							{
								// index = depiler pile_3'
								iIndex = t3PrimeStack.back();
								t3PrimeStack.pop_back();
								// marquer index comme I_BULGE
								m_tNucleotidesConfig[iIndex] = val::biocpp::na_secondary_structure::I_BULGE;
								++m_iBulgeBasesCount;
								// si index appariée
								if (0 != m_tNucleotidesPairs[iIndex])
								{
									// marquer index comme I_PSEUDOKNOT
									m_tNucleotidesConfigTertiaryStructure[iIndex] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
									++m_iPseudoknotBasesCount;
								// fin_si
								}
							// fin_pour
							}
						// fin_si
						}
					}
					// sinon // au moins 2 paires distinctes dans la boucle
					else
					{
						// tant_que index1 < index2
						while (iIndex1 < iIndex2)
						{
							// si index1_paire < index1
							if (m_tNucleotidesPairs[iIndex1] < iIndex1)
							{
								// marquer index1 comme I_MULTIBRANCHED_LOOP
								m_tNucleotidesConfig[iIndex1] = val::biocpp::na_secondary_structure::I_MULTIBRANCHED_LOOP;
								++m_iMultibrachedLoopBasesCount;
								// marquer index1 comme I_PSEUDOKNOT
								m_tNucleotidesConfigTertiaryStructure[iIndex1] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
								++m_iPseudoknotBasesCount;
								// index1++
								++iIndex1;
							}
							// sinon
							else
							{
								// empiler index1 dans pile_paires
								tStemsToParseStack.push_back(iIndex1);
								// index1 = index1_paire + 1
								iIndex1 = m_tNucleotidesPairs[iIndex1] + 1;
							// fin_si
							}
							// tant_que index1 non appariée et que index1 < index2
							while ((0 == m_tNucleotidesPairs[iIndex1]) && (iIndex1 < iIndex2))
							{
								// marquer index1 comme I_MULTIBRANCHED_LOOP
								m_tNucleotidesConfig[iIndex1] = val::biocpp::na_secondary_structure::I_MULTIBRANCHED_LOOP;
								++m_iMultibrachedLoopBasesCount;
								// index1++
								++iIndex1;
							// fin_tant_que
							}
						// fin_tant_que
						}
						// pour chaque base de la pile_5'
						while (false == t5PrimeStack.empty())
						{
							// index = depiler pile_5'
							iIndex = t5PrimeStack.back();
							t5PrimeStack.pop_back();
							// marquer index comme I_MULTIBRANCHED_LOOP
							m_tNucleotidesConfig[iIndex] = val::biocpp::na_secondary_structure::I_MULTIBRANCHED_LOOP;
							++m_iMultibrachedLoopBasesCount;
							// si index appariée
							if (0 != m_tNucleotidesPairs[iIndex])
							{
								// marquer index comme I_PSEUDOKNOT
								m_tNucleotidesConfigTertiaryStructure[iIndex] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
								++m_iPseudoknotBasesCount;
							// fin_si
							}
						// fin_pour
						}
						// pour chaque base de la pile_3'
						while (false == t3PrimeStack.empty())
						{
							// index = depiler pile_3'
							iIndex = t3PrimeStack.back();
							t3PrimeStack.pop_back();
							// marquer index comme I_MULTIBRANCHED_LOOP
							m_tNucleotidesConfig[iIndex] = val::biocpp::na_secondary_structure::I_MULTIBRANCHED_LOOP;
							++m_iMultibrachedLoopBasesCount;
							// si index appariée
							if (0 != m_tNucleotidesPairs[iIndex])
							{
								// marquer index comme I_PSEUDOKNOT
								m_tNucleotidesConfigTertiaryStructure[iIndex] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
								++m_iPseudoknotBasesCount;
							// fin_si
							}
						// fin_pour
						}
						// marquer stem_start comme un début de stem et index2_start + 1 comme une fin de stem
						m_tStemsAndStemLoopsNature.push_back(val::biocpp::na_secondary_structure::I_STEM_TYPE);
						m_tStemsAndStemLoopsStartingBases.push_back(std::pair<int, int>(iStemStart, iStartIndex2 + 1));
						++m_iStemsCount;
					// fin_si
					}
				// fin_si
				}
			// fin_si
			}
		// fin_tant_que
		}
	// fin_tant_que
	}
	//
	//+++ pseudo-code ci-dessous à revoir totalement
	// pour chaque baseP1 marquée I_PSEUDOKNOT
	//  chercher la 1ère base de paire (non I_PSEUDOKNOT ou I_KISSING_HAIRPIN ou I_KISSING_LOOP) en 5' et celle en 3'
	//  si ces 2 bases forment la même paire
	//   marquer baseP1 comme I_KISSING_HAIRPIN
	//   prendre la baseP2 appariée avec baseP1
	//   chercher la 1ère base de paire (non I_PSEUDOKNOT ou I_KISSING_HAIRPIN ou I_KISSING_LOOP) en 5' et celle en 3'
	//   si ces 2 bases forment la même paire
	//    marquer baseP2 comme I_KISSING_HAIRPIN
	//   sinon
	//    marquer baseP2 comme I_KISSING_LOOP
	//   fin_si
	//  fin_si
	// fin_pour
/*
	//+debug... 
	for (unsigned int j=1; j<m_tNucleotidesPairs.size();j++)
	{
		std::cout << j << ": pair=" << m_tNucleotidesPairs[j] << ", ti=" << m_tNucleotidesTertiaryInteractions[j] << ", cfg=" << m_tNucleotidesConfig[j] << ", tcfg=" << m_tNucleotidesConfigTertiaryStructure[j] << std::endl;
	}
	//...+debug
*/
	m_bIsStructureParsed = true;
}




/**
Returns the configuration of a nucleotide.
*/
int CNASecondaryStructure::GetNucleotideConfig(int iAtPosition) const
{
	TRACEI("int CNASecondaryStructure::GetNucleotideConfig(int iAtPosition) const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::GetNucleotideConfig: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	int iLength = (signed)m_tNucleotidesConfig.size();
	if ((iAtPosition <= 0) || (iAtPosition > iLength))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::GetNucleotideConfig(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}

	return m_tNucleotidesConfig[iAtPosition];
}


/**
Returns the tertiary configuration of a nucleotide.
*/
int CNASecondaryStructure::GetNucleotideTertiaryConfig(int iAtPosition) const
{
	TRACEI("int CNASecondaryStructure::GetNucleotideTertiaryConfig(int iAtPosition) const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::GetNucleotideTertiaryConfig: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	int iLength = (signed)m_tNucleotidesConfig.size();
	if ((iAtPosition <= 0) || (iAtPosition > iLength))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::GetNucleotideTertiaryConfig(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}

	return m_tNucleotidesConfigTertiaryStructure[iAtPosition];
}


/**
Tells if a nucleotide is an extern.
*/
bool CNASecondaryStructure::IsExternalBase(int iAtPosition) const
{
	TRACEI("bool CNASecondaryStructure::IsExternalBase(int iAtPosition) const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::IsExternalBase: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	int iLength = (signed)m_tNucleotidesConfig.size();
	if ((iAtPosition <= 0) || (iAtPosition > iLength))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::IsExternalBase(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	bool bRet = false;
	if (val::biocpp::na_secondary_structure::I_EXTERNAL_BASE == m_tNucleotidesConfig[iAtPosition])
	{
		bRet = true;
	}
	return bRet;
}


/**
Tells if a nucleotide is paired.
*/
bool CNASecondaryStructure::IsPaired(int iAtPosition) const
{
	TRACEI("bool CNASecondaryStructure::IsPaired(int iAtPosition) const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::IsPaired: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	int iLength = (signed)m_tNucleotidesConfig.size();
	if ((iAtPosition <= 0) || (iAtPosition > iLength))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::IsPaired(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	bool bRet = false;
	if (val::biocpp::na_secondary_structure::I_STACKED_PAIR == m_tNucleotidesConfig[iAtPosition])
	{
		bRet = true;
	}
	return bRet;
}


/**
Tells if a nucleotide is in a bulge.
*/
bool CNASecondaryStructure::IsInBulge(int iAtPosition) const
{
	TRACEI("bool CNASecondaryStructure::IsInBulge(int iAtPosition) const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::IsInBulge: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	int iLength = (signed)m_tNucleotidesConfig.size();
	if ((iAtPosition <= 0) || (iAtPosition > iLength))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::IsInBulge(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	bool bRet = false;
	if (val::biocpp::na_secondary_structure::I_BULGE == m_tNucleotidesConfig[iAtPosition])
	{
		bRet = true;
	}
	return bRet;
}


/**
Tells if a nucleotide is in an internal loop.
*/
bool CNASecondaryStructure::IsInInternalLoop(int iAtPosition) const
{
	TRACEI("bool CNASecondaryStructure::IsInInternalLoop(int iAtPosition) const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::IsInInternalLoop: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	int iLength = (signed)m_tNucleotidesConfig.size();
	if ((iAtPosition <= 0) || (iAtPosition > iLength))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::IsInInternalLoop(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	bool bRet = false;
	if (val::biocpp::na_secondary_structure::I_INTERNAL_LOOP == m_tNucleotidesConfig[iAtPosition])
	{
		bRet = true;
	}
	return bRet;
}


/**
Tells if a nucleotide is in a hairpin loop.
*/
bool CNASecondaryStructure::IsInHairpinLoop(int iAtPosition) const
{
	TRACEI("bool CNASecondaryStructure::IsInHairpinLoop(int iAtPosition) const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::IsInHairpinLoop: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	int iLength = (signed)m_tNucleotidesConfig.size();
	if ((iAtPosition <= 0) || (iAtPosition > iLength))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::IsInHairpinLoop(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	bool bRet = false;
	if (val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP == m_tNucleotidesConfig[iAtPosition])
	{
		bRet = true;
	}
	return bRet;
}


/**
Tells if a nucleotide is in a multibranched loop.
*/
bool CNASecondaryStructure::IsInMultibranchedLoop(int iAtPosition) const
{
	TRACEI("bool CNASecondaryStructure::IsInMultibranchedLoop(int iAtPosition) const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::IsInMultibranchedLoop: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	int iLength = (signed)m_tNucleotidesConfig.size();
	if ((iAtPosition <= 0) || (iAtPosition > iLength))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::IsInMultibranchedLoop(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	bool bRet = false;
	if (val::biocpp::na_secondary_structure::I_MULTIBRANCHED_LOOP == m_tNucleotidesConfig[iAtPosition])
	{
		bRet = true;
	}
	return bRet;
}


/**
Tells if a nucleotide is in a pseudoknot loop.
*/
bool CNASecondaryStructure::IsInPseudoknot(int iAtPosition) const
{
	TRACEI("bool CNASecondaryStructure::IsInPseudoknot(int iAtPosition) const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::IsInPseudoknot: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	int iLength = (signed)m_tNucleotidesConfig.size();
	if ((iAtPosition <= 0) || (iAtPosition > iLength))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::IsInPseudoknot(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	bool bRet = false;
	if (val::biocpp::na_secondary_structure::I_PSEUDOKNOT == m_tNucleotidesConfig[iAtPosition])
	{
		bRet = true;
	}
	return bRet;
}


/**
Tells if a nucleotide is in a kissing hairpin loop.
*/
bool CNASecondaryStructure::IsInKissingHairpin(int iAtPosition) const
{
	TRACEI("bool CNASecondaryStructure::IsInKissingHairpin(int iAtPosition) const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::IsInKissingHairpin: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	int iLength = (signed)m_tNucleotidesConfig.size();
	if ((iAtPosition <= 0) || (iAtPosition > iLength))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::IsInKissingHairpin(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	bool bRet = false;
	if (val::biocpp::na_secondary_structure::I_KISSING_HAIRPIN == m_tNucleotidesConfig[iAtPosition])
	{
		bRet = true;
	}
	return bRet;
}


/**
Tells if a nucleotide is in a kissing loop.
*/
bool CNASecondaryStructure::IsInKissingLoop(int iAtPosition) const
{
	TRACEI("bool CNASecondaryStructure::IsInKissingLoop(int iAtPosition) const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::IsInKissingLoop: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	int iLength = (signed)m_tNucleotidesConfig.size();
	if ((iAtPosition <= 0) || (iAtPosition > iLength))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::IsInKissingLoop(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	bool bRet = false;
	if (val::biocpp::na_secondary_structure::I_KISSING_LOOP == m_tNucleotidesConfig[iAtPosition])
	{
		bRet = true;
	}
	return bRet;
}


/**
Tells if a nucleotide is in a triple helix.
*/
bool CNASecondaryStructure::IsInTripleHelix(int iAtPosition) const
{
	TRACEI("bool CNASecondaryStructure::IsInTripleHelix(int iAtPosition) const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::IsInTripleHelix: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	int iLength = (signed)m_tNucleotidesConfig.size();
	if ((iAtPosition <= 0) || (iAtPosition > iLength))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::IsInTripleHelix(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	bool bRet = false;
	if (val::biocpp::na_secondary_structure::I_TRIPLE_HELIX == m_tNucleotidesConfig[iAtPosition])
	{
		bRet = true;
	}
	return bRet;
}


/**
Returns the number of external bases.
*/
int CNASecondaryStructure::GetExternalBasesCount() const
{
	TRACEI("int CNASecondaryStructure::GetExternalBasesCount() const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::GetExternalBasesCount: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	return m_iExternalBasesCount;
}


/**
Returns the number of bases in stacked pairs (a pair is 2 nucleotides).
*/
int CNASecondaryStructure::GetStackedPairBasesCount() const
{
	TRACEI("int CNASecondaryStructure::GetStackedPairBasesCount() const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::GetStackedPairBasesCount: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	return m_iStackedPairBasesCount;
}


/**
Returns the number of bulges.
*/
int CNASecondaryStructure::GetBulgeBasesCount() const
{
	TRACEI("int CNASecondaryStructure::GetBulgeBasesCount() const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::GetBulgeBasesCount: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	return m_iBulgeBasesCount;
}


/**
Returns the number of bases in internal loops.
*/
int CNASecondaryStructure::GetInternalLoopBasesCount() const
{
	TRACEI("int CNASecondaryStructure::GetInternalLoopBasesCount() const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::GetInternalLoopBasesCount: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	return m_iInternalLoopBasesCount;
}


/**
Returns the number of bases in hairpin loops.
*/
int CNASecondaryStructure::GetHairpinLoopBasesCount() const
{
	TRACEI("int CNASecondaryStructure::GetHairpinLoopBasesCount() const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::GetHairpinLoopBasesCount: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	return m_iHairpinLoopBasesCount;
}


/**
Returns the number of bases in multibranched loops.
*/
int CNASecondaryStructure::GetMultibranchedLoopBasesCount() const
{
	TRACEI("int CNASecondaryStructure::GetMultibranchedLoopBasesCount() const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::GetMultibranchedLoopBasesCount: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	return m_iMultibrachedLoopBasesCount;
}


/**
Returns the number of base in pseudoknot.
*/
int CNASecondaryStructure::GetPseudoknotBasesCount() const
{
	TRACEI("int CNASecondaryStructure::GetPseudoknotBasesCount() const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::GetPseudoknotBasesCount: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	return m_iPseudoknotBasesCount;
}


/**
Returns the number of base in kissing hairpin loops.
*/
int CNASecondaryStructure::GetKissingHairpinBasesCount() const
{
	TRACEI("int CNASecondaryStructure::GetKissingHairpinBasesCount() const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::GetKissingHairpinBasesCount: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	return m_iKissingHairpinBasesCount;
}


/**
Returns the number of base in kissing loops.
*/
int CNASecondaryStructure::GetKissingLoopBasesCount() const
{
	TRACEI("int CNASecondaryStructure::GetKissingLoopBasesCount() const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::GetKissingLoopBasesCount: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	return m_iKissingLoopBasesCount;
}


/**
Returns the number of base in triple helix.
*/
int CNASecondaryStructure::GetTripleHelixBasesCount() const
{
	TRACEI("int CNASecondaryStructure::GetTripleHelixBasesCount() const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::GetTripleHelixBasesCount: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	return m_iTripleHelixBasesCount;
}


/**
Returns the number of stem-loops.
*/
int CNASecondaryStructure::GetStemLoopsCount() const
{
	TRACEI("int CNASecondaryStructure::GetStemLoopsCount() const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::GetStemLoopsCount" + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	return m_iStemLoopsCount;
}


/**
Returns the number of stem (without hairpin loops).
*/
int CNASecondaryStructure::GetStemsCount() const
{
	TRACEI("int CNASecondaryStructure::GetStemsCount() const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::GetStemsCount: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	return m_iStemsCount;
}


/**
Returns the number of stem and stem-loops.
*/
int CNASecondaryStructure::GetStemsAndStemLoopsCount() const
{
	TRACEI("int CNASecondaryStructure::GetStemsAndStemLoopsCount() const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::GetStemsAndStemLoopsCount: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	return (signed)m_tStemsAndStemLoopsStartingBases.size();
}

/*
Tells if given stem/stem-loop index correspond to a stem.
1-based
*/
bool CNASecondaryStructure::IsStem(int iStemLoopNumber) const
{
	TRACEI("bool CNASecondaryStructure::IsStem(int iStemLoopNumber) const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::IsStem: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	if ((0 >= iStemLoopNumber) || ((signed)m_tStemsAndStemLoopsStartingBases.size() < iStemLoopNumber))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::IsStem(" << iStemLoopNumber << "): " << STR_ERROR_INVALID_STEM_OR_STEMLOOP << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	bool bRet = false;
	if (val::biocpp::na_secondary_structure::I_STEM_TYPE == m_tStemsAndStemLoopsNature[iStemLoopNumber-1])
	{
		bRet = true;
	}
	return bRet;
}

/*
Tells if given stem/stem-loop index correspond to a stem-loop.
1-based
*/
bool CNASecondaryStructure::IsStemLoop(int iStemLoopNumber) const
{
	TRACEI("bool CNASecondaryStructure::IsStemLoop(int iStemLoopNumber) const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::IsStemLoop: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	if ((0 >= iStemLoopNumber) || ((signed)m_tStemsAndStemLoopsStartingBases.size() < iStemLoopNumber))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::IsStemLoop(" << iStemLoopNumber << "): " << STR_ERROR_INVALID_STEM_OR_STEMLOOP << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	bool bRet = false;
	if (val::biocpp::na_secondary_structure::I_STEMLOOP_TYPE == m_tStemsAndStemLoopsNature[iStemLoopNumber - 1])
	{
		bRet = true;
	}
	return bRet;
}


/**
Returns the secondary structure of a stem-loop.
*/
CNASecondaryStructure CNASecondaryStructure::GetStemOrStemLoop(int iStemLoopNumber, int iCutMethod) const
{
	TRACEI("CNASecondaryStructure CNASecondaryStructure::GetStemOrStemLoop(int iStemLoopNumber, int iCutMethod) const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::GetStemOrStemLoop: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	if ((0 >= iStemLoopNumber) || ((signed)m_tStemsAndStemLoopsStartingBases.size() < iStemLoopNumber))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::GetStemOrStemLoop(" << iStemLoopNumber << ", " << iCutMethod << "): " << STR_ERROR_INVALID_STEM_OR_STEMLOOP << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}

	CNASecondaryStructure tStemLoop;
	if (val::biocpp::na_secondary_structure::I_STEM_TYPE == m_tStemsAndStemLoopsNature[iStemLoopNumber-1])
	{
		int iStartPosition1 = m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].first;
		int iEndPosition1 = m_tNucleotidesPairs[m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].second];
		int iStartPosition2 = m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].second;
		int iEndPosition2 = m_tNucleotidesPairs[m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].first];
		// stems or stem-loops must start with a pair
		if ((0 == iEndPosition1) || (iEndPosition1 < iStartPosition1)
		    || (0 == iEndPosition2) || (iEndPosition2 < iStartPosition2)
		    || (iStartPosition2 <= iEndPosition1))
		{
			throw std::runtime_error("CNASecondaryStructure::GetStemOrStemLoop: " + STR_ERROR_CORRUPTED_STRUCTURE_DATA);
		}

		// see where to cut
		switch (iCutMethod)
		{
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_5PRIME_FREE_BASES:
			{
				do
				{
					--iStartPosition1;
				}
				while ((0 < iStartPosition1) && (0 == m_tNucleotidesPairs[iStartPosition1]));
				++iStartPosition1;
				do
				{
					--iStartPosition2;
				}
				while ((0 < iStartPosition2) && (0 == m_tNucleotidesPairs[iStartPosition2]));
				++iStartPosition2;
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_5PRIME_FREE_BASES_WITH_LAST_3PRIME_FREE_BASES:
			{
				do
				{
					--iStartPosition1;
				}
				while ((0 < iStartPosition1) && (0 == m_tNucleotidesPairs[iStartPosition1]));
				++iStartPosition1;
				do
				{
					--iStartPosition2;
				}
				while ((0 < iStartPosition2) && (0 == m_tNucleotidesPairs[iStartPosition2]));
				++iStartPosition2;
				do
				{
					++iEndPosition2;
				}
				while (((signed)m_tNucleotidesPairs.size() > iEndPosition2) && (0 == m_tNucleotidesPairs[iEndPosition2]));
				if ((signed)m_tNucleotidesPairs.size() == iEndPosition2)
				{
					// end reached
					--iEndPosition2;
				}
				else
				{
					// end not reached
					iEndPosition2 = m_tNucleotidesPairs[m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].first];
				}
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_3PRIME_FREE_BASES:
			{
				do
				{
					++iEndPosition1;
				}
				while (((signed)m_tNucleotidesPairs.size() > iEndPosition1) && (0 == m_tNucleotidesPairs[iEndPosition1]));
				--iEndPosition1;
				do
				{
					++iEndPosition2;
				}
				while (((signed)m_tNucleotidesPairs.size() > iEndPosition2) && (0 == m_tNucleotidesPairs[iEndPosition2]));
				--iEndPosition2;
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_3PRIME_FREE_BASES_WITH_FIRST_5PRIME_FREE_BASES:
			{
				do
				{
					--iStartPosition1;
				}
				while ((0 < iStartPosition1) && (0 == m_tNucleotidesPairs[iStartPosition1]));
				if (0 == iStartPosition1)
				{
					// begining reached
					++iStartPosition1;
				}
				else
				{
					// begining not reached
					iStartPosition1 = m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].first;
				}
				do
				{
					++iEndPosition1;
				}
				while (((signed)m_tNucleotidesPairs.size() > iEndPosition1) && (0 == m_tNucleotidesPairs[iEndPosition1]));
				--iEndPosition1;
				do
				{
					++iEndPosition2;
				}
				while (((signed)m_tNucleotidesPairs.size() > iEndPosition2) && (0 == m_tNucleotidesPairs[iEndPosition2]));
				--iEndPosition2;
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_NEIGHBOR_FREE_BASES:
			{
				do
				{
					--iStartPosition1;
				}
				while ((0 < iStartPosition1) && (0 == m_tNucleotidesPairs[iStartPosition1]));
				++iStartPosition1;
				do
				{
					--iStartPosition2;
				}
				while ((0 < iStartPosition2) && (0 == m_tNucleotidesPairs[iStartPosition2]));
				++iStartPosition2;
				do
				{
					++iEndPosition1;
				}
				while (((signed)m_tNucleotidesPairs.size() > iEndPosition1) && (0 == m_tNucleotidesPairs[iEndPosition1]));
				--iEndPosition1;
				do
				{
					++iEndPosition2;
				}
				while (((signed)m_tNucleotidesPairs.size() > iEndPosition2) && (0 == m_tNucleotidesPairs[iEndPosition2]));
				--iEndPosition2;
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_HALF_FREE_BASES:
			{
				int iFreeBaseStart = iStartPosition1;
				int iFreeBaseEnd = iEndPosition1;
				do
				{
					--iFreeBaseStart;
				}
				while ((0 < iFreeBaseStart) && (0 == m_tNucleotidesPairs[iFreeBaseStart]));
				++iFreeBaseStart;
				do
				{
					++iFreeBaseEnd;
				}
				while (((signed)m_tNucleotidesPairs.size() > iFreeBaseEnd) && (0 == m_tNucleotidesPairs[iFreeBaseEnd]));
				--iFreeBaseEnd;
				// check if we reached the begining of the structure
				if (1 == iFreeBaseStart)
				{
					iStartPosition1 = 1;
				}
				else
				{
					iStartPosition1 = iStartPosition1 - (int)((float)(iStartPosition1 - iFreeBaseStart)/2. + 0.5);
				}
				iEndPosition1 = iFreeBaseEnd - (int)((float)(iFreeBaseEnd - iEndPosition1)/2. + 0.5);
				iFreeBaseStart = iStartPosition2;
				iFreeBaseEnd = iEndPosition2;
				do
				{
					--iFreeBaseStart;
				}
				while ((0 < iFreeBaseStart) && (0 == m_tNucleotidesPairs[iFreeBaseStart]));
				++iFreeBaseStart;
				do
				{
					++iFreeBaseEnd;
				}
				while (((signed)m_tNucleotidesPairs.size() > iFreeBaseEnd) && (0 == m_tNucleotidesPairs[iFreeBaseEnd]));
				--iFreeBaseEnd;
				iStartPosition2 = iStartPosition2 - (int)((float)(iStartPosition2 - iFreeBaseStart)/2. + 0.5);
				// check if we reached the end of the structure
				if ((signed)m_tNucleotidesPairs.size() == iFreeBaseEnd + 1)
				{
					iEndPosition2 = iFreeBaseEnd;
				}
				else
				{
					iEndPosition2 = iFreeBaseEnd - (int)((float)(iFreeBaseEnd - iEndPosition2)/2. + 0.5);
				}
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_AT_PAIR:
			{
				// no base to add
				break;
			}
			default:
			{
				break;
			}
		}

		// compute stem offset in current structure
		int iStartOffset1 = iStartPosition1 - 1;
		int iStartOffset2 = iStartPosition2 - 1;
		// get stem sub-sequence
		tStemLoop.AppendSequence(*this, iStartOffset1, iEndPosition1 - iStartOffset1);
		tStemLoop.m_iStemSplitPosition = tStemLoop.GetLength() + 1;
		tStemLoop.AppendSequence(*this, iStartOffset2, iEndPosition2 - iStartOffset2);
		iStartOffset2 -= tStemLoop.m_iStemSplitPosition - 1;
		// only one stem
		tStemLoop.m_iStemsCount = 1;
		tStemLoop.m_iStemLoopsCount = 0;
		tStemLoop.m_tStemsAndStemLoopsNature.push_back(val::biocpp::na_secondary_structure::I_STEM_TYPE);
		tStemLoop.m_tStemsAndStemLoopsStartingBases.push_back(std::pair<int, int>(m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].first - iStartOffset1, m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].second - iStartOffset2));

		// extract structural information
		tStemLoop.m_iExternalBasesCount = 0;
		tStemLoop.m_iStackedPairBasesCount = 0;
		tStemLoop.m_iBulgeBasesCount = 0;
		tStemLoop.m_iInternalLoopBasesCount = 0;
		tStemLoop.m_iHairpinLoopBasesCount = 0;
		tStemLoop.m_iMultibrachedLoopBasesCount = 0;
		tStemLoop.m_iPseudoknotBasesCount = 0;
		tStemLoop.m_iKissingHairpinBasesCount = 0;
		tStemLoop.m_iKissingLoopBasesCount = 0;
		tStemLoop.m_iTripleHelixBasesCount = 0;
	
		tStemLoop.m_tNucleotidesConfig.resize(tStemLoop.m_tNucleotidesPairs.size());
		tStemLoop.m_tNucleotidesConfigTertiaryStructure.resize(tStemLoop.m_tNucleotidesPairs.size());
		
		// work on first segment
		for (int iPosition=iStartPosition1; iPosition<=iEndPosition1; iPosition++)
		{
			// check if paired
			if (0 == m_tNucleotidesPairs[iPosition])
			{
				tStemLoop.m_tNucleotidesPairs[iPosition - iStartOffset1] = 0;
			}
			else
			{
				tStemLoop.m_tNucleotidesPairs[iPosition - iStartOffset1] = m_tNucleotidesPairs[iPosition] - iStartOffset2;
			}
			// update counters
			switch (m_tNucleotidesConfig[iPosition])
			{
				case val::biocpp::na_secondary_structure::I_EXTERNAL_BASE:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset1] = val::biocpp::na_secondary_structure::I_EXTERNAL_BASE;
					tStemLoop.m_iExternalBasesCount++;
					break;
				}
				case val::biocpp::na_secondary_structure::I_STACKED_PAIR:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset1] = val::biocpp::na_secondary_structure::I_STACKED_PAIR;
					tStemLoop.m_iStackedPairBasesCount++;
					break;
				}
				case val::biocpp::na_secondary_structure::I_BULGE:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset1] = val::biocpp::na_secondary_structure::I_BULGE;
					tStemLoop.m_iBulgeBasesCount++;
					break;
				}
				case val::biocpp::na_secondary_structure::I_INTERNAL_LOOP:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset1] = val::biocpp::na_secondary_structure::I_INTERNAL_LOOP;
					tStemLoop.m_iInternalLoopBasesCount++;
					break;
				}
				case val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP:
				{
					// replaces hairpin by external bases
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset1] = val::biocpp::na_secondary_structure::I_EXTERNAL_BASE;
					tStemLoop.m_iExternalBasesCount++;
					break;
				}
				case val::biocpp::na_secondary_structure::I_MULTIBRANCHED_LOOP:
				{
					// replaces multibranched by external bases
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset1] = val::biocpp::na_secondary_structure::I_EXTERNAL_BASE;
					tStemLoop.m_iExternalBasesCount++;
					break;
				}
				default:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset1] = val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG;
				}
			}

			// check if tertiary interaction is inside the stem
			if ((val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG != m_tNucleotidesConfigTertiaryStructure[iPosition])
			    && (m_tNucleotidesTertiaryInteractions[iPosition] >= iStartPosition1)
			    && (m_tNucleotidesTertiaryInteractions[iPosition] <= iEndPosition1))
			{
				tStemLoop.m_tNucleotidesTertiaryInteractions[iPosition - iStartOffset1] = m_tNucleotidesTertiaryInteractions[iPosition] - iStartOffset1;
				tStemLoop.m_tNucleotidesConfigTertiaryStructure[iPosition - iStartOffset1] = m_tNucleotidesConfigTertiaryStructure[iPosition];
				// update counters
				switch (m_tNucleotidesConfigTertiaryStructure[iPosition])
				{
					case val::biocpp::na_secondary_structure::I_PSEUDOKNOT:
					{
						tStemLoop.m_iPseudoknotBasesCount++;
						break;
					}
					case val::biocpp::na_secondary_structure::I_KISSING_HAIRPIN:
					{
						tStemLoop.m_iKissingHairpinBasesCount++;
						break;
					}
					case val::biocpp::na_secondary_structure::I_KISSING_LOOP:
					{
						tStemLoop.m_iKissingLoopBasesCount++;
						break;
					}
					case val::biocpp::na_secondary_structure::I_TRIPLE_HELIX:
					{
						tStemLoop.m_iTripleHelixBasesCount++;
						break;
					}
					default:
					{
					}
				}	
			}
			else if ((val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG != m_tNucleotidesConfigTertiaryStructure[iPosition])
			         && (m_tNucleotidesTertiaryInteractions[iPosition] >= iStartPosition2)
			         && (m_tNucleotidesTertiaryInteractions[iPosition] <= iEndPosition2))
			{
				tStemLoop.m_tNucleotidesTertiaryInteractions[iPosition - iStartOffset1] = m_tNucleotidesTertiaryInteractions[iPosition] - iStartOffset2;
				tStemLoop.m_tNucleotidesConfigTertiaryStructure[iPosition - iStartOffset1] = m_tNucleotidesConfigTertiaryStructure[iPosition];
				// update counters
				switch (m_tNucleotidesConfigTertiaryStructure[iPosition])
				{
					case val::biocpp::na_secondary_structure::I_PSEUDOKNOT:
					{
						tStemLoop.m_iPseudoknotBasesCount++;
						break;
					}
					case val::biocpp::na_secondary_structure::I_KISSING_HAIRPIN:
					{
						tStemLoop.m_iKissingHairpinBasesCount++;
						break;
					}
					case val::biocpp::na_secondary_structure::I_KISSING_LOOP:
					{
						tStemLoop.m_iKissingLoopBasesCount++;
						break;
					}
					case val::biocpp::na_secondary_structure::I_TRIPLE_HELIX:
					{
						tStemLoop.m_iTripleHelixBasesCount++;
						break;
					}
					default:
					{
					}
				}	
			}
			else

			{
				tStemLoop.m_tNucleotidesTertiaryInteractions[iPosition - iStartOffset1] = 0;
				tStemLoop.m_tNucleotidesConfigTertiaryStructure[iPosition - iStartOffset1] = val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG;
			}
		}
		// work on second segment
		for (int iPosition=iStartPosition2; iPosition<=iEndPosition2; iPosition++)
		{
			// check if paired
			if (0 == m_tNucleotidesPairs[iPosition])
			{
				tStemLoop.m_tNucleotidesPairs[iPosition - iStartOffset2] = 0;
			}
			else
			{
				tStemLoop.m_tNucleotidesPairs[iPosition - iStartOffset2] = m_tNucleotidesPairs[iPosition] - iStartOffset1;
			}
			// update counters
			switch (m_tNucleotidesConfig[iPosition])
			{
				case val::biocpp::na_secondary_structure::I_EXTERNAL_BASE:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset2] = val::biocpp::na_secondary_structure::I_EXTERNAL_BASE;
					tStemLoop.m_iExternalBasesCount++;
					break;
				}
				case val::biocpp::na_secondary_structure::I_STACKED_PAIR:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset2] = val::biocpp::na_secondary_structure::I_STACKED_PAIR;
					tStemLoop.m_iStackedPairBasesCount++;
					break;
				}
				case val::biocpp::na_secondary_structure::I_BULGE:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset2] = val::biocpp::na_secondary_structure::I_BULGE;
					tStemLoop.m_iBulgeBasesCount++;
					break;
				}
				case val::biocpp::na_secondary_structure::I_INTERNAL_LOOP:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset2] = val::biocpp::na_secondary_structure::I_INTERNAL_LOOP;
					tStemLoop.m_iInternalLoopBasesCount++;
					break;
				}
				case val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset2] = val::biocpp::na_secondary_structure::I_EXTERNAL_BASE;
					tStemLoop.m_iExternalBasesCount++;
					break;
				}
				case val::biocpp::na_secondary_structure::I_MULTIBRANCHED_LOOP:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset2] = val::biocpp::na_secondary_structure::I_EXTERNAL_BASE;
					tStemLoop.m_iExternalBasesCount++;
					break;
				}
				default:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset2] = val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG;
				}
			}

			// check if tertiary interaction is inside the stem
			if ((val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG != m_tNucleotidesConfigTertiaryStructure[iPosition])
			    && (m_tNucleotidesTertiaryInteractions[iPosition] >= iStartPosition1)
			    && (m_tNucleotidesTertiaryInteractions[iPosition] <= iEndPosition1))
			{
				tStemLoop.m_tNucleotidesTertiaryInteractions[iPosition - iStartOffset2] = m_tNucleotidesTertiaryInteractions[iPosition] - iStartOffset1;
				tStemLoop.m_tNucleotidesConfigTertiaryStructure[iPosition - iStartOffset2] = m_tNucleotidesConfigTertiaryStructure[iPosition];
				// update counters
				switch (m_tNucleotidesConfigTertiaryStructure[iPosition])
				{
					case val::biocpp::na_secondary_structure::I_PSEUDOKNOT:
					{
						tStemLoop.m_iPseudoknotBasesCount++;
						break;
					}
					case val::biocpp::na_secondary_structure::I_KISSING_HAIRPIN:
					{
						tStemLoop.m_iKissingHairpinBasesCount++;
						break;
					}
					case val::biocpp::na_secondary_structure::I_KISSING_LOOP:
					{
						tStemLoop.m_iKissingLoopBasesCount++;
						break;
					}
					case val::biocpp::na_secondary_structure::I_TRIPLE_HELIX:
					{
						tStemLoop.m_iTripleHelixBasesCount++;
						break;
					}
					default:
					{
					}
				}	
			}
			else if ((val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG != m_tNucleotidesConfigTertiaryStructure[iPosition])
			         && (m_tNucleotidesTertiaryInteractions[iPosition] >= iStartPosition2)
			         && (m_tNucleotidesTertiaryInteractions[iPosition] <= iEndPosition2))
			{
				tStemLoop.m_tNucleotidesTertiaryInteractions[iPosition - iStartOffset2] = m_tNucleotidesTertiaryInteractions[iPosition] - iStartOffset2;
				tStemLoop.m_tNucleotidesConfigTertiaryStructure[iPosition - iStartOffset2] = m_tNucleotidesConfigTertiaryStructure[iPosition];
				// update counters
				switch (m_tNucleotidesConfigTertiaryStructure[iPosition])
				{
					case val::biocpp::na_secondary_structure::I_PSEUDOKNOT:
					{
						tStemLoop.m_iPseudoknotBasesCount++;
						break;
					}
					case val::biocpp::na_secondary_structure::I_KISSING_HAIRPIN:
					{
						tStemLoop.m_iKissingHairpinBasesCount++;
						break;
					}
					case val::biocpp::na_secondary_structure::I_KISSING_LOOP:
					{
						tStemLoop.m_iKissingLoopBasesCount++;
						break;
					}
					case val::biocpp::na_secondary_structure::I_TRIPLE_HELIX:
					{
						tStemLoop.m_iTripleHelixBasesCount++;
						break;
					}
					default:
					{
					}
				}	
			}
			else
			{
				tStemLoop.m_tNucleotidesTertiaryInteractions[iPosition - iStartOffset2] = 0;
				tStemLoop.m_tNucleotidesConfigTertiaryStructure[iPosition - iStartOffset2] = val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG;
			}
		}
		// set structure as parsed
		tStemLoop.m_bIsStructureParsed = true;
	}
	else // I_STEMLOOP_TYPE
	{
		int iStartPosition = m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].first;
		int iEndPosition = m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].second;
		// stems or stem-loops must start with a pair
		if ((0 == iEndPosition) || (iEndPosition <= iStartPosition))
		{
			throw std::runtime_error("CNASecondaryStructure::GetStemOrStemLoop: " + STR_ERROR_CORRUPTED_STRUCTURE_DATA);
		}
	
		// see where to cut
		switch (iCutMethod)
		{
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_5PRIME_FREE_BASES:
			{
				do
				{
					--iStartPosition;
				}
				while ((0 < iStartPosition) && (0 == m_tNucleotidesPairs[iStartPosition]));
				++iStartPosition;
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_5PRIME_FREE_BASES_WITH_LAST_3PRIME_FREE_BASES:
			{
				do
				{
					--iStartPosition;
				}
				while ((0 < iStartPosition) && (0 == m_tNucleotidesPairs[iStartPosition]));
				++iStartPosition;
				do
				{
					++iEndPosition;
				}
				while (((signed)m_tNucleotidesPairs.size() > iEndPosition) && (0 == m_tNucleotidesPairs[iEndPosition]));
				if ((signed)m_tNucleotidesPairs.size() == iEndPosition)
				{
					--iEndPosition;
				}
				else
				{
					iEndPosition = m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].second;
				}
				break;				
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_3PRIME_FREE_BASES:
			{
				do
				{
					++iEndPosition;
				}
				while (((signed)m_tNucleotidesPairs.size() > iEndPosition) && (0 == m_tNucleotidesPairs[iEndPosition]));
				--iEndPosition;
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_3PRIME_FREE_BASES_WITH_FIRST_5PRIME_FREE_BASES:
			{
				do
				{
					--iStartPosition;
				}
				while ((0 < iStartPosition) && (0 == m_tNucleotidesPairs[iStartPosition]));
				if (0 == iStartPosition)
				{
					++iStartPosition;
				}
				else
				{
					iStartPosition = m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].first;
				}
				do
				{
					++iEndPosition;
				}
				while (((signed)m_tNucleotidesPairs.size() > iEndPosition) && (0 == m_tNucleotidesPairs[iEndPosition]));
				--iEndPosition;
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_NEIGHBOR_FREE_BASES:
			{
				do
				{
					--iStartPosition;
				}
				while ((0 < iStartPosition) && (0 == m_tNucleotidesPairs[iStartPosition]));
				++iStartPosition;
				do
				{
					++iEndPosition;
				}
				while (((signed)m_tNucleotidesPairs.size() > iEndPosition) && (0 == m_tNucleotidesPairs[iEndPosition]));
				--iEndPosition;
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_HALF_FREE_BASES:
			{
				int iFreeBaseStart = iStartPosition;
				int iFreeBaseEnd = iEndPosition;
				do
				{
					--iFreeBaseStart;
				}
				while ((0 < iFreeBaseStart) && (0 == m_tNucleotidesPairs[iFreeBaseStart]));
				++iFreeBaseStart;
				do
				{
					++iFreeBaseEnd;
				}
				while (((signed)m_tNucleotidesPairs.size() > iFreeBaseEnd) && (0 == m_tNucleotidesPairs[iFreeBaseEnd]));
				--iFreeBaseEnd;
				// check if we reached the begining of the structure
				if (1 == iFreeBaseStart)
				{
					iStartPosition = 1;
				}
				else
				{
					iStartPosition = iStartPosition - (int)((float)(iStartPosition - iFreeBaseStart)/2. + 0.5);
				}
				// check if we reached the end of the structure
				if ((signed)m_tNucleotidesPairs.size() == iFreeBaseEnd + 1)
				{
					iEndPosition = iFreeBaseEnd;
				}
				else
				{
					iEndPosition = iFreeBaseEnd - (int)((float)(iFreeBaseEnd - iEndPosition)/2. + 0.5);
				}
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_AT_PAIR:
			{
				// no base to add
				break;
			}
			default:
			{
				break;
			}
		}

		// compute stem-loop offset in current structure
		int iStartOffset = iStartPosition - 1;
		// get stem-loop sub-sequence
		tStemLoop.AppendSequence(*this, iStartOffset, iEndPosition - iStartOffset);
		// only one stem-loop
		tStemLoop.m_iStemsCount = 0;
		tStemLoop.m_iStemLoopsCount = 1;
		tStemLoop.m_tStemsAndStemLoopsNature.push_back(val::biocpp::na_secondary_structure::I_STEMLOOP_TYPE);
		tStemLoop.m_tStemsAndStemLoopsStartingBases.push_back(std::pair<int, int>(m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].first - iStartOffset, m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].second - iStartOffset));
		// extract structural information
		tStemLoop.m_iExternalBasesCount = 0;
		tStemLoop.m_iStackedPairBasesCount = 0;
		tStemLoop.m_iBulgeBasesCount = 0;
		tStemLoop.m_iInternalLoopBasesCount = 0;
		tStemLoop.m_iHairpinLoopBasesCount = 0;
		tStemLoop.m_iMultibrachedLoopBasesCount = 0;
		tStemLoop.m_iPseudoknotBasesCount = 0;
		tStemLoop.m_iKissingHairpinBasesCount = 0;
		tStemLoop.m_iKissingLoopBasesCount = 0;
		tStemLoop.m_iTripleHelixBasesCount = 0;
		tStemLoop.m_iStemSplitPosition = 0;
	
		tStemLoop.m_tNucleotidesConfig.resize(tStemLoop.m_tNucleotidesPairs.size());
		tStemLoop.m_tNucleotidesConfigTertiaryStructure.resize(tStemLoop.m_tNucleotidesPairs.size());
		
		for (int iPosition=iStartPosition; iPosition<=iEndPosition; iPosition++)
		{
			// check if paired
			if (0 == m_tNucleotidesPairs[iPosition])
			{
				tStemLoop.m_tNucleotidesPairs[iPosition - iStartOffset] = 0;
			}
			else
			{
				tStemLoop.m_tNucleotidesPairs[iPosition - iStartOffset] = m_tNucleotidesPairs[iPosition] - iStartOffset;
			}
			// update counters
			switch (m_tNucleotidesConfig[iPosition])
			{
				case val::biocpp::na_secondary_structure::I_EXTERNAL_BASE:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset] = val::biocpp::na_secondary_structure::I_EXTERNAL_BASE;
					tStemLoop.m_iExternalBasesCount++;
					break;
				}
				case val::biocpp::na_secondary_structure::I_STACKED_PAIR:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset] = val::biocpp::na_secondary_structure::I_STACKED_PAIR;
					tStemLoop.m_iStackedPairBasesCount++;
					break;
				}
				case val::biocpp::na_secondary_structure::I_BULGE:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset] = val::biocpp::na_secondary_structure::I_BULGE;
					tStemLoop.m_iBulgeBasesCount++;
					break;
				}
				case val::biocpp::na_secondary_structure::I_INTERNAL_LOOP:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset] = val::biocpp::na_secondary_structure::I_INTERNAL_LOOP;
					tStemLoop.m_iInternalLoopBasesCount++;
					break;
				}
				case val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset] = val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP;
					tStemLoop.m_iHairpinLoopBasesCount++;
					break;
				}
				case val::biocpp::na_secondary_structure::I_MULTIBRANCHED_LOOP:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset] = val::biocpp::na_secondary_structure::I_EXTERNAL_BASE;
					tStemLoop.m_iExternalBasesCount++;
					break;
				}
				default:
				{
					tStemLoop.m_tNucleotidesConfig[iPosition - iStartOffset] = val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG;
				}
			}

			// check if tertiary interaction is inside the stem-loop
			if ((val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG != m_tNucleotidesConfigTertiaryStructure[iPosition]) && (m_tNucleotidesTertiaryInteractions[iPosition] >= iStartPosition) && (m_tNucleotidesTertiaryInteractions[iPosition] <= iEndPosition))
			{
				tStemLoop.m_tNucleotidesTertiaryInteractions[iPosition - iStartOffset] = m_tNucleotidesTertiaryInteractions[iPosition] - iStartOffset;
				tStemLoop.m_tNucleotidesConfigTertiaryStructure[iPosition - iStartOffset] = m_tNucleotidesConfigTertiaryStructure[iPosition];
				// update counters
				switch (m_tNucleotidesConfigTertiaryStructure[iPosition])
				{
					case val::biocpp::na_secondary_structure::I_PSEUDOKNOT:
					{
						tStemLoop.m_iPseudoknotBasesCount++;
						break;
					}
					case val::biocpp::na_secondary_structure::I_KISSING_HAIRPIN:
					{
						tStemLoop.m_iKissingHairpinBasesCount++;
						break;
					}
					case val::biocpp::na_secondary_structure::I_KISSING_LOOP:
					{
						tStemLoop.m_iKissingLoopBasesCount++;
						break;
					}
					case val::biocpp::na_secondary_structure::I_TRIPLE_HELIX:
					{
						tStemLoop.m_iTripleHelixBasesCount++;
						break;
					}
					default:
					{
					}
				}	
			}
			else
			{
				tStemLoop.m_tNucleotidesTertiaryInteractions[iPosition - iStartOffset] = 0;
				tStemLoop.m_tNucleotidesConfigTertiaryStructure[iPosition - iStartOffset] = val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG;
			}
		}
		tStemLoop.m_bIsStructureParsed = true;
	}
	return tStemLoop;
}


/**
@fn int CNASecondaryStructure::GetStemSplitPosition() const
***********************************************************
Returns the position of the begining of the second segment of a stem or 0 for non-stem.

@return int:
 the index of the first base of the second segment of a stem or 0 if the
 structure is not a simple stem.
*/
int CNASecondaryStructure::GetStemSplitPosition() const
{
	TRACEI("int CNASecondaryStructure::GetStemSplitPosition() const");
	return m_iStemSplitPosition;
}


/**
@fn void CNASecondaryStructure::SetStemSplitPosition(int iPosition)
*******************************************************************
Sets the position of the begining of the second segment of a stem.
Throw an exception if the structure is not a stem-loop (before the call) or if
the split position is not on the last pair of the stem-loop or in the hairpin
loop.

@throw std::invalid_argument:\n
if the structure is not a stem-loop or the split position is not valid (not
hairpin loop or second base of the last pair)

*/
void CNASecondaryStructure::SetStemSplitPosition(int iPosition)
{
	TRACEI("void CNASecondaryStructure::SetStemSplitPosition(int iPosition)");
	// check if already parsed
	if (false == m_bIsStructureParsed)
	{
		// parse structure
		ParseStructure();
	}
	// check if we got a simple stem-loop
	if ((0 < GetStemsCount()) || (1 < GetStemLoopsCount()))
	{
		throw std::invalid_argument("ERROR: tried to split a non-stem-loop structure!");
	}
	// check split position
	if ((1 > iPosition) || (iPosition > GetLength()) || ((false == IsInHairpinLoop(iPosition)) && (false == IsInHairpinLoop(iPosition-1)) && (m_tNucleotidesPairs[iPosition-1] != iPosition) && (0 < m_iStackedPairBasesCount)))
	{
		std::ostringstream tOutputStream;
		tOutputStream << "ERROR: tried to split at an invalid position (" << iPosition << ")!";
		throw std::invalid_argument(tOutputStream.str());
	}

	m_iStemSplitPosition = iPosition;
	// now stem-loop 1 is a stem
	m_tStemsAndStemLoopsNature[0] = val::biocpp::na_secondary_structure::I_STEM_TYPE;
	m_iStemLoopsCount = 0;
	m_iStemsCount = 1;

	int iFixConfigPosition = 1;
	m_tStemsAndStemLoopsStartingBases[0].second = 0;
	while ((signed)m_tNucleotidesConfig.size() > iFixConfigPosition)
	{
		if (val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP == m_tNucleotidesConfig[iFixConfigPosition])
		{
			m_tNucleotidesConfig[iFixConfigPosition] = val::biocpp::na_secondary_structure::I_EXTERNAL_BASE;
		}
		// find second base of last pair
		if ((0 == m_tStemsAndStemLoopsStartingBases[0].second)
			&& (val::biocpp::na_secondary_structure::I_STACKED_PAIR == m_tNucleotidesConfig[iFixConfigPosition])
			&& (m_tNucleotidesPairs[iFixConfigPosition] < iFixConfigPosition))
		{
			m_tStemsAndStemLoopsStartingBases[0].second = iFixConfigPosition;
		}
		++iFixConfigPosition;
	}

	m_iExternalBasesCount += m_iHairpinLoopBasesCount;
	m_iHairpinLoopBasesCount = 0;
	
}


/**
@fn CNASecondaryStructure::SSubstructureData CNASecondaryStructure::GetStemOrStemLoopData(int iStemLoopNumber, int iCutMethod) const
************************************************************************************************************************************
For stem, returns index of starting and ending base of the 2 segments.
For stem-loop, returns starting base and ending base of the stem and the index
of the first pair.
*/
CNASecondaryStructure::SSubstructureData CNASecondaryStructure::GetStemOrStemLoopData(int iStemLoopNumber, int iCutMethod) const
{
	TRACEI("CNASecondaryStructure::SSubstructureData CNASecondaryStructure::GetStemOrStemLoopData(int iStemLoopNumber, int iCutMethod) const");
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CNASecondaryStructure::GetStemOrStemLoopData: " + STR_ERROR_NOT_PARSED_STRUCTURE);
	}
	if ((0 >= iStemLoopNumber) || ((signed)m_tStemsAndStemLoopsStartingBases.size() < iStemLoopNumber))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNASecondaryStructure::GetStemOrStemLoopData(" << iStemLoopNumber << ", " << iCutMethod << "): " << STR_ERROR_INVALID_STEM_OR_STEMLOOP << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}

	SSubstructureData tRetData;
	if (val::biocpp::na_secondary_structure::I_STEM_TYPE == m_tStemsAndStemLoopsNature[iStemLoopNumber-1])
	{
		tRetData.iStructureNature = val::biocpp::na_secondary_structure::I_STEM_TYPE;
		tRetData.iStartIndex1 = m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].first;
		tRetData.iEndIndex1 = m_tNucleotidesPairs[m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].second];
		tRetData.iStartIndex2 = m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].second;
		tRetData.iEndIndex2 = m_tNucleotidesPairs[m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].first];
		// stems or stem-loops must start with a pair
		if ((0 == tRetData.iEndIndex1) || (tRetData.iEndIndex1 < tRetData.iStartIndex1)
		    || (0 == tRetData.iEndIndex2) || (tRetData.iEndIndex2 < tRetData.iStartIndex2)
		    || (tRetData.iStartIndex2 <= tRetData.iEndIndex1))
		{
			throw std::runtime_error("CNASecondaryStructure::GetStemOrStemLoopData: " + STR_ERROR_CORRUPTED_STRUCTURE_DATA);
		}

		// see where to cut
		switch (iCutMethod)
		{
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_5PRIME_FREE_BASES:
			{
				do
				{
					--tRetData.iStartIndex1;
				}
				while ((0 < tRetData.iStartIndex1) && (0 == m_tNucleotidesPairs[tRetData.iStartIndex1]));
				++tRetData.iStartIndex1;
				do
				{
					--tRetData.iStartIndex2;
				}
				while ((0 < tRetData.iStartIndex2) && (0 == m_tNucleotidesPairs[tRetData.iStartIndex2]));
				++tRetData.iStartIndex2;
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_5PRIME_FREE_BASES_WITH_LAST_3PRIME_FREE_BASES:
			{
				do
				{
					--tRetData.iStartIndex1;
				}
				while ((0 < tRetData.iStartIndex1) && (0 == m_tNucleotidesPairs[tRetData.iStartIndex1]));
				++tRetData.iStartIndex1;
				do
				{
					--tRetData.iStartIndex2;
				}
				while ((0 < tRetData.iStartIndex2) && (0 == m_tNucleotidesPairs[tRetData.iStartIndex2]));
				++tRetData.iStartIndex2;
				do
				{
					++tRetData.iEndIndex2;
				}
				while (((signed)m_tNucleotidesPairs.size() > tRetData.iEndIndex2) && (0 == m_tNucleotidesPairs[tRetData.iEndIndex2]));
				if ((signed)m_tNucleotidesPairs.size() == tRetData.iEndIndex2)
				{
					// end reached
					--tRetData.iEndIndex2;
				}
				else
				{
					// end not reached
					tRetData.iEndIndex2 = m_tNucleotidesPairs[m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].first];
				}
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_3PRIME_FREE_BASES:
			{
				do
				{
					++tRetData.iEndIndex1;
				}
				while (((signed)m_tNucleotidesPairs.size() > tRetData.iEndIndex1) && (0 == m_tNucleotidesPairs[tRetData.iEndIndex1]));
				--tRetData.iEndIndex1;
				do
				{
					++tRetData.iEndIndex2;
				}
				while (((signed)m_tNucleotidesPairs.size() > tRetData.iEndIndex2) && (0 == m_tNucleotidesPairs[tRetData.iEndIndex2]));
				--tRetData.iEndIndex2;
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_3PRIME_FREE_BASES_WITH_FIRST_5PRIME_FREE_BASES:
			{
				do
				{
					--tRetData.iStartIndex1;
				}
				while ((0 < tRetData.iStartIndex1) && (0 == m_tNucleotidesPairs[tRetData.iStartIndex1]));
				if (0 == tRetData.iStartIndex1)
				{
					// begining reached
					++tRetData.iStartIndex1;
				}
				else
				{
					// begining not reached
					tRetData.iStartIndex1 = m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].first;
				}
				do
				{
					++tRetData.iEndIndex1;
				}
				while (((signed)m_tNucleotidesPairs.size() > tRetData.iEndIndex1) && (0 == m_tNucleotidesPairs[tRetData.iEndIndex1]));
				--tRetData.iEndIndex1;
				do
				{
					++tRetData.iEndIndex2;
				}
				while (((signed)m_tNucleotidesPairs.size() > tRetData.iEndIndex2) && (0 == m_tNucleotidesPairs[tRetData.iEndIndex2]));
				--tRetData.iEndIndex2;
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_NEIGHBOR_FREE_BASES:
			{
				do
				{
					--tRetData.iStartIndex1;
				}
				while ((0 < tRetData.iStartIndex1) && (0 == m_tNucleotidesPairs[tRetData.iStartIndex1]));
				++tRetData.iStartIndex1;
				do
				{
					--tRetData.iStartIndex2;
				}
				while ((0 < tRetData.iStartIndex2) && (0 == m_tNucleotidesPairs[tRetData.iStartIndex2]));
				++tRetData.iStartIndex2;
				do
				{
					++tRetData.iEndIndex1;
				}
				while (((signed)m_tNucleotidesPairs.size() > tRetData.iEndIndex1) && (0 == m_tNucleotidesPairs[tRetData.iEndIndex1]));
				--tRetData.iEndIndex1;
				do
				{
					++tRetData.iEndIndex2;
				}
				while (((signed)m_tNucleotidesPairs.size() > tRetData.iEndIndex2) && (0 == m_tNucleotidesPairs[tRetData.iEndIndex2]));
				--tRetData.iEndIndex2;
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_HALF_FREE_BASES:
			{
				int iFreeBaseStart = tRetData.iStartIndex1;
				int iFreeBaseEnd = tRetData.iEndIndex1;
				do
				{
					--iFreeBaseStart;
				}
				while ((0 < iFreeBaseStart) && (0 == m_tNucleotidesPairs[iFreeBaseStart]));
				++iFreeBaseStart;
				do
				{
					++iFreeBaseEnd;
				}
				while (((signed)m_tNucleotidesPairs.size() > iFreeBaseEnd) && (0 == m_tNucleotidesPairs[iFreeBaseEnd]));
				--iFreeBaseEnd;
				// check if we reached the begining of the structure
				if (1 == iFreeBaseStart)
				{
					tRetData.iStartIndex1 = 1;
				}
				else
				{
					tRetData.iStartIndex1 = tRetData.iStartIndex1 - (int)((float)(tRetData.iStartIndex1 - iFreeBaseStart)/2. + 0.5);
				}
				tRetData.iEndIndex1 = iFreeBaseEnd - (int)((float)(iFreeBaseEnd - tRetData.iEndIndex1)/2. + 0.5);
				iFreeBaseStart = tRetData.iStartIndex2;
				iFreeBaseEnd = tRetData.iEndIndex2;
				do
				{
					--iFreeBaseStart;
				}
				while ((0 < iFreeBaseStart) && (0 == m_tNucleotidesPairs[iFreeBaseStart]));
				++iFreeBaseStart;
				do
				{
					++iFreeBaseEnd;
				}
				while (((signed)m_tNucleotidesPairs.size() > iFreeBaseEnd) && (0 == m_tNucleotidesPairs[iFreeBaseEnd]));
				--iFreeBaseEnd;
				tRetData.iStartIndex2 = tRetData.iStartIndex2 - (int)((float)(tRetData.iStartIndex2 - iFreeBaseStart)/2. + 0.5);
				// check if we reached the end of the structure
				if ((signed)m_tNucleotidesPairs.size() == iFreeBaseEnd + 1)
				{
					tRetData.iEndIndex2 = iFreeBaseEnd;
				}
				else
				{
					tRetData.iEndIndex2 = iFreeBaseEnd - (int)((float)(iFreeBaseEnd - tRetData.iEndIndex2)/2. + 0.5);
				}
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_AT_PAIR:
			{
				// no base to add
				break;
			}
			default:
			{
				break;
			}
		}
	}
	else // I_STEMLOOP_TYPE
	{
		tRetData.iStructureNature = val::biocpp::na_secondary_structure::I_STEMLOOP_TYPE;
		tRetData.iStartIndex1 = tRetData.iStartIndex2 = m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].first;
		tRetData.iEndIndex1 = tRetData.iEndIndex2 = m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].second;
		// stems or stem-loops must start with a pair
		if ((0 == tRetData.iEndIndex1) || (tRetData.iEndIndex1 <= tRetData.iStartIndex1))
		{
			throw std::runtime_error("CNASecondaryStructure::GetStemOrStemLoop: " + STR_ERROR_CORRUPTED_STRUCTURE_DATA);
		}
	
		// see where to cut
		switch (iCutMethod)
		{
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_5PRIME_FREE_BASES:
			{
				do
				{
					--tRetData.iStartIndex1;
				}
				while ((0 < tRetData.iStartIndex1) && (0 == m_tNucleotidesPairs[tRetData.iStartIndex1]));
				++tRetData.iStartIndex1;
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_5PRIME_FREE_BASES_WITH_LAST_3PRIME_FREE_BASES:
			{
				do
				{
					--tRetData.iStartIndex1;
				}
				while ((0 < tRetData.iStartIndex1) && (0 == m_tNucleotidesPairs[tRetData.iStartIndex1]));
				++tRetData.iStartIndex1;
				do
				{
					++tRetData.iEndIndex1;
				}
				while (((signed)m_tNucleotidesPairs.size() > tRetData.iEndIndex1) && (0 == m_tNucleotidesPairs[tRetData.iEndIndex1]));
				if ((signed)m_tNucleotidesPairs.size() == tRetData.iEndIndex1)
				{
					--tRetData.iEndIndex1;
				}
				else
				{
					tRetData.iEndIndex1 = m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].second;
				}
				break;				
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_3PRIME_FREE_BASES:
			{
				do
				{
					++tRetData.iEndIndex1;
				}
				while (((signed)m_tNucleotidesPairs.size() > tRetData.iEndIndex1) && (0 == m_tNucleotidesPairs[tRetData.iEndIndex1]));
				--tRetData.iEndIndex1;
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_3PRIME_FREE_BASES_WITH_FIRST_5PRIME_FREE_BASES:
			{
				do
				{
					--tRetData.iStartIndex1;
				}
				while ((0 < tRetData.iStartIndex1) && (0 == m_tNucleotidesPairs[tRetData.iStartIndex1]));
				if (0 == tRetData.iStartIndex1)
				{
					++tRetData.iStartIndex1;
				}
				else
				{
					tRetData.iStartIndex1 = m_tStemsAndStemLoopsStartingBases[iStemLoopNumber-1].first;
				}
				do
				{
					++tRetData.iEndIndex1;
				}
				while (((signed)m_tNucleotidesPairs.size() > tRetData.iEndIndex1) && (0 == m_tNucleotidesPairs[tRetData.iEndIndex1]));
				--tRetData.iEndIndex1;
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_NEIGHBOR_FREE_BASES:
			{
				do
				{
					--tRetData.iStartIndex1;
				}
				while ((0 < tRetData.iStartIndex1) && (0 == m_tNucleotidesPairs[tRetData.iStartIndex1]));
				++tRetData.iStartIndex1;
				do
				{
					++tRetData.iEndIndex1;
				}
				while (((signed)m_tNucleotidesPairs.size() > tRetData.iEndIndex1) && (0 == m_tNucleotidesPairs[tRetData.iEndIndex1]));
				--tRetData.iEndIndex1;
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_INCLUDING_HALF_FREE_BASES:
			{
				int iFreeBaseStart = tRetData.iStartIndex1;
				int iFreeBaseEnd = tRetData.iEndIndex1;
				do
				{
					--iFreeBaseStart;
				}
				while ((0 < iFreeBaseStart) && (0 == m_tNucleotidesPairs[iFreeBaseStart]));
				++iFreeBaseStart;
				do
				{
					++iFreeBaseEnd;
				}
				while (((signed)m_tNucleotidesPairs.size() > iFreeBaseEnd) && (0 == m_tNucleotidesPairs[iFreeBaseEnd]));
				--iFreeBaseEnd;
				// check if we reached the begining of the structure
				if (1 == iFreeBaseStart)
				{
					tRetData.iStartIndex1 = 1;
				}
				else
				{
					tRetData.iStartIndex1 = tRetData.iStartIndex1 - (int)((float)(tRetData.iStartIndex1 - iFreeBaseStart)/2. + 0.5);
				}
				// check if we reached the end of the structure
				if ((signed)m_tNucleotidesPairs.size() == iFreeBaseEnd + 1)
				{
					tRetData.iEndIndex1 = iFreeBaseEnd;
				}
				else
				{
					tRetData.iEndIndex1 = iFreeBaseEnd - (int)((float)(iFreeBaseEnd - tRetData.iEndIndex1)/2. + 0.5);
				}
				break;
			}
			case val::biocpp::na_secondary_structure::I_CUT_AT_PAIR:
			{
				// no base to add
				break;
			}
			default:
			{
				break;
			}
		}
	}
	return tRetData;
}





//##typedef val::biocpp::tools::CTreeNode<int> CIntegerTreeNode;

/*
void PrintRecursive( CHelixIndexTreeNode::CTreeNodePointer treeNode )
{
	CHelixIndexTreeNode::CTreeNodePointerList suffixList_Children = treeNode->GetDirectChildrenList();
	CHelixIndexTreeNode::CTreeNodePointerList::iterator tIter1 = suffixList_Children.begin();

	while (tIter1 != suffixList_Children.end())
    {

		std::cout << (*tIter1)->GetObject().ToString() << ","; 
		if ((*tIter1)->GetDirectChildrenCount()==2)    // knot m
		{
			std::cout << "(,";
			PrintRecursive(*tIter1);
			std::cout << "),";
		}
		else                                           // leaf
		{
		    PrintRecursive(*tIter1);
		}
		
		tIter1++;
	}
}
* */
typedef val::biocpp::tools::CTreeNode<hishape::comp::CHelixIndex> CHelixIndexTreeNode;
/**
Initializes nodes information (call it before using nodes info methods)
*/
CHelixIndexTreeNode::CTreeNodePointer CNASecondaryStructure::ParseStructureWithTree()
{
	TRACEI("void CNASecondaryStructure::Structure2Hishape()");
	// check if already parsed
	//if (true == m_bIsStructureParsed)
	//{return;}

	m_iExternalBasesCount = 0;
	m_iStackedPairBasesCount = 0;
	m_iBulgeBasesCount = 0;
	m_iInternalLoopBasesCount = 0;
	m_iHairpinLoopBasesCount = 0;
	m_iMultibrachedLoopBasesCount = 0;
	m_iPseudoknotBasesCount = 0;
	m_iKissingHairpinBasesCount = 0;
	m_iKissingLoopBasesCount = 0;
	m_iTripleHelixBasesCount = 0;
	m_iStemSplitPosition = 0;
	
//	m_tNucleotidesConfig.clear();
	m_tNucleotidesConfig.resize(m_tNucleotidesPairs.size());
//	m_tNucleotidesConfigTertiaryStructure.clear();
	m_tNucleotidesConfigTertiaryStructure.resize(m_tNucleotidesPairs.size());
	for (unsigned int i=0; i < m_tNucleotidesPairs.size(); i++)
	{
		m_tNucleotidesConfig[i] = m_tNucleotidesConfigTertiaryStructure[i] = val::biocpp::na_secondary_structure::I_UNKNOWN_CONFIG;
	}

	m_iStemsCount = 0;
	m_iStemLoopsCount = 0;
	m_tStemsAndStemLoopsNature.clear();
	m_tStemsAndStemLoopsStartingBases.clear();

	// index1 = base_5'
	int iIndex1 = 1;
	std::deque<int> tStemsToParseStack;
	
    std::deque<int> tStemIndexDeque;
    std::deque<CHelixIndexTreeNode::CTreeNodePointer> tStemPointerDeque;


    //// std::cout << "m_tNucleotidesPairs[5]=" << m_tNucleotidesPairs[5] << std::endl;
    //// std::cout << "m_tNucleotidesPairs[15]=" << m_tNucleotidesPairs[15] << std::endl;
	// the size is the length of the sequence, if it is paired (4,11) or (5,0), (external, Null)
	
	
	hishape::comp::CHelixIndex tHelixIndex(0.0f, 'N');
	CHelixIndexTreeNode::CTreeNodePointer ptRoot = new CHelixIndexTreeNode(tHelixIndex);
	
	// it returns the number of external bases
	while (iIndex1 < (signed)m_tNucleotidesPairs.size())
	{
		// in first round, it finds all external bases
		while ((iIndex1 < (signed)m_tNucleotidesPairs.size()) && (0 == m_tNucleotidesPairs[iIndex1]))
		{
			m_tNucleotidesConfig[iIndex1] = val::biocpp::na_secondary_structure::I_EXTERNAL_BASE;
			++m_iExternalBasesCount;
			// index1++
			++iIndex1;
		}
		/*
		 * m_iExternalBasesCount=0
			iIndex1=1
			m_iExternalBasesCount=4
			iIndex1=19
			m_iExternalBasesCount=6
			iIndex1=151
			* */
			
		// in second round, 1--this while-->14--first while-->19--this while-->148
		while ((iIndex1 < (signed)m_tNucleotidesPairs.size()) && (0 != m_tNucleotidesPairs[iIndex1]))
		{
			if (m_tNucleotidesPairs[iIndex1] > iIndex1)
			{
				tStemIndexDeque.push_back(iIndex1);
				tStemPointerDeque.push_back(ptRoot);
				iIndex1 = m_tNucleotidesPairs[iIndex1] + 1;
			}
			else
			{
				m_tNucleotidesConfig[iIndex1] = val::biocpp::na_secondary_structure::I_EXTERNAL_BASE;
				m_tNucleotidesConfig[m_tNucleotidesPairs[iIndex1]] = val::biocpp::na_secondary_structure::I_EXTERNAL_BASE;
				m_iExternalBasesCount += 2;
				// marquer index1 et index1_paire comme I_PSEUDOKNOT
				m_tNucleotidesConfigTertiaryStructure[iIndex1] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
				m_tNucleotidesConfigTertiaryStructure[m_tNucleotidesPairs[iIndex1]] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
				m_iPseudoknotBasesCount += 2;
				// index1++
				++iIndex1;
				
			}
		}
		
	}
	
	
	CHelixIndexTreeNode::CTreeNodePointer currentNode = NULL;
	//
	int iIndex = 0;
	// tant_que pile_paires non vide
	while (false == tStemIndexDeque.empty())  // only two elements, later on, it should use deque
	{
		// deque first input first output
		int iStemStart = iIndex1 = tStemIndexDeque.front();
		tStemIndexDeque.pop_front();
		CHelixIndexTreeNode::CTreeNodePointer parentNode = tStemPointerDeque.front();
		tStemPointerDeque.pop_front();

	    
		//##currentNode = ptRoot->AddChild(iIndex1);
	    hishape::comp::CHelixIndex tHelixIndex(iIndex1*1.0f, 'N');
	    currentNode = parentNode->AddChild(tHelixIndex);
								
		// index2 is the pair position of index1 ==> (index1,index2)
		int iIndex2 = m_tNucleotidesPairs[iIndex1];

 	
		while (iIndex1 < iIndex2)
		{
	        std::ostringstream outs;    
			do
			{
				// from (1,14) ==> if (paired) further
				m_tNucleotidesConfig[iIndex1] = val::biocpp::na_secondary_structure::I_STACKED_PAIR;
				m_tNucleotidesConfig[iIndex2] = val::biocpp::na_secondary_structure::I_STACKED_PAIR;
				m_iStackedPairBasesCount += 2;
				++iIndex1;
				iIndex2--;
			  // if next (++iIndex1,index2--) are also paired, coutining
			} while ((iIndex2 == m_tNucleotidesPairs[iIndex1]) && (iIndex1 < iIndex2));
			
			// index1_start = index1
			int iStartIndex1 = iIndex1;
			// index2_start = index2
			int iStartIndex2 = iIndex2;
			
			//std::cout << "case0: iStartIndex1=" << iStartIndex1 << "iStartIndex2=" << iStartIndex2 << std::endl;			
			std::deque<int> t5PrimeStack, t3PrimeStack;
			if (iIndex1 > iIndex2)
			{
				//std::cout << "case1: index1=" << iIndex1 << "index2=" << iIndex2 << std::endl;

				m_tStemsAndStemLoopsNature.push_back(val::biocpp::na_secondary_structure::I_STEMLOOP_TYPE);
				m_tStemsAndStemLoopsStartingBases.push_back(std::pair<int, int>(iStemStart, m_tNucleotidesPairs[iStemStart]));
				++m_iStemLoopsCount;
			}
			else if (iIndex1 == iIndex2)
			{
				//std::cout << "case2: index1=" << iIndex1 << "index2=" << iIndex2 << std::endl;
				
				// marquer index1 comme I_HAIRPIN_LOOP
				m_tNucleotidesConfig[iIndex1] = val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP;
				++m_iHairpinLoopBasesCount;
				if (0 != m_tNucleotidesPairs[iIndex1])
				{
					//+FIXME: retirer l'appariement et le mettre dans les interactions tertiaires
					// marquer index1 comme I_PSEUDOKNOT
					m_tNucleotidesConfigTertiaryStructure[iIndex1] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
					++m_iPseudoknotBasesCount;
				// fin_si
				}
				// marquer stem_start comme un d¿but de stem-loop I_STEMLOOP_TYPE
				m_tStemsAndStemLoopsNature.push_back(val::biocpp::na_secondary_structure::I_STEMLOOP_TYPE);
				m_tStemsAndStemLoopsStartingBases.push_back(std::pair<int, int>(iStemStart, m_tNucleotidesPairs[iStemStart]));
				++m_iStemLoopsCount;
			}
			else
			{
				//std::cout << "case3: index1=" << iIndex1 << "index2=" << iIndex2 << std::endl;
				// tant_que index1 <= index2 et que (index1 non appari¿e ou que index1_paire < index1_start ou que index1_paire > index2_start)
				while ((iIndex1 <= iIndex2) && ((0 == m_tNucleotidesPairs[iIndex1]) || (m_tNucleotidesPairs[iIndex1] < iStartIndex1) || (m_tNucleotidesPairs[iIndex1] > iStartIndex2)))
				{
					// empiler index1 dans pile_5'
					t5PrimeStack.push_back(iIndex1);
					// index1++
					++iIndex1;
				// fin_tant_que
				}
				

				// tant_que index1 < index2 et que (index2 non appari¿e ou que index2_paire < index1_start ou que index2_paire > index2_start)
				if (iIndex1 < iIndex2)
				{
					//##std::cout << "case3_1: index1=" << iIndex1 << "index2=" << iIndex2 << std::endl;
					while ((0 == m_tNucleotidesPairs[iIndex2]) || (m_tNucleotidesPairs[iIndex2] < iStartIndex1) || (m_tNucleotidesPairs[iIndex2] > iStartIndex2))
					{
						// empiler index2 dans pile_3'
						t3PrimeStack.push_back(iIndex2);
						// index2--
						--iIndex2;
					// fin_tant_que
					}
				}
				// si index1 > index2 // aucune paire rencontr¿e de index1 vers index2
				if (iIndex1 > iIndex2)
				{
					//##std::cout << "case3_2: index1=" << iIndex1 << "index2=" << iIndex2 << std::endl;
					// pour chaque base de la pile_5'
					while (false == t5PrimeStack.empty())
					{
						// index = depiler pile_5'
						iIndex = t5PrimeStack.back();
						t5PrimeStack.pop_back();
						//  marquer index comme I_HAIRPIN_LOOP
						m_tNucleotidesConfig[iIndex] = val::biocpp::na_secondary_structure::I_HAIRPIN_LOOP;
						++m_iHairpinLoopBasesCount;
						//##hishape::comp::CHelixIndex tHelixIndex((iIndex1+iIndex2)/2.0f, 'h');
				        //##currentNode = ptRoot->AddChild(new CHelixIndexTreeNode(tHelixIndex));
						
						if (0 != m_tNucleotidesPairs[iIndex])
						{
							// marquer index comme I_PSEUDOKNOT
							m_tNucleotidesConfigTertiaryStructure[iIndex] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
							++m_iPseudoknotBasesCount;
						// fin_si
						}
					// fin_pour
					}
					// marquer stem_start comme un d¿but de stem-loop
					m_tStemsAndStemLoopsNature.push_back(val::biocpp::na_secondary_structure::I_STEMLOOP_TYPE);
					m_tStemsAndStemLoopsStartingBases.push_back(std::pair<int, int>(iStemStart, m_tNucleotidesPairs[iStemStart]));
					++m_iStemLoopsCount;
				}
				else
				{
                    //##std::cout << "case3_3: (index1,index2)=" << iIndex2 << std::endl;
					if (m_tNucleotidesPairs[iIndex1] == iIndex2)
					{
						
						//  t5PrimeStack.empty + t3PrimeStack.empty
						if ((false == t5PrimeStack.empty()) && (false == t3PrimeStack.empty()))
						{
							// d¿piler et marquer chaque base de pile_5' et de pile_3' comme I_INTERNAL_LOOP
							while (false == t5PrimeStack.empty())
							{
								// index = depiler pile_5'
								iIndex = t5PrimeStack.back();
								t5PrimeStack.pop_back();
								// marquer index comme I_BULGE
								m_tNucleotidesConfig[iIndex] = val::biocpp::na_secondary_structure::I_INTERNAL_LOOP;
								++m_iInternalLoopBasesCount;
								// si index appari¿e
								if (0 != m_tNucleotidesPairs[iIndex])
								{
									// marquer index comme I_PSEUDOKNOT
									m_tNucleotidesConfigTertiaryStructure[iIndex] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
									++m_iPseudoknotBasesCount;
								// fin_si
								}
							}
							while (false == t3PrimeStack.empty())
							{
								// index = depiler pile_3'
								iIndex = t3PrimeStack.back();
								t3PrimeStack.pop_back();
								// marquer index comme I_BULGE
								m_tNucleotidesConfig[iIndex] = val::biocpp::na_secondary_structure::I_INTERNAL_LOOP;
								++m_iInternalLoopBasesCount;
								// si index appari¿e
								if (0 != m_tNucleotidesPairs[iIndex])
								{
									// marquer index comme I_PSEUDOKNOT
									m_tNucleotidesConfigTertiaryStructure[iIndex] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
									++m_iPseudoknotBasesCount;
								// fin_si
								}
							}
						}
						// t5PrimeStack.empty || t3PrimeStack.empty
						else
						{
							// pour chaque base de la pile_5'
							while (false == t5PrimeStack.empty())
							{
								// index = depiler pile_5'
								iIndex = t5PrimeStack.back();
								t5PrimeStack.pop_back();
								// marquer index comme I_BULGE
								m_tNucleotidesConfig[iIndex] = val::biocpp::na_secondary_structure::I_BULGE;
								++m_iBulgeBasesCount;
								// si index appari¿e
								if (0 != m_tNucleotidesPairs[iIndex])
								{
									// marquer index comme I_PSEUDOKNOT
									m_tNucleotidesConfigTertiaryStructure[iIndex] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
									++m_iPseudoknotBasesCount;
								// fin_si
								}
							// fin_pour
							}
							// pour chaque base de la pile_3'
							while (false == t3PrimeStack.empty())
							{
								// index = depiler pile_3'
								iIndex = t3PrimeStack.back();
								t3PrimeStack.pop_back();
								// marquer index comme I_BULGE
								m_tNucleotidesConfig[iIndex] = val::biocpp::na_secondary_structure::I_BULGE;
								++m_iBulgeBasesCount;
								// si index appari¿e
								if (0 != m_tNucleotidesPairs[iIndex])
								{
									// marquer index comme I_PSEUDOKNOT
									m_tNucleotidesConfigTertiaryStructure[iIndex] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
									++m_iPseudoknotBasesCount;
								// fin_si
								}
							// fin_pour
							}
						// fin_si
						}
					}
					//std::cout << "case3_3: (index1,index2)!=" << iIndex2 << std::endl;
					else
					{
						while (iIndex1 < iIndex2)
						{    
							// si index1_paire < index1
							if (m_tNucleotidesPairs[iIndex1] < iIndex1)
							{
								// marquer index1 comme I_MULTIBRANCHED_LOOP
								m_tNucleotidesConfig[iIndex1] = val::biocpp::na_secondary_structure::I_MULTIBRANCHED_LOOP;
								++m_iMultibrachedLoopBasesCount;
								// marquer index1 comme I_PSEUDOKNOT
								m_tNucleotidesConfigTertiaryStructure[iIndex1] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
								++m_iPseudoknotBasesCount;
								// index1++
								++iIndex1;
							}
							// sinon
							else
							{
								// empiler index1 dans pile_paires
								tStemIndexDeque.push_back(iIndex1);
								//######################################
								//## save currentNode as parent point ##
								//######################################
                                tStemPointerDeque.push_back(currentNode);  // add pointer of parent
                                				
								// index1 = index1_paire + 1
								iIndex1 = m_tNucleotidesPairs[iIndex1] + 1;
							// fin_si
							}
							// tant_que index1 non appari¿e et que index1 < index2
							while ((0 == m_tNucleotidesPairs[iIndex1]) && (iIndex1 < iIndex2))
							{
								// marquer index1 comme I_MULTIBRANCHED_LOOP
								m_tNucleotidesConfig[iIndex1] = val::biocpp::na_secondary_structure::I_MULTIBRANCHED_LOOP;
								++m_iMultibrachedLoopBasesCount;
								// index1++
								++iIndex1;
							// fin_tant_que
							}
						// fin_tant_que
						}
						// pour chaque base de la pile_5'
						while (false == t5PrimeStack.empty())
						{
							// index = depiler pile_5'
							iIndex = t5PrimeStack.back();
							t5PrimeStack.pop_back();
							// marquer index comme I_MULTIBRANCHED_LOOP
							m_tNucleotidesConfig[iIndex] = val::biocpp::na_secondary_structure::I_MULTIBRANCHED_LOOP;
							++m_iMultibrachedLoopBasesCount;
							// si index appari¿e
							if (0 != m_tNucleotidesPairs[iIndex])
							{
								// marquer index comme I_PSEUDOKNOT
								m_tNucleotidesConfigTertiaryStructure[iIndex] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
								++m_iPseudoknotBasesCount;
							// fin_si
							}
						// fin_pour
						}
						// pour chaque base de la pile_3'
						while (false == t3PrimeStack.empty())
						{
							// index = depiler pile_3'
							iIndex = t3PrimeStack.back();
							t3PrimeStack.pop_back();
							// marquer index comme I_MULTIBRANCHED_LOOP
							m_tNucleotidesConfig[iIndex] = val::biocpp::na_secondary_structure::I_MULTIBRANCHED_LOOP;
							++m_iMultibrachedLoopBasesCount;
							// si index appari¿e
							if (0 != m_tNucleotidesPairs[iIndex])
							{
								// marquer index comme I_PSEUDOKNOT
								m_tNucleotidesConfigTertiaryStructure[iIndex] = val::biocpp::na_secondary_structure::I_PSEUDOKNOT;
								++m_iPseudoknotBasesCount;
							// fin_si
							}
						// fin_pour
						}
						// marquer stem_start comme un d¿but de stem et index2_start + 1 comme une fin de stem
						m_tStemsAndStemLoopsNature.push_back(val::biocpp::na_secondary_structure::I_STEM_TYPE);
						m_tStemsAndStemLoopsStartingBases.push_back(std::pair<int, int>(iStemStart, iStartIndex2 + 1));
						++m_iStemsCount;
					// fin_si
					}
				// fin_si
				}
			// fin_si
			}
		// fin_tant_que
		}
	// fin_tant_que
	}		
	m_bIsStructureParsed = true;
	return ptRoot;
} 


}; // namespace biocpp
}; // namespace val
