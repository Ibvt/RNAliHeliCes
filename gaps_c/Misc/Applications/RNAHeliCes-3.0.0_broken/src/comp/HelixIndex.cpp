//*****************************************************************************
#include <iostream>
#include <sstream>
#include "HelixIndex.hpp"
#include "DebugTools.hpp"

namespace hishape
{
namespace comp
{
/**
Default Constructor
*******************
Creates a new initilized instance of @c CHelixIndex.

@param const nucleic_acids::SNucleotideCode &tFirstNucleotide:\n
 first nucleotide code value.\n
 Default value: @c nucleic_acids::Unknown

@param const nucleic_acids::SNucleotideCode &tSecondNucleotide:\n
 second nucleotide code value.\n
 Default value: @c nucleic_acids::Unknown

@return CHelixIndex:
 a new initilized instance of CHelixIndex.

@see nucleic_acids::SNucleotideCode
@see nucleic_acids::Null
*/
CHelixIndex::CHelixIndex(const float hiNumber, const char hiType/*, const int firstPositionOfPair*/)
{
	//TRACEI("CHelixIndex::CHelixIndex(const nucleic_acids::SNucleotideCode &tNucleotide, const nucleic_acids::SNucleotideCode &tNucleotide)");
	m_hiNumber = hiNumber;
	m_hiType = hiType;
	//m_firstPositionOfPair = firstPositionOfPair;
	////std::cout << "CHelixIndex::CHelixIndex(" << hiNumber << "," << hiType << /*"," << firstPositionOfPair <<*/ ")" << std::endl;
}
CHelixIndex::CHelixIndex()
{
	//TRACEI("CHelixIndex::CHelixIndex(const nucleic_acids::SNucleotideCode &tNucleotide, const nucleic_acids::SNucleotideCode &tNucleotide)");
	m_hiNumber = 0.0f;
	m_hiType = '_';
	////std::cout << "CHelixIndex::CHelixIndex(void)" << std::endl;
}


/**
Destructor
**********
A standard virtual destructor.
*/
CHelixIndex::~CHelixIndex(void)
{
	TRACEI("CHelixIndex::~CHelixIndex(void)");
}




float CHelixIndex::GetHiNumber() const
{
    //TRACEI("int CHelixIndex::GetIntCode() const");
    return m_hiNumber;
}
//!
void CHelixIndex::SetHiNumber(const float hiNumber)
{
    m_hiNumber = hiNumber;
}

char CHelixIndex::GetHiType() const
{
    //TRACEI("int CHelixIndex::GetIntCode() const");
    return m_hiType;
}
//!
void CHelixIndex::SetHiType(const char &hiType)
{
    m_hiType = hiType;
    ////std::cout << "SetHiType(" << hiType << ")" << std::endl;
}

/**
@fn int CHelixIndex::GetIntCode() const
*******************************************
Returns the integer codes of the two nucleotide codes into the same integer.
The high part is the first nucleotide, the low part is the second nucleotide.

@return int:\n
 the pair integer code.
*/
int CHelixIndex::GetIntCode() const
{
	TRACEI("int CHelixIndex::GetIntCode() const");
	//return (int)m_hiNumber;//m_iIntCode;
	int intCode = -1;
	switch (m_hiType)
	{
	    case 'N':
	      intCode = N;
	      break;
	    case 'i':
	      intCode = I;
	      break;
	    case 'b':
	      intCode = B;
	      break;
	    case 'h':
	      intCode = H;
	      break;
	    case 'm':
	      intCode = M;
	      break; 
	    default :
	      std::cout << "Not valid hishape type!" << std::endl;
	}
	return intCode;
}

int CHelixIndex::GetGapCode() const
{
	TRACEI("int CHelixIndex::GetGapCode() const");
	return INDEL;
}

/*
int CHelixIndex::GetFristPositionOfPair() const
{
    return m_firstPositionOfPair;
}

void CHelixIndex::SetFristPositionOfPair(const int firstPositionOfPair)
{
    m_firstPositionOfPair = firstPositionOfPair;
}*/

/**
@fn std::string CHelixIndex::ToString() const
*************************************************
Returns the base codes of this pair of nucleotide in a string. The returned
string should be upper case characters if no custom nucleotides were used.

@note Does the same thing that GetPairBaseCodes() does.

@return std::string:
 the base codes of this paire of nucleotide in a string.
*/
// TODO
std::string CHelixIndex::ToString() const
{
	TRACEI("std::string CHelixIndex::ToString() const");
	//std::ostringstream out;
        //out << "[" << m_hiNumber << "," << m_hiType << "]";
	//return out.str();
        std::ostringstream outs;    // Declare an output string stream.
        // TODO
        ////version outs << "[" << m_hiNumber << "," << m_hiType << "]";           // Convert value into a string.
	outs << m_hiNumber << m_hiType;
        return outs.str();     // Get the created string from the output stream.
}


}; // namespace comp
}; // namespace hishape
