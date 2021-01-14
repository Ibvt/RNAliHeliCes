//*****************************************************************************
#include "NucleotidePair.hpp"
#include "DebugTools.hpp"
#include "BioCPP.hpp"

namespace val
{
namespace biocpp
{
/**
Default Constructor
*******************
Creates a new initilized instance of @c CNucleotidePair.

@param const nucleic_acids::SNucleotideCode &tFirstNucleotide:\n
 first nucleotide code value.\n
 Default value: @c nucleic_acids::Unknown

@param const nucleic_acids::SNucleotideCode &tSecondNucleotide:\n
 second nucleotide code value.\n
 Default value: @c nucleic_acids::Unknown

@return CNucleotidePair:
 a new initilized instance of CNucleotidePair.

@see nucleic_acids::SNucleotideCode
@see nucleic_acids::Null
*/
CNucleotidePair::CNucleotidePair(const nucleic_acids::SNucleotideCode &tFirstNucleotide, const nucleic_acids::SNucleotideCode &tSecondNucleotide)
{
	TRACEI("CNucleotidePair::CNucleotidePair(const nucleic_acids::SNucleotideCode &tNucleotide, const nucleic_acids::SNucleotideCode &tNucleotide)");
	m_atNucleotideCodes[0] = tFirstNucleotide;
	m_atNucleotideCodes[1] = tSecondNucleotide;
	m_iIntCode = SNucleotideCodeToIntCode(tFirstNucleotide)*I_NUCLEOTIDES_COUNT + SNucleotideCodeToIntCode(tSecondNucleotide);
}


/**
Copy Constructor
****************
Constructs a new @c CNucleotidePair object that has the same nucleotide codes
as the given @c CNucleotidePair.

@param const CNucleotidePair &tNucleotidePair:
 the nucleotide pair to copy.

@return CNucleotidePair:
 a new instance of CNucleotidePair with the same nucleotide codes as
 @c tNucleotidePair.

*/
CNucleotidePair::CNucleotidePair(const CNucleotidePair &tNucleotidePair)
{
	TRACEI("CNucleotidePair::CNucleotidePair(const CNucleotidePair &tNucleotidePair)");
	m_atNucleotideCodes[0] = tNucleotidePair.m_atNucleotideCodes[0];
	m_atNucleotideCodes[1] = tNucleotidePair.m_atNucleotideCodes[1];
	m_iIntCode = tNucleotidePair.m_iIntCode;
}


/**
Destructor
**********
A standard virtual destructor.
*/
CNucleotidePair::~CNucleotidePair(void)
{
	TRACEI("CNucleotidePair::~CNucleotidePair(void)");
}


/**
@fn std::string CNucleotidePair::GetPairBaseCodes() const
*********************************************************
Returns the base codes of this pair of nucleotide in a string. The returned
string should be upper case characters if no custom nucleotides were used.

@note Does the same thing that ToString() does.

@return std::string:
 the base codes of this pair of nucleotide in a string.
*/
std::string CNucleotidePair::GetPairBaseCodes() const
{
	TRACEI("std::string CNucleotidePair::GetPairBaseCodes() const");
	return std::string((const char*)m_atNucleotideCodes, 2);
}


/**
@fn CNucleotidePair& CNucleotidePair::operator=(const CNucleotidePair &tRightNucleotidePair)
********************************************************************************************
Standard assignation operator. It only copies the nucleotide codes.

@param const CNucleotidePair &tRightNucleotidePair:
 the nucleotides pair that sets the new value.

@return CNucleotidePair&:
 (*this) @c CNucleotidePair with the same nucleotide codes as
 @c tRightNucleotidePair.
*/
CNucleotidePair& CNucleotidePair::operator=(const CNucleotidePair &tRightNucleotidePair)
{
	TRACEI("CNucleotidePair& CNucleotidePair::operator=(const CNucleotidePair &tRightNucleotidePair)");
	m_atNucleotideCodes[0] = tRightNucleotidePair.m_atNucleotideCodes[0];
	m_atNucleotideCodes[1] = tRightNucleotidePair.m_atNucleotideCodes[1];
	m_iIntCode = tRightNucleotidePair.m_iIntCode;
	return *this;
}


/**
@fn bool CNucleotidePair::operator==(const CNucleotidePair &tRightNucleotidePair) const
***************************************************************************************
Tells if two @c CNucleotidePair have the same nucleotide codes in the same
order. (case sensitive)

@param const CNucleotidePair &tRightNucleotidePair:
 the other nucleotide pair for comparison.

@return bool:
 true if pair nucleotides are equal 1 by 1, false otherwise.
*/
bool CNucleotidePair::operator==(const CNucleotidePair &tRightNucleotidePair) const
{
	TRACEI("bool CNucleotidePair::operator==(const CNucleotidePair &tRightNucleotidePair) const");
	return 	(m_iIntCode == tRightNucleotidePair.m_iIntCode);
}


/**
@fn bool CNucleotidePair::operator!=(const CNucleotidePair &tRightNucleotidePair) const
***************************************************************************************
Tells if two @c CNucleotidePair have different nucleotide codes.
(case sensitive)

@param const CNucleotidePair &tRightNucleotidePair:
 the other nucleotide pair for comparison.

@return bool:
 true if at least one of the pair nucleotide is different from its
 corresponding nucleotide in the other pair, false otherwise.
*/
bool CNucleotidePair::operator!=(const CNucleotidePair &tRightNucleotidePair) const
{
	TRACEI("bool CNucleotidePair::operator!=(const CNucleotidePair &tRightNucleotidePair) const");
	return (m_iIntCode != tRightNucleotidePair.m_iIntCode);
}


/**
@fn void CNucleotidePair::SetNucleotide(const ENucleotideIndex &eNucleotideIndex, const nucleic_acids::SNucleotideCode& tNucleotide)
************************************************************************************************************************************
Set the specified nucleotide code of this pair to the new given value.

@param const ENucleotideIndex &eNucleotideIndex:\n
 the nucleotide index. See @c ENucleotideIndex for possible values.

@param const nucleic_acids::SNucleotideCode& tNucleotide:\n
 the new value to use.

@see ENucleotideIndex
*/
void CNucleotidePair::SetNucleotide(const ENucleotideIndex &eNucleotideIndex, const nucleic_acids::SNucleotideCode& tNucleotide)
{
	TRACEI("void CNucleotidePair::SetNucleotide(const ENucleotideIndex &eNucleotideIndex, const nucleic_acids::SNucleotideCode& tNucleotide)");
	m_atNucleotideCodes[eNucleotideIndex] = tNucleotide;
	m_iIntCode = SNucleotideCodeToIntCode(m_atNucleotideCodes[0])*I_NUCLEOTIDES_COUNT + SNucleotideCodeToIntCode(m_atNucleotideCodes[1]);

}


/**
@fn void CNucleotidePair::SetPair(const nucleic_acids::SNucleotideCode &tFirstNucleotide, const nucleic_acids::SNucleotideCode &tSecondNucleotide)
**************************************************************************************************************************************************
Set the nucleotide codes of this pair to the new given values.

@param const nucleic_acids::SNucleotideCode& tFirstNucleotide:\n
 the first new value to use (left nucleotide).

@param const nucleic_acids::SNucleotideCode& tSecondNucleotide:\n
 the second new value to use (right nucleotide).
*/
void CNucleotidePair::SetPair(const nucleic_acids::SNucleotideCode &tFirstNucleotide, const nucleic_acids::SNucleotideCode &tSecondNucleotide)
{
	TRACEI("void CNucleotidePair::SetPair(const nucleic_acids::SNucleotideCode &tFirstNucleotide, const nucleic_acids::SNucleotideCode &tSecondNucleotide)");
	m_atNucleotideCodes[0] = tFirstNucleotide;
	m_atNucleotideCodes[1] = tSecondNucleotide;
	m_iIntCode = SNucleotideCodeToIntCode(m_atNucleotideCodes[0])*I_NUCLEOTIDES_COUNT + SNucleotideCodeToIntCode(m_atNucleotideCodes[1]);
}


/**
@fn const nucleic_acids::SNucleotideCode& CNucleotidePair::operator[](const ENucleotideIndex &eNucleotideIndex) const
*********************************************************************************************************************
Returns the specified nucleotide code of this pair.

@param const ENucleotideIndex &eNucleotideIndex:\n
 the nucleotide index. See @c ENucleotideIndex for possible values.

@return const nucleic_acids::SNucleotideCode&:
 a reference to the wanted nucleotide code (read only).

@see ENucleotideIndex
*/
const nucleic_acids::SNucleotideCode& CNucleotidePair::operator[](const ENucleotideIndex &eNucleotideIndex) const
{
	TRACEI("const nucleic_acids::SNucleotideCode& CNucleotidePair::operator[](const ENucleotideIndex &eNucleotideIndex) const");
	return m_atNucleotideCodes[eNucleotideIndex];
}


/**
@fn int CNucleotidePair::GetIntCode() const
*******************************************
Returns the integer codes of the two nucleotide codes into the same integer.
The high part is the first nucleotide, the low part is the second nucleotide.

@return int:\n
 the pair integer code.
*/
int CNucleotidePair::GetIntCode() const
{
	TRACEI("int CNucleotidePair::GetIntCode() const");
	return m_iIntCode;
}


/**
@fn int CNucleotidePair::GetFirstHalfIntCode() const
****************************************************
Returns the integer code of the first nucleotide code of the pair.

@return int:\n
 the first nucleotide integer code.
*/
int CNucleotidePair::GetFirstHalfIntCode() const
{
	TRACEI("int CNucleotidePair::GetFirstHalfIntCode() const");
	return SNucleotideCodeToIntCode(m_atNucleotideCodes[0]);
}


/**
@fn int CNucleotidePair::GetSecondHalfIntCode() const
*****************************************************
Returns the integer code of the second nucleotide code of the pair.

@return int:\n
 the second nucleotide integer code.
*/
int CNucleotidePair::GetSecondHalfIntCode() const
{
	TRACEI("int CNucleotidePair::GetSecondHalfIntCode() const");
	return SNucleotideCodeToIntCode(m_atNucleotideCodes[1]);
}


/**
@fn val::biocpp::nucleic_acids::SNucleotideCode CNucleotidePair::GetFirstHalfNucleotideCode() const
********************************************************************************
Returns the first nucleotide code of the pair.

@return val::biocpp::nucleic_acids::SNucleotideCode:\n
 the first nucleotide code.
*/
val::biocpp::nucleic_acids::SNucleotideCode CNucleotidePair::GetFirstHalfNucleotideCode() const
{
	TRACEI("val::biocpp::nucleic_acids::SNucleotideCode CNucleotidePair::GetFirstHalfNucleotideCode() const");
	return m_atNucleotideCodes[0];
}


/**
@fn val::biocpp::nucleic_acids::SNucleotideCode CNucleotidePair::GetSecondHalfNucleotideCode() const
*************************************************************************************
Returns the second nucleotide code of the pair.

@return val::biocpp::nucleic_acids::SNucleotideCode:\n
 the second nucleotide code.
*/
val::biocpp::nucleic_acids::SNucleotideCode CNucleotidePair::GetSecondHalfNucleotideCode() const
{
	TRACEI("val::biocpp::nucleic_acids::SNucleotideCode CNucleotidePair::GetSecondHalfNucleotideCode() const");
	return m_atNucleotideCodes[1];
}


/**
@fn std::string CNucleotidePair::ToString() const
*************************************************
Returns the base codes of this pair of nucleotide in a string. The returned
string should be upper case characters if no custom nucleotides were used.

@note Does the same thing that GetPairBaseCodes() does.

@return std::string:
 the base codes of this paire of nucleotide in a string.
*/
std::string CNucleotidePair::ToString() const
{
	TRACEI("std::string CNucleotidePair::ToString() const");
	return std::string((const char*)m_atNucleotideCodes, 2);
}


/**
@fn CNucleotidePair CNucleotidePair::FusionNucleotides(const CNucleotidePair &tFirstHalf, const CNucleotidePair &tSecondHalf)
*****************************************************************************************************************************
Returns a pair of nucleotide made of the first base of the first given pair
and the second base of the second given pair.
For example, it the first pair is "GC" and the second pair is "AU", returned
pair will be "GU".
@note: This is a static function.

@param const CNucleotidePair &tFirstHalf:
 the pair containing the first nucleotide that will be used to create the new
 pair.

@param const CNucleotidePair &tSecondHalf:
 the pair containing the second nucleotide that will be used to create the
 new pair.

@return CNucleotidePair:
 the new pair mode from 1 base of each provided pair.

*/
CNucleotidePair CNucleotidePair::FusionNucleotides(const CNucleotidePair &tFirstHalf, const CNucleotidePair &tSecondHalf)
{
	CNucleotidePair tNewPair;
	tNewPair.SetPair(tFirstHalf.GetFirstHalfNucleotideCode(), tSecondHalf.GetSecondHalfNucleotideCode());
	return tNewPair;
}

}; // namespace biocpp
}; // namespace val
