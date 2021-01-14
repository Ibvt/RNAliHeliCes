//*****************************************************************************
#include "NAPrimaryStructure.hpp"
#include <stdexcept>
#include <iterator>
#include <list>
#include <deque>
#include <sstream>
#include "DebugTools.hpp"


namespace val
{
namespace biocpp
{
/**
Default Constructor
*******************
Constructs a new @c CNAPrimaryStructure object that contains an empty sequence
of nucleotides.

@return CNAPrimaryStructure:\n
 a new instance of @c CNAPrimaryStructure.
*/
CNAPrimaryStructure::CNAPrimaryStructure(void)
: m_strDescription("")
, m_strName("Unnamed")
{
	TRACEI("CNAPrimaryStructure::CNAPrimaryStructure(void)");
}


/**
Copy Constructor
****************
Constructs a new @c CNAPrimaryStructure object that contains a copy of the
given sequence of nucleotides.

@return CNAPrimaryStructure:\n
 a new instance of CNAPrimaryStructure wich is a copy of the given sequence.
*/
CNAPrimaryStructure::CNAPrimaryStructure(const CNAPrimaryStructure &tPrimaryStructure)
: m_strDescription(tPrimaryStructure.m_strDescription)
, m_strName(tPrimaryStructure.m_strName)
, m_tSequence(tPrimaryStructure.m_tSequence)
{
	TRACEI("CNAPrimaryStructure::CNAPrimaryStructure(const CNAPrimaryStructure &tPrimaryStructure)");
}


/**
@fn CNAPrimaryStructure::CNAPrimaryStructure(const std::string &strSequence, int iFromOffset, int iLength)
**********************************************************************************************************
Use the given string @c strSequence as the sequence. Only @c iLength base
codes are read from @c iFromOffset (0-based) and added to the primary
structure. If the sequence length is omitted or negative, all the base codes
from @c iFromOffset to the end of the string are added. If the offset is
negative or larger than the string length, the sequence will be empty

@note if the string contains invalid base codes, they will be added as
 @c nucleic_acids::Unknown.

@param const std::string &strSequence:\n
 a string of base codes in the 5' to 3' order (case unsensitive).

@param int iFromOffset:\n
 base codes of @c strSequence will be read starting from this offset.\n
 Default value: 0

@param int iLength:\n
 count of base codes to use for the sequence. Zero means "empty sequence".
 A negative value means the whole string should be used (size unknown).\n
 Default value: -1

@return CNAPrimaryStructure:\n
 a new instance of CNAPrimaryStructure.
*/
CNAPrimaryStructure::CNAPrimaryStructure(const std::string &strSequence, int iFromOffset, int iLength)
{
	TRACEI("CNAPrimaryStructure::CNAPrimaryStructure(const std::string &strSequence, int iFromOffset, int iLength)");
	int iEndLimit = iFromOffset;
	if (0 < iLength)
	{
		// sequence length specified
		iEndLimit = iFromOffset + iLength;
		if (iEndLimit > (signed)strSequence.length())
		{
			iEndLimit = (signed)strSequence.length();
		}
	}
	else if (0 > iLength)
	{
		// sequence length unknown
		iEndLimit = (signed)strSequence.length();
	}

	// check offset (do nothing if invalid)
	if (iFromOffset >= 0)
	{
		for (int i = iFromOffset; i < iEndLimit; i++)
		{
			m_tSequence.push_back(CharToSNucleotideCodeSafe(strSequence[i]));
		}
	}
}


/**
Destructor
**********
Virtual destructor.
*/
CNAPrimaryStructure::~CNAPrimaryStructure(void)
{
	TRACEI("CNAPrimaryStructure::~CNAPrimaryStructure(void)");
	ClearSequence();
}


/**
@fn void CNAPrimaryStructure::ClearSequence(void)
*************************************************
Clears sequence (remove all nucleotides).
*/
void CNAPrimaryStructure::ClearSequence(void)
{
	TRACEI("void CNAPrimaryStructure::ClearSequence(void)");
	m_tSequence.clear();
}


/**
@fn CNAPrimaryStructure::CStructureIterator CNAPrimaryStructure::GetIteratorAt(int iAtPosition)
***********************************************************************************************
Returns an iterator on the sequence at the given postion.

@param int iAtPosition:\n
 requested position of the iterator (1-based).

@return std::deque<nucleic_acids::SNucleotideCode>::iterator:\n
 an iterator on the sequence at the given postion.

@throw std::out_of_range:\n
 if the position is outside the sequence.
*/
CNAPrimaryStructure::CStructureIterator CNAPrimaryStructure::GetIteratorAt(int iAtPosition)
{
	TRACEI("CNAPrimaryStructure::CStructureIterator CNAPrimaryStructure::GetIteratorAt(int iAtPosition)");
	int iSequenceSize = (signed) m_tSequence.size();
	// keep in mind that iAtPosition is not a 0-based index
	// ie. first position is 1
	if ((iAtPosition <= 0) || (iSequenceSize < iAtPosition))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNAPrimaryStructure::GetIteratorAt(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	return (m_tSequence.begin() + (iAtPosition-1));
}


/**
@fn const CNAPrimaryStructure::CStructureConstIterator CNAPrimaryStructure::GetConstIteratorAt(int iAtPosition) const
****************************************************************************************************************************
Returns a const_iterator on the sequence at the given postion.

@param int iAtPosition:\n
 requested position of the iterator.

@return const std::deque<nucleic_acids::SNucleotideCode>::const_iterator:\n
 a const_iterator on the sequence at the given postion (1-based).

@throw std::out_of_range:\n
 if the position is outside the sequence.
*/
const CNAPrimaryStructure::CStructureConstIterator CNAPrimaryStructure::GetConstIteratorAt(int iAtPosition) const
{
	TRACEI("const CNAPrimaryStructure::CStructureConstIterator CNAPrimaryStructure::GetConstIteratorAt(int iAtPosition) const");
	int iSequenceSize = (signed) m_tSequence.size();
	// keep in mind that iPosition is not a 0-based index, first position is 1
	if ((iAtPosition <= 0) || (iSequenceSize < iAtPosition))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNAPrimaryStructure::GetConstIteratorAt(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	return (m_tSequence.begin() + (iAtPosition-1));
}


/**
@fn void CNAPrimaryStructure::AddOn3Prime(const nucleic_acids::SNucleotideCode &tNucleotide)
********************************************************************************************
Appends a nucleotide on the 3' end of the sequence (end of the sequence).
This method does the same as @c Append.

@param const nucleic_acids::SNucleotideCode &tNucleotide:\n
 the nucleotide to append.

@see nucleic_acids::SNucleotideCode
@see Append
*/
void CNAPrimaryStructure::AddOn3Prime(const nucleic_acids::SNucleotideCode &tNucleotide)
{
	TRACEI("void CNAPrimaryStructure::AddOn3Prime(const nucleic_acids::SNucleotideCode &tNucleotide)");
	m_tSequence.push_back(tNucleotide);
}


/**
@fn void CNAPrimaryStructure::AddOn5Prime(const nucleic_acids::SNucleotideCode &tNucleotide)
*************************************************************************
Adds a nucleotide on the 5' end of the sequence (start of the sequence).

@param const nucleic_acids::SNucleotideCode &tNucleotide:\n
 the nucleotide to add.

@see nucleic_acids::SNucleotideCode
*/
void CNAPrimaryStructure::AddOn5Prime(const nucleic_acids::SNucleotideCode &tNucleotide)
{
	TRACEI("void CNAPrimaryStructure::AddOn5Prime(const nucleic_acids::SNucleotideCode &tNucleotide)");
	m_tSequence.push_front(tNucleotide);
}


/**
@fn void CNAPrimaryStructure::Append(const nucleic_acids::SNucleotideCode &tNucleotide)
********************************************************************
Appends a nucleotide on the 3' end of the sequence (end of the sequence).
This method does the same as @c AddOn3Prime.

@param const nucleic_acids::SNucleotideCode &tNucleotide:\n
 the nucleotide code to append.

@see nucleic_acids::SNucleotideCode
*/
void CNAPrimaryStructure::Append(const nucleic_acids::SNucleotideCode &tNucleotide)
{
	TRACEI("void CNAPrimaryStructure::Append(const nucleic_acids::SNucleotideCode &tNucleotide)");
	m_tSequence.push_back(tNucleotide);
}


/**
@fn int CNAPrimaryStructure::AppendSequence(const std::string &strSequence, int iFromOffset, int iAppendLength)
***************************************************************************************************************
Appends the given sequence @c strSequence to this one on the 3' end. Only
@c iAppendLength base codes are read from @c iFromOffset (0-based) and added
to the primary structure. If the sequence length is omitted or negative, all
the base codes from @c iFromOffset to the end of the string are added.

@note if the string contains invalid base codes, they will be added as
 @c nucleic_acids::Unknown.

@param const std::string &strSequence:\n
 a string of base codes in the 5' to 3' order (case unsensitive).

@param int iFromOffset:\n
 base codes of @c strSequence will be read starting from this offset.\n
 Default value: 0

@param int iAppendLength:\n
 count of base codes to append. Zero means "append nothing".
 A negative value means the whole sequence should be append (size unknown).\n
 Default value: -1

@return int:
 count of nucleotides appended.
*/
int CNAPrimaryStructure::AppendSequence(const std::string &strSequence, int iFromOffset, int iAppendLength)
{
	TRACEI("int CNAPrimaryStructure::AppendSequence(const std::string &strSequence, int iFromOffset, int iAppendLength)");
	int iAppended = iFromOffset;
	if (0 <= iFromOffset)
	{
		int iEndLimit = 0;
		if (0 < iAppendLength)
		{ // sequence length specified
			iEndLimit = iFromOffset + iAppendLength;
			if (iEndLimit > (signed)strSequence.length())
			{
				iEndLimit = (signed)strSequence.length();
			}
		}
		else if (0 > iAppendLength)
		{ // sequence length unknown
			iEndLimit = (signed)strSequence.length();
		}
		while (iAppended < iEndLimit)
		{
			m_tSequence.push_back(CharToSNucleotideCodeSafe(strSequence[iAppended++]));
		}
	}
	return iAppended - iFromOffset;
}


/**
@fn int CNAPrimaryStructure::AppendSequence(const CNAPrimaryStructure &tSequence, int iSequenceOffset, int iAppendLength)
****************************************************************************************************************************************
Appends the given sequence @c tSequence to this one on the 3' end.

@param const CNAPrimaryStructure &tSequence:\n
 a sequence of nucleotides to append. 

@param int iFromOffset:\n
 offset from where @c tSequence should be started.\n
 Default value: 0

@param int iAppendLength:\n
 length of the copy to append. Zero means "append nothing". A negative value
 means the whole sequence should be append (unknown size).\n
 Default value: -1

@return int:
 count of nucleotides appended.

@throw std::out_of_range:\n
 if the given offset and length leads out of @c tSequence.
*/
int CNAPrimaryStructure::AppendSequence(const CNAPrimaryStructure &tSequence, int iFromOffset, int iAppendLength)
{
	TRACEI("int CNAPrimaryStructure::AppendSequence(const CNAPrimaryStructure &tSequence)");
	int iBeginPosition = iFromOffset;
	int iEndLimit = 0;
	if ((0 <= iBeginPosition) && (tSequence.GetLength() > iBeginPosition))
	{
		if (0 < iAppendLength)
		{ // sequence length specified
			iEndLimit = iFromOffset + iAppendLength;
			if (iEndLimit > tSequence.GetLength())
			{
				iEndLimit = tSequence.GetLength();
			}
		}
		else if (0 > iAppendLength)
		{ // sequence length unknown
			iEndLimit = tSequence.GetLength();
		}
		else
		{
			return  0;
		}
		CStructureConstIterator tBeginPosition = tSequence.GetConstIteratorAt(iBeginPosition+1);
		CStructureConstIterator tEndPosition = tSequence.GetConstIteratorAt(iEndLimit) + 1;
		m_tSequence.insert(m_tSequence.end(), tBeginPosition, tEndPosition);
		return iEndLimit - iFromOffset;
	}
	else
	{
		return 0;
	}
}


/**
@fn void CNAPrimaryStructure::InsertAt(const nucleic_acids::SNucleotideCode &tNucleotide, int iAtPosition, int iCopyCount)
**************************************************************************************************************************
Inserts @c iCopyCount copies of the nucleotide @c tBaseCode at the position
@c iPosition in this sequence (before the nucleotide that has the given
position).

@note The first position in the sequence is @b 1 and
 @code InsertAt(nucleic_acids::Adenine, 1) @endcode will do the same as
 @code AddOn5Prime(nucleic_acids::Adenine) @endcode

@param const nucleic_acids::SNucleotideCode &tNucleotide:\n
 The nucleotide code to insert.

@param int iAtPosition:\n
 position where the nucleotide should be inserted (1-based).

@param int iCopyCount:\n
 number of copies of the nucleotide to insert.\n
 Default value: 1.

@throw std::out_of_range:\n
 if the position is not in the sequence. It often appends when trying to use a
 position of 0 while the first position in the sequence is 1.

@see nucleic_acids::SNucleotideCode
*/
void CNAPrimaryStructure::InsertAt(const nucleic_acids::SNucleotideCode &tNucleotide, int iAtPosition, int iCopyCount)
{
	TRACEI("void CNAPrimaryStructure::InsertAt(const nucleic_acids::SNucleotideCode &tNucleotide, int iAtPosition, int iCopyCount)");

	CStructureIterator tInsertPosition = GetIteratorAt(iAtPosition);
	// creates a temporary deque to insert all the elements at once
	std::deque<nucleic_acids::SNucleotideCode> tTempSequence;
	for (int i = 0; i < iCopyCount; i++)
	{
		tTempSequence.push_back(tNucleotide);
	}
	// inserts the new elements
	m_tSequence.insert(tInsertPosition, tTempSequence.begin(), tTempSequence.end());
}


/**
@fn void CNAPrimaryStructure::InsertSequenceAt(const std::string &strSequence, int iAtPosition, int iSequenceOffset, int iInsertLength, int iCopyCount)
*******************************************************************************************************************************************************
Inserts @c iCopyCount copies of the given nucleotides sequence (with the given
offset and length) at the position @c iAtPosition in this sequence.

@note The first position in the sequence is 1 (1-based).

@param const std::string &strSequence:\n
 a sequence of @c iSequenceLength nucleotides.

@param  int iAtPosition:\n
 the position in current sequence where the new sequence should be inserted.

@param int iFromOffset:\n
 offset (0-based) from where @c strSequence should be read.
 Default value: 0.

@param int iInsertLength:\n
 the length of one copy of the sequence to insert. Zero means "append nothing".
 A negative value means the whole sequence should be inserted (unknown size).
 Default value: -1.

@param int iCopyCount:\n
 number of copies of the new sequence that has to be inserted.\n
 Default value: 1.

@return int:
 count of nucleotides inserted (copies included).
*/
int CNAPrimaryStructure::InsertSequenceAt(const std::string &strSequence, int iAtPosition, int iFromOffset, int iInsertLength, int iCopyCount)
{
	TRACEI("void CNAPrimaryStructure::InsertSequenceAt(const std::string &strSequence, int iAtPosition, int iFromOffset, int iInsertLength, int iCopyCount)");
	int iRealInsertLength = iInsertLength;
	if (0 > iRealInsertLength)
	{
		iRealInsertLength = (signed)strSequence.length() - iFromOffset;
	}
	// compute insertion size
	int iInsertCount = iCopyCount*iRealInsertLength;
	if (0 >= iInsertCount)
	{
		iInsertCount = 0;
	}
	else
	{
		// creates the temporary sequence to insert
		std::deque<nucleic_acids::SNucleotideCode> tTempSequence;
		// gets the position
		CStructureIterator tInsertPosition = GetIteratorAt(iAtPosition);
		// let iInsertLength be the end insertion limit
		iRealInsertLength += iFromOffset;
		// insert iCopyCount copies
		for (int i = 0; i < iCopyCount; i++)
		{
			for (int j = iFromOffset; j < iRealInsertLength; j++)
			{
				tTempSequence.push_back(CharToSNucleotideCodeSafe(strSequence[j]));
			}
			m_tSequence.insert(tInsertPosition, tTempSequence.begin(), tTempSequence.end());
			tTempSequence.clear();
		}
	}
	return iInsertCount;
}


/**
@fn void CNAPrimaryStructure::InsertSequenceAt(const CNAPrimaryStructure &tSequence, int iFromOffset, int iInsertLength, int iAtPosition, int iCopyCount)
*****************************************************************************************************************************************************************************************
Inserts @c iCopyCount copies of the given nucleotides sequence at the position
@c iAtPosition in this sequence.

@note The first position in the sequence is 1.

@param const CNAPrimaryStructure &tSequence:\n
 a sequence of nucleotides.

@param int iFromOffset:\n
 offset (0-based) from where @c tSequence should be read.

@param int iInsertLength:\n
 the length of the sequence to insert.

@param int iAtPosition:\n
 the position in current sequence where the new sequence should be inserted.
 (1-based)

@param int iCopyCount:\n
 number of copies of the new sequence that has to be inserted.\n
 Default value: 1.
*/
void CNAPrimaryStructure::InsertSequenceAt(const CNAPrimaryStructure &tSequence, int iFromOffset, int iInsertLength, int iAtPosition, int iCopyCount)
{
	TRACEI("void CNAPrimaryStructure::InsertSequenceAt(const CNAPrimaryStructure &tSequence, int iSequenceOffset, int iInsertLength, int iAtPosition, int iCopyCount)");
	int iRealInsertLength = iInsertLength;
	if (0 > iRealInsertLength)
	{ // real size unknown, compute it
		iRealInsertLength = (signed)tSequence.m_tSequence.size();
	}
	if (0 < iRealInsertLength)
	{
		// gets the position
		CStructureIterator tInsertPosition = GetIteratorAt(iAtPosition);
		// creates the sequence to insert
		std::deque<nucleic_acids::SNucleotideCode> tTempSequence;
		// insert iCopyCount copies
		for (int i = 0; i < iCopyCount; i++)
		{
			for (int j = iFromOffset; j < iFromOffset + iRealInsertLength; j++)
			{
				tTempSequence.push_back(tSequence.m_tSequence[j]);
			}
			m_tSequence.insert(tInsertPosition, tTempSequence.begin(), tTempSequence.end());
			tTempSequence.clear();
		}
	}
}


/**
@fn void CNAPrimaryStructure::RemoveOn3Prime(int iNucleotideCount)
*************************************************************************
Removes @c iNucleotideCount nucleotides starting from the 3' end (removes the
lasts @c iNucleotideCount nucleotides from the end of the sequence).

@param int iNucleotideCount:\n
 number of nucleotides to remove.\n
 Default value: 1.
*/
void CNAPrimaryStructure::RemoveOn3Prime(int iNucleotideCount)
{
	TRACEI("void CNAPrimaryStructure::RemoveOn3Prime(int iNucleotideCount)");
	for (int i=0; i < iNucleotideCount; i++)
	{
		m_tSequence.pop_back();
	}
}


/**
@fn void CNAPrimaryStructure::RemoveOn5Prime(int iNucleotideCount)
*************************************************************************
Removes @c iNucleotideCount nucleotides starting from the 5' end (removes the
lasts @c iNucleotideCount nucleotides from the end of the sequence).

@param int iNucleotideCount:\n
 number of nucleotides to remove.\n
 Default value: 1.

@throw std::invalid_argument:\n
 if the given length is invalid (ie. negative).
*/
void CNAPrimaryStructure::RemoveOn5Prime(int iNucleotideCount)
{
	TRACEI("void CNAPrimaryStructure::RemoveOn5Prime(int iNucleotideCount)");
	if (0 > iNucleotideCount)
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNAPrimaryStructure::RemoveOn5Prime(" << iNucleotideCount << "): " << STR_ERROR_INVALID_LENGTH << std::endl;
		throw std::invalid_argument(stmErrorData.str());
	}
	for (int i=0; i <iNucleotideCount; i++)
	{
		m_tSequence.pop_front();
	}
}


/**
@fn void CNAPrimaryStructure::RemoveAt(int iAtPosition, int iNucleotideCount)
*****************************************************************************
Removes @c iNucleotideCount nucleotides starting from the position
@c iAtPosition. Nucleotides are removed in the 5' to 3' way.

@note The first position in the sequence is 1.

@param int iAtPosition:\n
 position of the first nucleotide to be removed.

@param int iNucleotideCount:\n
 number of nucleotides to remove.\n
 Default value: 1.

@throw std::invalid_argument:\n
 if the given length is invalid (ie. negative).
*/
void CNAPrimaryStructure::RemoveAt(int iAtPosition, int iNucleotideCount)
{
	TRACEI("void CNAPrimaryStructure::RemoveAt(int iAtPosition, int iNucleotideCount)");
	if (0 > iNucleotideCount)
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNAPrimaryStructure::RemoveAt(" << iAtPosition << ", " << iNucleotideCount << "): " << STR_ERROR_INVALID_LENGTH << std::endl;
		throw std::invalid_argument(stmErrorData.str());
	}
	// removes the pointers from the deque
	// 1-based indexing
	CStructureIterator tIterStart = GetIteratorAt(iAtPosition);
	CStructureIterator tIterEnd = tIterStart + iNucleotideCount;
	m_tSequence.erase(tIterStart, tIterEnd);
}


/**
@fn void CNAPrimaryStructure::RemoveBefore(int iAtPosition, int iNucleotideCount)
***********************************************************************************************
Removes @c iNucleotideCount nucleotides before the position @c iAtPosition.

@note The first position in the sequence is 1.

@note The position of the nucleotide at @c iPosition will be changed after the
 removal!

@param int iAtPosition:\n
 position of the first nucleotide after the nucleotides to be removed in current
 sequence (position before removal).

@param int iNucleotideCount:\n
 number of nucleotides to remove.\n
 Default value: 1.

@throw std::invalid_argument:\n
 if the given length is invalid (ie. negative).
*/
void CNAPrimaryStructure::RemoveBefore(int iAtPosition, int iNucleotideCount)
{
	TRACEI("void CNAPrimaryStructure::RemoveBefore(int iAtPosition, int iNucleotideCount)");
	if (0 > iNucleotideCount)
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNAPrimaryStructure::RemoveBefore(" << iAtPosition << ", " << iNucleotideCount << "): " << STR_ERROR_INVALID_LENGTH << std::endl;
		throw std::invalid_argument(stmErrorData.str());
	}
	// removes the pointers from the deque
	// 1-based indexing
	CStructureIterator tIterStart = GetIteratorAt(iAtPosition - iNucleotideCount);
	CStructureIterator tIterEnd = tIterStart + iNucleotideCount;
	m_tSequence.erase(tIterStart, tIterEnd);
}


/**
@fn void CNAPrimaryStructure::RemoveAfter(int iPosition, int iNucleotideCount)
******************************************************************************
Removes @c iNucleotideCount nucleotides after the position @c iAtPosition.

@note The first position in the sequence is 1.

@param int iAtPosition:\n
 position of the first nucleotide before the nucleotides to be removed in
 current sequence (position before removal).

@param int iNucleotideCount:\n
 number of nucleotides to remove.\n
 Default value: 1.

@throw std::invalid_argument:\n
 if the given length is invalid (ie. negative).
*/
void CNAPrimaryStructure::RemoveAfter(int iAtPosition, int iNucleotideCount)
{
	TRACEI("void CNAPrimaryStructure::RemoveAfter(int iAtPosition, int iNucleotideCount)");
	if (0 > iNucleotideCount)
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNAPrimaryStructure::RemoveAfter(" << iAtPosition << ", " << iNucleotideCount << "): " << STR_ERROR_INVALID_LENGTH << std::endl;
		throw std::invalid_argument(stmErrorData.str());
	}
	// removes the pointers from the deque
	// 1-based indexing
	CStructureIterator tIterStart = GetIteratorAt(iAtPosition+1);
	CStructureIterator tIterEnd = tIterStart + iNucleotideCount;
	m_tSequence.erase(tIterStart, tIterEnd);
}


/**
@fn nucleic_acids::SNucleotideCode& CNAPrimaryStructure::GetAt(int iAtPosition) const
*************************************************************************
Returns the nucleotide code of the nucleotide at the given position.

@param int iAtPosition:\n
 Position of the nucleotide which base code is needed.
 
@return nucleic_acids::SNucleotideCode&:\n
 the code of the nucleotide at this position

@throw std::out_of_range:\n
 if the given position is outside the sequence.

@see nucleic_acids::SNucleotideCode
*/
const nucleic_acids::SNucleotideCode& CNAPrimaryStructure::GetAt(int iAtPosition) const
{
	TRACEI("const nucleic_acids::SNucleotideCode& CNAPrimaryStructure::GetAt(int iAtPosition) const");
	std::deque<nucleic_acids::SNucleotideCode>::size_type tSequenceSize = m_tSequence.size();
	// keep in mind that iPosition is not a 0-based index, first position is 1
	if ((iAtPosition <= 0) || ((signed)tSequenceSize < iAtPosition))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNAPrimaryStructure::GetAt(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	return m_tSequence[iAtPosition-1];
}


/**
@fn void CNAPrimaryStructure::SetAt(int iAtPosition, const nucleic_acids::SNucleotideCode& tNewNucleotide)
**********************************************************************************************
Changes the nucleotide code at the given position to the new value.

@param int iAtPosition:\n
 Position of the nucleotide to replace (1-based).
 
@param const nucleic_acids::SNucleotideCode& tNewNucleotide:\n
 the new nucleotide to set at that position.

@throw std::out_of_range:\n
 if the given position is outside the sequence.

@see nucleic_acids::SNucleotideCode
*/
void CNAPrimaryStructure::SetAt(int iAtPosition, const nucleic_acids::SNucleotideCode& tNewNucleotide)
{
	TRACEI("void CNAPrimaryStructure::SetAt(int iAtPosition, const nucleic_acids::SNucleotideCode& tNewNucleotide)");
	std::deque<nucleic_acids::SNucleotideCode>::size_type tSequenceSize = m_tSequence.size();
	// keep in mind that iPosition is not a 0-based index, first position is 1
	if ((iAtPosition <= 0) || ((signed)tSequenceSize < iAtPosition))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNAPrimaryStructure::SetAt(" << iAtPosition << ", SNucleotideCode): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	m_tSequence[iAtPosition-1] = tNewNucleotide;
}


/**
@fn int CNAPrimaryStructure::GetLength(void) const
**************************************************
Returns the length of the sequence (number of nucleotides).

@return int:\n
 count of nucleotides in the sequence (including gaps and other special codes).
*/
int CNAPrimaryStructure::GetLength(void) const
{
	TRACEI("int CNAPrimaryStructure::GetLength(void) const");
	return (signed)m_tSequence.size();
}


/**
@fn void CNAPrimaryStructure::TranscriptToDNA(void)
***************************************************
Replaces all the Uridines in current sequence by Thymines (replaces U by T).
*/
void CNAPrimaryStructure::TranscriptToDNA(void)
{
	TRACEI("void CNAPrimaryStructure::TranscriptToDNA(void)");
	for (CStructureIterator tIter = m_tSequence.begin(); tIter != m_tSequence.end(); tIter++)
	{
		if (nucleic_acids::Uridine == (*tIter))
		{
			*tIter = nucleic_acids::Thymine;
		}
	}
}


/**
@fn void CNAPrimaryStructure::TranscriptToRNA(void)
***************************************************
Replaces all the Thymines in current sequence by Uridines (replaces T by U).
*/
void CNAPrimaryStructure::TranscriptToRNA(void)
{
	TRACEI("void CNAPrimaryStructure::TranscriptToRNA(void)");
	for (CStructureIterator tIter = m_tSequence.begin(); tIter != m_tSequence.end(); tIter++)
	{
		if (nucleic_acids::Thymine == (*tIter))
		{
			*tIter = nucleic_acids::Uridine;
		}
	}
}


/**
@fn CNAPrimaryStructure CNAPrimaryStructure::RevertSequence(void) const
***********************************************************************
Returns a new sequence wich is the exact copy of this sequence except it is in
inverted order: the first nucleotide of current sequence is the last one on the
new sequence and vice-versa. In biological terms, the resulting sequence has
respectively 5' and 3' ends where this sequence has 3' and 5' ends.

@return CNAPrimaryStructure:\n
 a new sequence wich is in inverted order.
*/
CNAPrimaryStructure CNAPrimaryStructure::RevertSequence(void) const
{
	TRACEI("CNAPrimaryStructure CNAPrimaryStructure::RevertSequence(void) const");
	CNAPrimaryStructure tRevertedSequence;
	for (std::deque<nucleic_acids::SNucleotideCode>::const_reverse_iterator tIter = m_tSequence.rend(); tIter < m_tSequence.rbegin(); tIter++)
	{
		tRevertedSequence.AddOn5Prime(*tIter);
	}
	return tRevertedSequence;
}


/**
@fn CNAPrimaryStructure CNAPrimaryStructure::ComplementSequence(void) const
***************************************************************************
Returns a new sequence wich complements current sequence. The new sequence is a
base-per-base complement of the current sequence and can be seen, in the case
of DNA, as the single-stranded complementary DNA. Does the same as @c operator~.

@note in case of base codes of groups of bases, the group is complemented by
 a group that contains all the complement bases. For example, V is considered as
 the complement of B. See @c GetComplementnucleic_acids::SNucleotideCode for details.
 
@return CNAPrimaryStructure:\n
 a new sequence wich complements current sequence.

@see GetComplementnucleic_acids::SNucleotideCode
@see operator~(void)
*/
CNAPrimaryStructure CNAPrimaryStructure::ComplementSequence(void) const
{
	TRACEI("CNAPrimaryStructure CNAPrimaryStructure::ComplementSequence(void) const");
	CNAPrimaryStructure tComplementSequence;
	for (CStructureConstIterator tIter = m_tSequence.begin(); tIter != m_tSequence.end(); tIter++)
	{
		tComplementSequence.AddOn3Prime(GetComplementSNucleotideCode(*tIter));
	}
	return tComplementSequence;
}


/**
@fn std::string CNAPrimaryStructure::ToString(void) const
*********************************************************
Returns a string that contains the full sequence as characters (in capital
letters).

@return std::string:\n
 a new independant string that contains the sequence.
*/
std::string CNAPrimaryStructure::ToString(void) const
{
	TRACEI("std::string CNAPrimaryStructure::ToString(void) const");
	std::string strSeq;
	for (CStructureConstIterator tIter = m_tSequence.begin(); tIter != m_tSequence.end(); tIter++)
	{
		strSeq += SNucleotideCodeToChar(*tIter);
	}
	return strSeq;
}


/**
@fn CNAPrimaryStructure& CNAPrimaryStructure::operator=(const CNAPrimaryStructure &tSequence)
*********************************************************************************************
Assignation operator. Copies the given sequence into this one.

@param const CNAPrimaryStructure &tSequence:\n
 the sequence to assign to the current one.

@return CNAPrimaryStructure&:\n
 returns (*this) object.
*/
CNAPrimaryStructure& CNAPrimaryStructure::operator=(const CNAPrimaryStructure &tSequence)
{
	TRACEI("CNAPrimaryStructure& CNAPrimaryStructure::operator=(const CNAPrimaryStructure &tSequence)");
	// check if it's a self-assignation
	if (this == &tSequence)
	{
		return *this;
	}
	ClearSequence();
	m_tSequence = tSequence.m_tSequence;
	m_strDescription = tSequence.m_strDescription;
	m_strName = tSequence.m_strName;

	return *this;
}


/**
@fn const nucleic_acids::SNucleotideCode& CNAPrimaryStructure::operator[](int iAtPosition) const
************************************************************************************
Returns the nucleotide code at the given position.

@param int iAtPosition:\n
 position (1-based) of the nucleotide to know.

@return const nucleic_acids::SNucleotideCode&:\n
 the nucleotide code at that position.

@throw std::out_of_range:\n
 if the given position is outside the sequence.

@see nucleic_acids::SNucleotideCode
*/
const nucleic_acids::SNucleotideCode& CNAPrimaryStructure::operator[](int iAtPosition) const
{
	TRACEI("const nucleic_acids::SNucleotideCode& CNAPrimaryStructure::GetAt(int iAtPosition) const");
	std::deque<nucleic_acids::SNucleotideCode>::size_type tSequenceSize = m_tSequence.size();
	// keep in mind that iPosition is not a 0-based index, first position is 1
	if ((iAtPosition <= 0) || ((signed)tSequenceSize < iAtPosition))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNAPrimaryStructure::operator[](" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	return m_tSequence[iAtPosition-1];
}


/**
@fn CNAPrimaryStructure CNAPrimaryStructure::operator+(const CNAPrimaryStructure &tRightSequence) const
*******************************************************************************************************
Returns the concatenation of the two sequences in a new sequence.

@param const CNAPrimaryStructure &tRightSequence:\n
 the second sequence to append after the current one.

@return CNAPrimaryStructure:\n
 a new concatenation sequence of the two sequences.
*/
CNAPrimaryStructure CNAPrimaryStructure::operator+(const CNAPrimaryStructure &tRightSequence) const
{
	TRACEI("CNAPrimaryStructure CNAPrimaryStructure::operator+(const CNAPrimaryStructure &tRightSequence) const");
	CNAPrimaryStructure tNewSequence(*this);
	tNewSequence.AppendSequence(tRightSequence);
	return tNewSequence;
}


/**
@fn CNAPrimaryStructure& CNAPrimaryStructure::operator+=(const CNAPrimaryStructure &tRightSequence)
***************************************************************************************************
Concatenate the given sequence to the current one.

@param const CNAPrimaryStructure &tRightSequence:\n
 the second sequence to append after the current one.

@return CNAPrimaryStructure&:\n
 (*this) current sequence with @c tRightSequence sequence concatenated after.
*/
CNAPrimaryStructure& CNAPrimaryStructure::operator+=(const CNAPrimaryStructure &tRightSequence)
{
	TRACEI("CNAPrimaryStructure& CNAPrimaryStructure::operator+=(const CNAPrimaryStructure &tRightSequence)");
	AppendSequence(tRightSequence);
	return *this;
}


/**
@fn CNAPrimaryStructure CNAPrimaryStructure::operator~(void) const
******************************************************************
Returns the complementary sequence of current one just like
@c ComplementSequence(void). Nucleotides are complemented using
@c GetComplementSNucleotideCode function.

@note in case of base codes of groups of bases, the group is complemented by
 a group that contains all the complement bases. For example, V is considered as
 the complement of B. See @c GetComplementSNucleotideCode for details.

@return CNAPrimaryStructure:\n
 the complementary sequence.
 
@see GetComplementSNucleotideCode
@see ComplementSequence(void)
*/
CNAPrimaryStructure CNAPrimaryStructure::operator~(void) const
{
	TRACEI("CNAPrimaryStructure CNAPrimaryStructure::operator~(void) const");
	CNAPrimaryStructure tComplementSequence;
	for (CStructureConstIterator tIter = m_tSequence.begin(); tIter != m_tSequence.end(); tIter++)
	{
		tComplementSequence.m_tSequence.push_front(GetComplementSNucleotideCode(*tIter));
	}
	return tComplementSequence;
}


/**
@fn bool CNAPrimaryStructure::operator==(const CNAPrimaryStructure &tSequence) const
************************************************************************************
Tells if two sequences are identical.

@param const CNAPrimaryStructure &tSequence:\n
 a sequence to compare.

@return bool:\n
 true if the sequences have the same nucleotides at the same place, false
 otherwise.
*/
bool CNAPrimaryStructure::operator==(const CNAPrimaryStructure &tSequence) const
{
	TRACEI("bool CNAPrimaryStructure::operator==(const CNAPrimaryStructure &tSequence) const");
	// check sizes
	std::deque<nucleic_acids::SNucleotideCode>::size_type tSequenceSize = m_tSequence.size();
	if (tSequence.m_tSequence.size() != tSequenceSize)
	{
		return false;
	}
	bool bRet = true;
	std::deque<nucleic_acids::SNucleotideCode>::size_type tIndex=0;
	while ((true == bRet) && (tIndex<tSequenceSize))
	{
		bRet = (m_tSequence[tIndex] == (tSequence.m_tSequence[tIndex]));
		++tIndex;
	}
	return bRet;
}


/**
@fn bool CNAPrimaryStructure::operator!=(const CNAPrimaryStructure &tSequence) const
************************************************************************************
Tells if two sequences are different.

@param const CNAPrimaryStructure &tSequence:\n
 a sequence to compare.

@return bool:\n
 true if the sequences have different nucleotides at the same position, false
 otherwise.
*/
bool CNAPrimaryStructure::operator!=(const CNAPrimaryStructure &tSequence) const
{
	TRACEI("bool CNAPrimaryStructure::operator!=(const CNAPrimaryStructure &tSequence) const");
	// check sizes
	std::deque<nucleic_acids::SNucleotideCode>::size_type tSequenceSize = m_tSequence.size();
	if (tSequence.m_tSequence.size() != tSequenceSize)
	{
		return true;
	}
	bool bRet = true;
	std::deque<nucleic_acids::SNucleotideCode>::size_type tIndex=0;
	while ((true == bRet) && (tIndex<tSequenceSize))
	{
		bRet = (m_tSequence[tIndex] != (tSequence.m_tSequence[tIndex]));
		++tIndex;
	}
	return bRet;
}


/**
@fn void CNAPrimaryStructure::SetDescription(const std::string &strDescription)
*******************************************************************************
Sets the description of the structure.

@param const std::string &strDescription:\n
 the description string of the structure (can be multilines).

*/
void CNAPrimaryStructure::SetDescription(const std::string &strDescription)
{
	TRACEI("void CNAPrimaryStructure::SetDescription(const std::string &strDescription)");
	m_strDescription = strDescription;
}


/**
@fn std::string CNAPrimaryStructure::GetDescription(void) const
***************************************************************
Returns the description of the structure.

@return std::string:\n
 the description string of the structure (can be multilines).

*/
std::string CNAPrimaryStructure::GetDescription(void) const
{
	TRACEI("std::string CNAPrimaryStructure::GetDescription(void) const");
	return m_strDescription;
}


/**
@fn void CNAPrimaryStructure::SetName(const std::string &strName)
*****************************************************************
Sets the name of the structure.

@param const std::string &strName:\n
 the name of the structure (one-line string).

@throw std::invalid_argument:\n
 if the name contains more than a line.
*/
void CNAPrimaryStructure::SetName(const std::string &strName)
{
	TRACEI("void CNAPrimaryStructure::SetName(const std::string &strName)");
	// check for more than one line
	if (std::string::npos != strName.find("\n"))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CNAPrimaryStructure::SetName(\"" << strName << "\"): given name contains an end-of-line character!" << std::endl;
		throw std::invalid_argument(stmErrorData.str());
	}
	m_strName = strName;
}


/**
@fn std::string CNAPrimaryStructure::GetName(void) const
***************************************************************
Returns the name of the structure.

@return std::string:\n
 the name of the structure (one-line guaranteed).

*/
std::string CNAPrimaryStructure::GetName(void) const
{
	TRACEI("std::string CNAPrimaryStructure::GetName(void) const");
	return m_strName;
}


}; // namespace biocpp
}; // namespace val
