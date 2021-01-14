//*****************************************************************************
#ifndef __NAPRIMARYSTRUCUTRE_H
#define __NAPRIMARYSTRUCUTRE_H 1

/**
@file NAPrimaryStructure.h
**************************
@brief NAPrimaryStructure.h is the header file for
 val::biocpp::CNAPrimaryStructure class.

PrimaryStructure.h contains val::biocpp::CNAPrimaryStructure class declaration
and methods prototypes.

@see val::biocpp::CNAPrimaryStructure

@author Valentin GUIGNON
@version 2.0
@date 08/02/2007
*/

#include <deque>
#include <string>
#include <iterator>
#include "BioCPP.hpp"

namespace val
{
namespace biocpp
{
//! Error message when the user specified an invalid sequence length.
const std::string STR_ERROR_INVALID_LENGTH = std::string("Error: The sequence length specified is not valid!\n");

//! Error message when the user specified a position outside the sequence.
const std::string STR_ERROR_INVALID_POSITION = std::string("Error: The sequence position specified is not valid!\n");

/**
@class CNAPrimaryStructure   NAPrimaryStructure.h "NAPrimaryStructure.hpp"
**************************
@brief CNAPrimaryStructure implements a nucleic acids chain primary structure
(DNA, RNA) and its basic operations.

This class implements a nucleic acids chain primary structure (DNA, RNA) and
its basic operations. Nucleic acids used by the class are nucleic_acids::SNucleotideCode and any
nucleic acid of the structure can be accessed by its position using the []
operator.

@note when the term "position" is used, it means a 1-based index (items are
 indexed starting from 1); when the term "offset" is used, it means a 0-based
 index.

@see nucleic_acids::SNucleotideCode
*/
//+FIXME: put a better description and add code examples
//+FIXME: counting methods (GetNucleotideCount(base code), GetSubSequenceCount(sub sequence))
class CNAPrimaryStructure
{
public:
	//! Default constructor.
	CNAPrimaryStructure(void);
	//! Copy constructor.
	CNAPrimaryStructure(const CNAPrimaryStructure &tPrimaryStructure);
	//! Constructor from a string.
	CNAPrimaryStructure(const std::string &strSequence, int iFromOffset=0, int iLength=-1);
	//! Destructor.
	virtual ~CNAPrimaryStructure(void);
	//! Adds a nucleotide on the 3' end of the sequence.
	virtual void AddOn3Prime(const nucleic_acids::SNucleotideCode &tNucleotide);
	//! Adds a nucleotide on the 5' end of the sequence.
	virtual void AddOn5Prime(const nucleic_acids::SNucleotideCode &tNucleotide);
	//! Appends a nucleotide on the 3' end of the sequence.
	virtual void Append(const nucleic_acids::SNucleotideCode &tNucleotide);
	//! Appends a sequence to this sequence on the 3' end.
	virtual int AppendSequence(const std::string &strSequence, int iFromOffset=0, int iAppendLength=-1);
	//! Appends a sequence to this sequence on the 3' end.
	virtual int AppendSequence(const CNAPrimaryStructure &tSequence, int iSequenceOffset=0, int iAppendLength=-1);
	//! Inserts iCopyCount copies of tBaseCode after the given position.
	virtual void InsertAt(const nucleic_acids::SNucleotideCode &tNucleotide, int iAtPosition, int iCopyCount=1);
	//! Inserts a sequence of nucleotides at the given position.
	virtual int InsertSequenceAt(const std::string &strSequence, int iAtPosition, int iFromOffset=0, int iInsertLength=-1, int iCopyCount=1);
	//! Inserts a sequence of nucleotides at the given position.
	virtual void InsertSequenceAt(const CNAPrimaryStructure &tSequence, int iAtPosition, int iFromOffset=0, int iInsertLength=-1, int iCopyCount=1);
	//! Removes iNucleotideCount nucleotides on the 3' end.
	virtual void RemoveOn3Prime(int iNucleotideCount=1);
	//! Removes iNucleotideCount nucleotides on the 5' end.
	virtual void RemoveOn5Prime(int iNucleotideCount=1);
	//! Removes iNucleotideCount nucleotides at the given position.
	virtual void RemoveAt(int iAtPosition, int iNucleotideCount=1);
	//! Removes iNucleotideCount nucleotides before the given nucleotide.
	virtual void RemoveBefore(int iAtPosition, int iNucleotideCount=1);
	//! Removes iNucleotideCount nucleotides after the given nucleotide.
	virtual void RemoveAfter(int iAtPosition, int iNucleotideCount=1);
	//! Returns the iAtPosition-th nucleotide (starts from 1).
	virtual const nucleic_acids::SNucleotideCode& GetAt(int iAtPosition) const;
	//! Changes the nucleotide at the given position.
	virtual void SetAt(int iAtPosition, const nucleic_acids::SNucleotideCode& tNewNucleotide);
	//! Returns the sequence length (nucleotides count).
	virtual int GetLength(void) const;
	//! Replace all the Us in the sequence by Ts.
	virtual void TranscriptToDNA(void);
	//! Replace all the Ts in the sequence by Us.
	virtual void TranscriptToRNA(void);
	//! Reverts 3' and 5' ends.
	virtual CNAPrimaryStructure RevertSequence(void) const;
	//! Returns a nucleotide-complementary sequence of this sequence.
	virtual CNAPrimaryStructure ComplementSequence(void) const;
	//! Returns a string that contains the sequence as an array of characters
	virtual std::string ToString(void) const; 
	//! Returns the nucleotide at the given position (1-based index).
	virtual const nucleic_acids::SNucleotideCode& operator[](int iAtPosition) const;
	//! Assignation operator.
	virtual CNAPrimaryStructure& operator=(const val::biocpp::CNAPrimaryStructure &tSequence);
	//! Concatenates 2 sequences
	virtual CNAPrimaryStructure operator+(const val::biocpp::CNAPrimaryStructure &tRightSequence) const;
	//! Concatenates the second sequence to the first one
	virtual CNAPrimaryStructure& operator+=(const val::biocpp::CNAPrimaryStructure &tRightSequence);
	//! Returns the complementary sequence (same as calling ComplementSequence)
	virtual CNAPrimaryStructure operator~(void) const;
	//! tells if 2 sequences are identical (same bases at the same positions)
 	virtual bool operator==(const val::biocpp::CNAPrimaryStructure &tSequence) const;
	//! tells if 2 sequences are different
 	virtual bool operator!=(const val::biocpp::CNAPrimaryStructure &tSequence) const;
	//! returns structure description
	std::string GetDescription(void) const;
	//! Sets the description of the structure
	void SetDescription(const std::string &strDescription);
	//! returns structure name
	std::string GetName(void) const;
	//! Sets the name of the structure
	void SetName(const std::string &strName);


protected:
	//! type of the const_iterator on the array of nucleotides.
	typedef std::deque<nucleic_acids::SNucleotideCode>::const_iterator CStructureConstIterator;
	//! returns a const_iterator at the given position.
	const CStructureConstIterator GetConstIteratorAt(int iAtPosition) const;
	//! Description of the structure
	std::string m_strDescription;
	//! Name of the structure
	std::string m_strName;

private:
 	//! stores the nucleotides sequence as a double-ended queue of nucleic_acids::SNucleotideCode*.
	std::deque<nucleic_acids::SNucleotideCode> m_tSequence;
	//! Deletes all the nucleic_acids::SNucleotideCode allocated by this instance.
	void ClearSequence(void);
	//! type of the iterator on the array of nucleotides.
	typedef std::deque<nucleic_acids::SNucleotideCode>::iterator CStructureIterator;
	//! returns an iterator at the given position.
	CStructureIterator GetIteratorAt(int iAtPosition);
}; // class CNAPrimaryStructure
}; // namespace biocpp
}; // namespace val
#endif //ifndef __NAPRIMARYSTRUCUTRE_H
