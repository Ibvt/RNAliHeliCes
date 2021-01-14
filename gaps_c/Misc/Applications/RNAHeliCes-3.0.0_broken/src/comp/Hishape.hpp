//*****************************************************************************
#ifndef __HISHAPE_H
#define __HISHAPE_H 1

/**
@file NASecondaryStructure.h
****************************
@brief NASecondaryStructure.h is the header file for
@c hishape::comp::CHishape class.

SecondaryStructure.h contains @c hishape::comp::CHishape class
declaration and methods prototypes.

@see hishape::comp::CHishape

*/

#include <map>
#include <list>
#include <vector>
//#include "NAPrimaryStructure.hpp"
#include <deque>
#include <string>

namespace hishape
{
namespace comp
{

//! Error message when the user specified an invalid stem or stem-loop number.
const std::string STR_ERROR_INVALID_STEM_OR_STEMLOOP = std::string("Error: The stem or stem-loop number specified is not valid!\n");

//! Error message when the user specified an invalid stem or stem-loop number.
const std::string STR_ERROR_CORRUPTED_STRUCTURE_DATA = std::string("Error: Stem or stem-loop data is not valid! Maybe the structure was not parsed correctly.\n");

//! Error message when a function requiers the structure to be parsed while the structure was not parsed.
const std::string STR_ERROR_NOT_PARSED_STRUCTURE = std::string("Error: Secondary structre has not been parsed before this call!\n");

//! Error message when the user specified a position outside the sequence.
const std::string STR_ERROR_INVALID_POSITION = std::string("Error: The sequence position specified is not valid!\n");

namespace na_secondary_structure
{
        // placeholder for constants
	//! The nucleotide is in a stacked pair.
	const int I_STACKED_PAIR             = -1;
	/**
	@name String format constants
	*****************************
	The following constants are used to select a string format of the secondary
	structure when calling @c ToString.
        */

}; // end namespace CHishape

/**
@class CHishape   NASecondaryStructure.h "NASecondaryStructure.hpp"
****************************
@brief CHishape implements a nucleic acids chain secondary
structure (DNA, RNA) and its basic operations.

This class works with @c C3LinksNucleotide nucleotides.

*/
class CHishape//: public CNAPrimaryStructure
{
public:
	/**
	@struct SSubstructureData
	*****************
	@brief SSubstructureData contains stem or stem-loop data.
	
	@c SSubstructureData contains stem or stem-loop data.

	*/
	/*
	struct SSubstructureData
	{
		//! Contains I_STEM_TYPE or I_STEMLOOP_TYPE (see CHishape constants).
		int iStructureNature;
		//! First segment start index.
		int iStartIndex1;
		//! Second segment start index.
		int iStartIndex2;
		//! First segment end index.
		int iEndIndex1;
		//! Second segment en index.
		int iEndIndex2;
	};
*/
	//! Default constructor.
	CHishape(void);
	CHishape(std::vector<std::string> &helixIndices);
	//CHishape(const CHishape& tHishape);
	//! Destructor.
	virtual ~CHishape(void);
	//! Assignation operator.
	virtual CHishape& operator=(const CHishape &tHishape);
	////! Inherited assignation operator.
	//virtual CNAPrimaryStructure& operator=(const val::comp::CNAPrimaryStructure &tSequence);
//	//! Concatenates 2 sequences
//	virtual CHishape operator+(const val::comp::CHishape &tRightStructure) const;



	//! Pairs two nucleotides together.
	virtual void Pair(int iFirstPosition, int iSecondPosition);
	//! Returns the pair index if a nucleotide is paired.
	virtual int GetPairedBaseIndex(int iAtPosition) const;
	//! Tells if a nucleotide is paired.
	virtual bool IsPaired(int iAtPosition) const;

	virtual void PrintPairedBaseIndices(int tNucleotidesPairsLength) const;
	
	//! returns hishape length
	virtual int GetLength(void) const;
	//! Sets the length of the hishape
	virtual void SetLength(const int iLength);
	
	virtual std::vector<std::string> GetHelixIndices(void) const;
	
//protected:
  public:
	//! Length of the hishape
	int m_iLength;
	/**
	m_tNucleotidesPairs
	*******************
	m_tNucleotidesPairs contains "sequence length + 1" elements. The first
	element is not used and should be always set to 0. m_tNucleotidesPairs
	indexation is 1-based, which means base 4 pairing data is stored in
	m_tNucleotidesPairs[4]. If that data is 0, the base 4 is not paired
	otherwise, the value is the index of the pair (still 1-based) in the
	sequence. For example, if base 4 is paired with base 17 then
	m_tNucleotidesPairs[4] == 17 and m_tNucleotidesPairs[17] == 4. If base 4
	is not paired then m_tNucleotidesPairs[4] == 0.
	
	Note: m_tNucleotidesPairs must always be GetLength() + 1 even if the
	      structure has not been parsed yet.
	*/
	std::deque<int> m_tNucleotidesPairs;
	//int m_tNucleotidesPairs[5];   
	
	std::vector<std::string> m_helixIndices;
	

	/**
	m_tNucleotidesConfig
	********************
	1-based array. First element is not used (wasted). The length of this
	array is set in ParseStructure() so it can only be used when
	m_bIsStructureParsed is true.
	*/
	std::deque<int> m_tNucleotidesConfig;


	//! type of the iterator on the array of nucleotides.
	//typedef std::deque<int>::iterator CDequeIntIterator;
	//! returns an iterator at the given position.
	//CDequeIntIterator GetPairsDataIteratorAt(int iAtPosition);

}; // class CHishape
}; // namespace comp
}; // namespace hishape
#endif //ifndef __HISHAPE_H
