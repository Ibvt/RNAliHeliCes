//*****************************************************************************
#ifndef __NUCLEOTIDEPAIR_H
#define __NUCLEOTIDEPAIR_H 1

/**
@file NucleotidePair.h
**********************
@brief NucleotidePair.h is the header file for val::biocpp::CNucleotidePair
 class.

NucleotidePair.h contains val::biocpp::CNucleotidePair class declaration and
methods prototypes.

@see val::biocpp::CNucleotidePair

@author Valentin GUIGNON
@version 1.0
@date 24/06/2007
*/

#include "BioCPP.hpp"
#include <string>

namespace val
{
namespace biocpp
{
/**
@class CNucleotidePair   NucleotidePair.h "NucleotidePair.hpp"
**********************
@brief CNucleotidePair implements a nucleotide pair and its basic operations.

The class @c CNucleotidePair implements a nucleotide pair that just contains 2
base codes.

@see CNucleotide
@see nucleic_acids::SNucleotideCode

@author Valentin GUIGNON
@version 1.0
@date 22/08/2004
*/
class CNucleotidePair
{
public:
	/**
	@enum ENucleotideIndex
	**********************
	@brief A simple enumeration for nucleotide indexing inside the pair.

	A simple enumeration for nucleotide indexing inside the pair. There are
	only two allowed values: @c FirstNucleotide and @c SecondNucleotide.
	@c FirstNucleotide has the value 0 and @c SecondNucleotide has the value 1.

	@see operator[]
	*/
	enum ENucleotideIndex
	{
		//! The first (left) nucleotide of the pair
		FirstNucleotide = 0,
		//! the second (right) nucleotide of the pair
		SecondNucleotide = 1
	};

protected:
	//! The nucleotide codes of this nucleotide.
	nucleic_acids::SNucleotideCode m_atNucleotideCodes[2];
	//! The associated int code for this pair
	int m_iIntCode;

public:
	//! Default constructor.
	CNucleotidePair(const nucleic_acids::SNucleotideCode &tFirstNucleotide=nucleic_acids::Unknown, const nucleic_acids::SNucleotideCode &tSecondNucleotide=nucleic_acids::Unknown);
	//! Copy constructor.
	CNucleotidePair(const CNucleotidePair &tNucleotidePair);
	//! Destructor.
	virtual ~CNucleotidePair(void);
	//! Returns the base codes as a const array of 2 unsigned char.
	std::string GetPairBaseCodes() const;
	//! Assignation operator.
	CNucleotidePair& operator=(const CNucleotidePair &tRightNucleotidePair);
	//! == comparison operator.
	bool operator==(const CNucleotidePair &tRightNucleotidePair) const;
	//! != comparison operator.
	bool operator!=(const CNucleotidePair &tRightNucleotidePair) const;
	//! Sets a nucleotide of the pair.
	void SetNucleotide(const ENucleotideIndex &eNucleotideIndex, const nucleic_acids::SNucleotideCode& tNucleotide);
	//! Sets the nucleotides of the pair.
	void SetPair(const nucleic_acids::SNucleotideCode &tFirstNucleotide, const nucleic_acids::SNucleotideCode &tSecondNucleotide);
	//! [] operator to access to the bases value (const version).
	const nucleic_acids::SNucleotideCode& operator[](const ENucleotideIndex &eNucleotideIndex) const;
	//! Returns the integer code associated to the pair.
	int GetIntCode() const;
	//! Returns the integer code of the first nucleotide code.
	int GetFirstHalfIntCode() const;
	//! Returns the integer code of the second nucleotide code.
	int GetSecondHalfIntCode() const;
	//! Returns the first nucleotide code.
	val::biocpp::nucleic_acids::SNucleotideCode GetFirstHalfNucleotideCode() const;
	//! Returns the second nucleotide code.
	val::biocpp::nucleic_acids::SNucleotideCode GetSecondHalfNucleotideCode() const;
	//! Returns the 2 nucleotides into a 2 characters string.
	std::string ToString() const;
	//! Returns a pair of nucleotides made from 1 nucleotide of each given pair.
	static CNucleotidePair FusionNucleotides(const CNucleotidePair &tFirstHalf, const CNucleotidePair &tSecondHalf);

}; // class CNucleotidePair


/**
@namespace val::biocpp::nucleotide_pair
***************************************
@brief Contains nucleotide pairs constants.

Namespace that contains nucleotide pairs (@c CNucleotidePair) names constants.

@see CNucleotidePair
@see SNucleotideCode
@see BioCPP.h
*/
namespace nucleotide_pair
{
	//! Gap pair (--).
	const CNucleotidePair GapPair = CNucleotidePair(val::biocpp::nucleic_acids::Gap, val::biocpp::nucleic_acids::Gap);
}; // namespace nucleotide_pair
}; // namespace biocpp
}; // namespace val
#endif //ifndef __NUCLEOTIDEPAIR_H
