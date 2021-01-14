//*****************************************************************************
#ifndef __HELIXINDEX_H
#define __HELIXINDEX_H 1

/**
@file HelixIndex.h
**********************
@brief HelixIndex.h is the header file for hishape::comp::CHelixIndex
 class.

HelixIndex.h contains hishape::comp::CHelixIndex class declaration and
methods prototypes.

@see hishape::comp::CHelixIndex
*/

#include <string>

namespace hishape
{
namespace comp
{
/**
@class CHelixIndex   HelixIndex.h "HelixIndex.hpp"
**********************
@brief CHelixIndex implements a nucleotide pair and its basic operations.

The class @c CHelixIndex implements a nucleotide pair that just contains 2
base codes.

@see CNucleotide
@see nucleic_acids::SNucleotideCode
*/
class CHelixIndex
{
public:

protected:
	float m_hiNumber;
	char m_hiType;
	//int m_firstPositionOfPair;
	
	
	// replacement costs
	static const int N				= 0;
	static const int I				= 1;
	static const int B				= 2;
	static const int H				= 3;
	static const int M				= 4;
	// delete or insert
	static const int INDEL				= 5;
	//static const int SCOREMATRIX_W                  = 6;
	//static const int SCOREMATRIX_H                  = 6;
	
	
	enum {
	  /*
	    CODE_ERROR = -1,
	    CODE_FIRST,

	    CODE_ACCEPT = CODE_FIRST,  // 0
	    CODE_AVERAGE,              // 1
	    CODE_COST,                 // 2
	    CODE_NUMLATHES,            // 3
	    CODE_DONE,                 // 4

	    CODE_LAST = CODE_DONE*/
	};


public:
	//! Default constructor.
	CHelixIndex(const float hiNumber, const char hiType/*, const int firstPositionOfPair*/);
	CHelixIndex();
	//! Destructor.
	virtual ~CHelixIndex(void);
	/*
	//! Returns the base codes as a const array of 2 unsigned char.
	std::string GetPairBaseCodes() const;
	//! Assignation operator.
	CHelixIndex& operator=(const CHelixIndex &tRightHelixIndex);
	//! == comparison operator.
	bool operator==(const CHelixIndex &tRightHelixIndex) const;
	//! != comparison operator.
	bool operator!=(const CHelixIndex &tRightHelixIndex) const;
	*/
	
	//! Sets a nucleotide of the pair.
	//void SetNucleotide(const ENucleotideIndex &eNucleotideIndex, const nucleic_acids::SNucleotideCode& tNucleotide);

	
	//! Returns the integer code associated to the pair.
	// it can be defined as virtual function or normal function
	float GetHiNumber() const;
	char GetHiType() const;
        void SetHiNumber(const float hiNumber);
        void SetHiType(const char &hiType);
	//! Returns the integer code associated to the pair.
	int GetIntCode() const;
	int GetGapCode() const;
	//int GetFristPositionOfPair() const;
	//void SetFristPositionOfPair(const int firstPositionOfPair);

	
/*	//! Returns the first nucleotide code.
	hishape::comp::nucleic_acids::SNucleotideCode GetFirstHalfNucleotideCode() const;
	//! Returns the second nucleotide code.
	hishape::comp::nucleic_acids::SNucleotideCode GetSecondHalfNucleotideCode() const;
*/	//! Returns the 2 nucleotides into a 2 characters string.
	std::string ToString() const;

}; // class CHelixIndex



/**
@namespace hishape::comp::helix_index
***************************************
@brief Contains nucleotide pairs constants.

Namespace that contains nucleotide pairs (@c CHelixIndex) names constants.

@see CHelixIndex
*/
/*
namespace helix_index
{
	//! Gap pair (--).
	const CHelixIndex GapPair = CHelixIndex(hishape::comp::nucleic_acids::Gap, hishape::comp::nucleic_acids::Gap);
}; // namespace helix_index  */
}; // namespace comp
}; // namespace hishape
#endif //ifndef __NUCLEOTIDEPAIR_H
