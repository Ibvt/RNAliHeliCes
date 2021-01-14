//*****************************************************************************
#include "Hishape.hpp"
#include <algorithm>
#include <list>
#include <sstream>
#include <stdexcept>

#include <string>


#include "DebugTools.hpp"

namespace hishape
{
namespace comp
{
/**
Default constructor
*******************
Returns a new instance of @c CHishape.

@return CHishape:\n
 a new instance of @c CHishape.
*/
CHishape::CHishape(void)//: CNAPrimaryStructure()
{
	TRACEI("CHishape::CHishape(void)");
	std::cout << "CHishape::CHishape(void)" << std::endl;
}

// store the pieces of helix indices into a vector
CHishape::CHishape(std::vector<std::string> &helixIndices)//: CNAPrimaryStructure()
{
	TRACEI("CHishape::CHishape(std::vector<std::string> &helixIndices)");
	m_helixIndices = helixIndices;
	std::cout << "CHishape::CHishape(std::vector<std::string> &helixIndices)" << std::endl;
}

/**
copy constructor: it needs operator= and destructor
****************
Returns a new instance of @c CNASecondaryStructure which is a copy of the
given instance.

@param const CNASecondaryStructure& tSecondaryStructure:\n
 a secondary structure to duplicate.

@return CNASecondaryStructure:\n
 a new instance of @c CNASecondaryStructure.
*/
/*
CHishape::CHishape(const CHishape& tHishape)
        ,m_iLength(tHishape.m_iLength)
	,m_tNucleotidesPairs(tHishape.m_tNucleotidesPairs)
	,m_helixIndices(tHishape.m_helixIndices)	
	,m_tNucleotidesConfig(tHishape.m_tNucleotidesConfig)
{
	TRACEI("CHishape::CHishape(const CHishape& tHishape)");
}
*/



/**
Destructor
**********
A standard virtual destructor.
*/
CHishape::~CHishape(void)
{
	TRACEI("CHishape::~CHishape(void)");
}


CHishape& CHishape::operator=(const CHishape &tHishape)
{
	TRACEI("CHishape& CHishape::operator=(const CHishape &tHishape)");

	// check if it's a self-assignation
	//if (this == &tHishape)
	//{
	//	return *this;
	//}
	// call inherited = operator
	//CNAPrimaryStructure::operator=(tHishape);
	
	// TODO: learn more about the knowledge about shallow and deep copy in C++, it waste me 2 days!!!
	// copy fields
        m_iLength = tHishape.m_iLength;
	m_tNucleotidesPairs = tHishape.m_tNucleotidesPairs;  // deep copy
        m_helixIndices = tHishape.m_helixIndices;
        //m_tNucleotidesConfig(tHishape.m_tNucleotidesConfig);
	
	// learn from http://www.cplusplus.com/forum/articles/18749/
	
        for( int i = 0; i < tHishape.GetLength(); ++i )
	{
            m_tNucleotidesPairs[i] = tHishape.m_tNucleotidesPairs[i];
	    m_helixIndices[i] = tHishape.m_helixIndices[i];	
	   // m_tNucleotidesConfig[i] = tHishape.m_tNucleotidesConfig[i];
	}
	
         

	//std::swap( m_tNucleotidesPairs, tHishape.m_tNucleotidesPairs );
        //std::swap( m_helixIndices, tHishape.m_helixIndices );
        //std::swap( m_tNucleotidesConfig, tHishape.m_tNucleotidesConfig);
	
	return *this;
}



/*
1-based
*/
/*
CHishape::CDequeIntIterator CHishape::GetPairsDataIteratorAt(int iAtPosition)
{
	TRACEI("CDequeIntIterator CHishape::GetPairsDataIteratorAt(int iAtPosition)");
	// keep in mind that iAtPosition is a 0-based index
	// ie. first position is 0
	int iLength = (signed)m_helixIndices.size();
	if ((iAtPosition < 0) || (iAtPosition >= iLength) )
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CHishape::GetPairsDataIteratorAt(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	// TODO: check if m_tNucleotidesPairs is the correct number
	return (m_tNucleotidesPairs.begin() + iAtPosition);
}
*/


/**
@fn void CHishape::Pair(int iFirstPosition, int iSecondPosition)
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
void CHishape::Pair(int iFirstPosition, int iSecondPosition)
{
	TRACEI("void CHishape::Pair(int iFirstPosition, int iSecondPosition)");
	
	int iLength = (signed)m_helixIndices.size();
	/*std::cout << "iLength=" << iLength << std::endl;  */
	// check given index    keep in mind that iSecondPosition can be zero
	// -1 means that iFirstPosition does not have paring base
	if ((iFirstPosition < 0) || (iSecondPosition < -1) || (iFirstPosition >= iLength) || (iSecondPosition >= iLength))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CHishape::Pair(" << iFirstPosition << ", " << iSecondPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}

	////m_bIsStructureParsed = false;
	// check if bases were already paired
	/*
	if (0 != m_tNucleotidesPairs[iFirstPosition])
	{
		// break other pair
		m_tNucleotidesPairs[m_tNucleotidesPairs[iFirstPosition]] = 0;
	}
	if (0 != m_tNucleotidesPairs[iSecondPosition])
	{
		// break other pair
		m_tNucleotidesPairs[m_tNucleotidesPairs[iSecondPosition]] = 0;
	}*/
	m_tNucleotidesPairs[iFirstPosition] = iSecondPosition;
	//debuged on Sep. 29 2011: m_tNucleotidesPairs[iSecondPosition] = iFirstPosition;
}




/**
Returns 0 if there are no pair
*/
int CHishape::GetPairedBaseIndex(int iAtPosition) const
{
	TRACEI("int CHishape::GetPairedBaseIndex(int iAtPosition) const");
	int iLength = (signed)m_helixIndices.size();
	// check given index
	if ((iAtPosition < 0) || (iAtPosition >= iLength))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CHishape::GetPairedBaseIndex(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	return m_tNucleotidesPairs[iAtPosition];
}


/**
Tells if a nucleotide is paired.
*/
bool CHishape::IsPaired(int iAtPosition) const
{
	TRACEI("bool CHishape::IsPaired(int iAtPosition) const");
	/*
	if (false == m_bIsStructureParsed)
	{
		throw std::logic_error("CHishape::IsPaired: " + STR_ERROR_NOT_PARSED_STRUCTURE);  
	}*/
	int iLength = (signed)m_helixIndices.size();
	if ((iAtPosition < 0) || (iAtPosition >= iLength))
	{
		std::ostringstream stmErrorData;
		stmErrorData << "CHishape::IsPaired(" << iAtPosition << "): " << STR_ERROR_INVALID_POSITION << std::endl;
		throw std::out_of_range(stmErrorData.str());
	}
	
	//std::cout << "Hishape::I_STACKED_PAIR=" << hishape::comp::na_secondary_structure::I_STACKED_PAIR << "  GetPairedBaseIndex(iAtPosition)=" << GetPairedBaseIndex(iAtPosition) << "  iAtPosition=" << iAtPosition << std::endl;
	bool bRet = false;
	if (hishape::comp::na_secondary_structure::I_STACKED_PAIR != m_tNucleotidesPairs[iAtPosition])
	{       
		bRet = true;
	}
	return bRet;
}


/**
Tells if a nucleotide is paired.
*/
// TODO
void CHishape::PrintPairedBaseIndices(int tNucleotidesPairsLength) const
{
  std::cout << "Hishape::mydeque contains: " << m_helixIndices.size() << std::endl;
  //int i;
  //for (i=0; i<m_helixIndices.size(); i++)
  //    std::cout << "i=" << i << "value="<< m_tNucleotidesPairs.at(i) << std::endl;
}


	//! returns hishape length
	int CHishape::GetLength(void) const
	{
	    return m_helixIndices.size(); 
	}
	//! Sets the length of the hishape
	void CHishape::SetLength(const int iLength)
	{
	    m_iLength = iLength;
	}

	//! Sets the length of the hishape
	std::vector<std::string> CHishape::GetHelixIndices(void) const
	{
	    return m_helixIndices;
	}
	
	





}; // namespace comp
}; // namespace hishape
