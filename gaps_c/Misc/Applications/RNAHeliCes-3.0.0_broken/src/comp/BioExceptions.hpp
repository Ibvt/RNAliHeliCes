//*****************************************************************************
#ifndef __BIOEXCEPTIONS_H
#define __BIOEXCEPTIONS_H 1

/**
@file BioExceptions.h
*********************
@brief BioExceptions.h is a library of biological exceptions.

BioExceptions.h library contains exceptions thrown by BioCPP tools. These
exceptions are here to bring more details when a biological exception occurs.

@see val::biocpp::bioexceptions

@author Valentin GUIGNON
@version 1.0
@date 30/06/2004
*/

#include <stdexcept>
#include <string>
#include <sstream>

namespace val
{
namespace biocpp
{
/**
@namespace val::biocpp::bioexceptions
*************************************
@brief bioexceptions contains biological exceptions.

bioexceptions contains exceptions types that deal with BioCPP-related
problems.

@see BioExceptions.h
*/
namespace bioexceptions
{
	/**
	@class CInvalidLink   BioExceptions.h "BioExceptions.hpp"
	*******************
	@brief CInvalidLink is thrown on invalid links between molecules.

	CInvalidLink exception occurs when something tries to link together things
	that can't be linked.
	*/
	class CInvalidLink: public std::runtime_error
	{
	protected:
		//! Pointer to the first element involved in the link try
		void* m_pvFirstElement;
		//! Pointer to the second element involved in the link try
		void* m_pvSecondElement;
		//! exception message string
		std::string m_strMessage;
	public:
		CInvalidLink(const std::string& strMessage, void* pvFirstElement, void* pvSecondElement)
		: std::runtime_error(strMessage)
		, m_pvFirstElement(pvFirstElement)
		, m_pvSecondElement(pvSecondElement)
		, m_strMessage(strMessage)
		{
			std::ostringstream stmPointerData;
			stmPointerData << m_pvFirstElement;
			m_strMessage += stmPointerData.str();
			m_strMessage += ", ";
			stmPointerData << m_pvSecondElement;
			m_strMessage += stmPointerData.str();
		}
		virtual ~CInvalidLink() throw()
		{
		}
		virtual const char *what() const throw()
		{
			return m_strMessage.c_str();
		}
	}; // class CInvalidLink


	/**
	@class CInvalidCode   BioExceptions.h "BioExceptions.hpp"
	*******************
	@brief CInvalidCode is thrown when an invalid code is met.

	CInvalidCode exception reports that a code was found (at a 
	specific position in a sequence) with an invalid value.
	The invalid code can be obtained using the getter @c GetCode.
	The position can be obtained using the getter @c GetPosition.
	*/
	class CInvalidCode: public std::range_error
	{
	protected:
		//! Position where the invalid code was found
		int m_iPosition;
		//! Invalid code value met
		char m_cCode;
		//! exception message string
		std::string m_strMessage;
	public:
		CInvalidCode(const std::string& strMessage, const char &cCode='\0', const int &iAtPosition=-1)
		: std::range_error(strMessage)
		, m_iPosition(iAtPosition)
		, m_cCode(cCode)
		, m_strMessage(strMessage)
		{
			if (isprint(cCode))
			{
				m_strMessage += " ";
				m_strMessage += m_cCode;
				if (0 <= iAtPosition)
				{
					std::ostringstream stmPosition;
					m_strMessage += ", ";
					stmPosition << m_iPosition;
					m_strMessage += stmPosition.str();
				}
			}
		}
		virtual ~CInvalidCode() throw()
		{
		}
		virtual const char *what() const throw()
		{
			return m_strMessage.c_str();
		}
		/**
		@fn int GetPosition() const throw()
		***********************************
		Returns the invalid code (char) position or -1 if unknown.

		@return char:
		 the invalid code (char) position.
		*/
		int GetPosition() const throw()
		{
			return m_iPosition;
		}
		/**
		@fn char GetCode() const throw()
		********************************
		Returns the invalid code (char) or '\0' if unknown.

		@return char:
		 the invalid code.
		*/
		char GetCode() const throw()
		{
			return m_cCode;
		}
	}; // class CInvalidCode


	/**
	@class CInvalidNucleotide   BioExceptions.h "BioExceptions.hpp"
	*************************
	@brief CInvalidNucleotide is thrown when an invalid nucleotide is met.

	CInvalidNucleotide exception reports that a nucleotide was found (at a 
	specific position) with an invalid value.
	The invalid nucleotide code can be obtained using the getter @c GetCode.
	The position can be obtained using the getter @c GetPosition.
	*/
	class CInvalidNucleotide: public CInvalidCode
	{
	public:
		CInvalidNucleotide(const std::string& message, const char &cNucleotideCode='\0', const int &iAtPosition=-1)
		: CInvalidCode(message, cNucleotideCode, iAtPosition)
		{
		}
		virtual ~CInvalidNucleotide() throw()
		{
		}
	}; // class CInvalidNucleotide


	/**
	@class CInvalidAminoAcid   BioExceptions.h "BioExceptions.hpp"
	************************
	@brief CInvalidAminoAcid is thrown when an invalid amino acid is met.

	CInvalidAminoAcid exception reports that an amino acid was found (at a 
	specific position) with an invalid value.
	The invalid amino acid code can be obtained using the getter @c GetCode.
	The position can be obtained using the getter @c GetPosition.
	*/
	class CInvalidAminoAcid: public CInvalidCode
	{
	public:
		CInvalidAminoAcid(const std::string& message, const char &cNucleotideCode='\0', const int &iAtPosition=-1)
		: CInvalidCode(message, cNucleotideCode, iAtPosition)
		{
		}
		virtual ~CInvalidAminoAcid() throw()
		{
		}
	}; // class CInvalidAminoAcid


	/**
	@class CAlgorithmLimitReached   BioExceptions.h "BioExceptions.hpp"
	*****************************
	@brief CAlgorithmLimitReached is thrown when an algorithm reached its limits.

	CAlgorithmLimitReached is thrown when an algorithm detects its computation
	needs will exceed its fixed limit.
	*/
	class CAlgorithmLimitReached: public std::runtime_error
	{
	protected:
		//! exception message string
		std::string m_strMessage;
	public:
		CAlgorithmLimitReached(const std::string& strMessage)
		: std::runtime_error(strMessage)
		, m_strMessage(strMessage)
		{
		}
		virtual ~CAlgorithmLimitReached() throw()
		{
		}
		virtual const char *what() const throw()
		{
			return m_strMessage.c_str();
		}
	}; // class CAlgorithmLimitReached


}; // namespace bioexceptions
}; // namespace biocpp
}; // namespace val
#endif //ifndef __BIOEXCEPTIONS_H
