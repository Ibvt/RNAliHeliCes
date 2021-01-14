//*****************************************************************************
#ifndef __BIOCPP_H
#define __BIOCPP_H

/*
BioCPP library.
Copyright (C) 2004, Valentin GUIGNON

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Contact:
 guignonv@yahoo.fr

Address:
 Laboratoire de Biologie informatique et Théorique
 DIRO, Université de Montréal,
 C.P. 6128, succursale Centre-ville,
 Montreal (Québec)
 Canada, H3C 3J7
*/


/**
@file BioCPP.h
**************
@brief BioCPP.h contains tools, consts and typedefs for the BioCPP environment.

BioCPP.h is a header file that contains the main BioCPP constant definitions,
type definitions, structures, macros and also provides some basic tools such as
converters for various types and more. It defines the namespace val::biocpp.

@see val::biocpp

@author Valentin GUIGNON
@version 1.0
@date 18/07/2004
*/
#include <string>
#include "BioExceptions.hpp"

// defines NULL if needed
#ifndef NULL
 #define NULL 0
#endif


/**
@namespace val
**************
@brief Val's namespace.

val namespace contains the tools developed by Valentin GUIGNON.
*/
namespace val
{


/**
@namespace val::biocpp
**********************
@brief biocpp is a set of bioinformatics tools.

biocpp is a namespace that contains various bioinformatics tools such as
nucleotides and sequences implementation.

@see BioCPP.h
@see nucleic_acids
@see amino_acids
*/
namespace biocpp
{
	//! Error message when an invalid nucleotide is found.
	const std::string STR_ERROR_INVALID_NUCLEOTIDE = std::string("Error: an invalid nucleotide code was found (code, offset): \n");

	//! Error message when an invalid amino acid is found.
	const std::string STR_ERROR_INVALID_AMINOACID = std::string("Error: an invalid amino acid code was found (code, offset): \n");


	/**
	@namespace val::biocpp::nucleic_acids
	*************************************
	@brief Contains nucleic acids definitions.

	Namespace that contains nucleic acids names and groups of nucleic acids.

	@see SNucleotideCode
	@see BioCPP.h
	*/
	namespace nucleic_acids
	{
		/**
		@struct SNucleotideCode
		***********************
		@brief A 1 char structure that represents a simple nucleotide type.

		The nucleotide code structure contains only one member variable which
		is of type char. Basically, it means a @c SNucleotideCode can be
		considered just like being a char. However the definition and the use
		of this structure prevents mistyping errors on nucleotides codes. Using
		this structure and its associated constants is SAFER than using chars
		(and doesn't bring any additionnal cost) because the compiler will be
		able to detect many logic errors on nucleotide codes!

		@see AI_BASE_TO_INT
		@see AT_INT_TO_BASE

		@version 1.0
		@date 08/07/2004
		*/
/*
//+FIXME:
Nucléotides d'ARN dérivés:
Pseudouridine (W) (psi)
Dihydrouridine (D)
Ribothymidine (T)
4-Thiouridine (s^4U)
3-Méthylecytidine (m^3C)
N^4-Acétylcytidine (ac^4C)
Lysidine (L)
1-Méthyladénosine (m^1A)
N^6-Isopentényladénosine (i^6A)
Inosine (I)
N^7-Méthylguanosine (m^7G)
N^2,N^2-Diméthyleguanosine (m^2_2G)
Wyosine (Wyo)
*/
		struct SNucleotideCode
		{
			// WARNING: do not add any other variable member to this class!
			//          And do not add any constructor, destructor or virtual
			//          methodes to insure that one instance will always be
			//          only 1 byte in memory (== sizeof(char)).
			//          Changing the size of the instance would make this
			//          structure unusable by most of the functions defined
			//          in this file and would also make this structure
			//          useless. Basically... DON'T TOUCH THIS STRUCTURE! :-)
			//! The base code as an unsigned char.
			unsigned char m_cBaseCode;
			// Some basic operations that don't need to be explained...
			inline SNucleotideCode& operator=(const nucleic_acids::SNucleotideCode &tRightNucleotide)
			{
				m_cBaseCode = tRightNucleotide.m_cBaseCode;
				return *this;
			};
			inline bool operator==(const val::biocpp::nucleic_acids::SNucleotideCode &tRightNucleotide)
			{
				return (m_cBaseCode == tRightNucleotide.m_cBaseCode);
			};
			inline bool operator!=(const val::biocpp::nucleic_acids::SNucleotideCode &tRightNucleotide)
			{
				return (m_cBaseCode != tRightNucleotide.m_cBaseCode);
			};
			inline bool operator==(const val::biocpp::nucleic_acids::SNucleotideCode &tRightNucleotide) const
			{
				return (m_cBaseCode == tRightNucleotide.m_cBaseCode);
			};
			inline bool operator!=(const val::biocpp::nucleic_acids::SNucleotideCode &tRightNucleotide) const
			{
				return (m_cBaseCode != tRightNucleotide.m_cBaseCode);
			};
		};

		const SNucleotideCode	Adenine		= {'A'}; //!< Adenine; letter code: @b A
		const SNucleotideCode	Cytosine	= {'C'}; //!< Cytosine; letter code: @b C
		const SNucleotideCode	Guanine		= {'G'}; //!< Guanine; letter code: @b G
		const SNucleotideCode	Thymine		= {'T'}; //!< Thymine; letter code: @b T
		const SNucleotideCode	Uridine		= {'U'}; //!< Uridine; letter code: @b U
		const SNucleotideCode	Purine		= {'R'}; //!< A or G; letter code: @b R
		const SNucleotideCode	Pyrimidine	= {'Y'}; //!< C, T or U; letter code: @b Y
		const SNucleotideCode	Amino		= {'M'}; //!< A or C; letter code: @b M
		const SNucleotideCode	Keto		= {'K'}; //!< G, T or U; letter code: @b K
		const SNucleotideCode	Strong		= {'S'}; //!< G or C (3 H-bounds); letter code: @b S
		const SNucleotideCode	Weak		= {'W'}; //!< A, T or U (2 H-bounds); letter code: @b W
		const SNucleotideCode	NotA		= {'B'}; //!< C, G, T or U but not A; letter code: @b B (B follows A)
		const SNucleotideCode	NotC		= {'D'}; //!< G, A, T or U but not C; letter code: @b D (D follows C)
		const SNucleotideCode	NotG		= {'H'}; //!< A, C, T or U but not G; letter code: @b H (H follows G)
		const SNucleotideCode	NotT		= {'V'}; //!< G, C or A but not T nor U; letter code: @b V (V follows T and U)
		const SNucleotideCode	NotU		= {'V'}; //!< G, C or A but not T nor U; letter code: @b V (V follows U)
		const SNucleotideCode	Any			= {'N'}; //!< A, C, G, T or U; letter code: @b N (aNy)
		const SNucleotideCode	Gap			= {'-'}; //!< Gap; letter code: @b -
		const SNucleotideCode	Unknown		= {'?'}; //!< Unknown base type; letter code: @b ?
		const SNucleotideCode	Null		= {'\0'}; //!< End sequence character; letter code: @b NULL
	}; // namespace nucleic_acids


    // placeholder

	/**
	@var I_ASCII_CHAR_COUNT
	***********************
	@brief Number of ASCII characters.

	Number of ASCII characters. This constant is used by convertion tables that
	work with characters.

	@see AI_BASE_TO_INT
	@see AT_INT_TO_BASE
	*/
	const int I_ASCII_CHAR_COUNT	= 0x0100;


	/**
	@var I_NUCLEOTIDES_COUNT
	************************
	@brief Number of nucleotide codes.

	Number of nucleotide codes. The count includes groups of nucleotides (like
	Keto, etc.). This constant is used by convertion tables that work with
	nucleotide bases. The count reflects the "maximum value + 1" a nucleotide
	can take when it is converted to an integer code (NOT an index!).

	@see AI_BASE_TO_INT
	@see AT_INT_TO_BASE
	*/
	const int I_NUCLEOTIDES_COUNT	= 0x0010;


	//! Converts a base code into a nucleotide code.
	nucleic_acids::SNucleotideCode& CharToSNucleotideCode(char &cBaseCode);


	//! Converts an array of char into an array of nucleotide codes.
	nucleic_acids::SNucleotideCode* CharArrayToSNucleotideCodeArray(char *acBaseCodes);


	/**
	@var AT_CHAR_TO_SNUC
	********************
	@brief Converts a base code into a nucleotide (case unsensitive).

	Converts a single character representing a base code into its associated
	nucleotide code. This convertion table is case unsensitive. If the given
	base code is unknown, it will return nucleic_acids::Unknown.

	Example:
	@code
	nucleic_acids::SNucleotideCode tMyBase;
	tMyBase = AT_CHAR_TO_SNUC['a']; // lower case
	assert(nucleic_acids::Adenine == tMyBase); // ok
	std::cout << tMyBase.m_cBaseCode << std::endl; // displays "A", upper case
	@endcode

	@note @c CharToSNucleotideCodeSafe function is safer to use.

	@see nucleic_acids::SNucleotideCode
	@see CharToSNucleotideCodeSafe
	@see CharToSNucleotideCode
	@see CharArrayToSNucleotideCodeArray
	@see CheckBaseCodes
	@see CheckAndUpcaseBaseCodes
	@see SNucleotideCodeToChar
	@see SNucleotideCodeArrayToCharArray

	@hideinitializer
	*/
	const nucleic_acids::SNucleotideCode AT_CHAR_TO_SNUC[I_ASCII_CHAR_COUNT] =
	{
		nucleic_acids::Null, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // 00h-0Fh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // 10h-1Fh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Gap, nucleic_acids::Unknown, nucleic_acids::Unknown, // 20h-2Fh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // 30h-3Fh
		nucleic_acids::Unknown, nucleic_acids::Adenine, nucleic_acids::NotA, nucleic_acids::Cytosine, nucleic_acids::NotC, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Guanine, nucleic_acids::NotG, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Keto, nucleic_acids::Unknown, nucleic_acids::Amino, nucleic_acids::Any, nucleic_acids::Unknown, // 40h-4Fh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Purine, nucleic_acids::Strong, nucleic_acids::Thymine, nucleic_acids::Uridine, nucleic_acids::NotT, nucleic_acids::Weak, nucleic_acids::Unknown, nucleic_acids::Pyrimidine, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // 50h-5Fh
		nucleic_acids::Unknown, nucleic_acids::Adenine, nucleic_acids::NotA, nucleic_acids::Cytosine, nucleic_acids::NotC, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Guanine, nucleic_acids::NotG, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Keto, nucleic_acids::Unknown, nucleic_acids::Amino, nucleic_acids::Any, nucleic_acids::Unknown, // 40h-4Fh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Purine, nucleic_acids::Strong, nucleic_acids::Thymine, nucleic_acids::Uridine, nucleic_acids::NotT, nucleic_acids::Weak, nucleic_acids::Unknown, nucleic_acids::Pyrimidine, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // 50h-5Fh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // 80h-8Fh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // 90h-9Fh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // A0h-AFh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // B0h-BFh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // C0h-CFh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // D0h-DFh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // E0h-EFh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown  // F0h-FFh
	};


	//! Converts a base code into a nucleotide code (case unsensitive).
	const nucleic_acids::SNucleotideCode& CharToSNucleotideCodeSafe(const char &cBaseCode);


	//! Checks if the chars can be converted or rise an exception.
	void CheckBaseCodes(const char *acBaseArray, const int &iArraySize=-1);


	//! Upcases the array and check it can be converted or rise an exception.
	void CheckAndUpcaseBaseCodes(char *acBaseArray, const int &iArraySize=-1);
	//+FIXME: desc
	void CheckAndUpcaseBaseCodes(std::string strBaseArray, const int &iArraySize=-1);

	//! Converts a nucleotide code into a base code (char).
	char& SNucleotideCodeToChar(nucleic_acids::SNucleotideCode& tNucleotideCode);


	//! Converts a nucleotide code into a base code (const char).
	const char& SNucleotideCodeToChar(const nucleic_acids::SNucleotideCode& tNucleotideCode);


	//! Converts a nucleotide code into an index for lookup tables.
	const unsigned char& SNucleotideCodeToTableIndex(const nucleic_acids::SNucleotideCode& tNucleotideCode);


	//! Converts an array of nucleotide codes into an array of char.
	char* SNucleotideCodeArrayToCharArray(nucleic_acids::SNucleotideCode* atNucleotideCode);


	/**
	@var AI_BASE_TO_INT
	*******************
	@brief Converts a base code code into a specific integer.

	This lookup table can convert a base code or an @c SNucleotideCode
	(converted to a char) into an integer value. It is case unsensitive	and
	faster than any	procedure for the same job. Each of the 4 bases has an
	associated bit as follow:
	- G = 1     = 0001b (bit 0)
	- A = 2     = 0010b (bit 1)
	- C = 4     = 0100b (bit 2)
	- T = U = 8 = 1000b (bit 3)

	Using these bit organisation, there are @c I_NUCLEOTIDES_COUNT different
	codes. Associated decimal values for each base code are:
	- Gap         (-) = 0
	- Unknown     (?) = 0
	- Guanine     (G) = 1
	- Adenine     (A) = 2
	- Purine      (R) = 3
	- Cytosine    (C) = 4
	- Strong      (S) = 5
	- Amino       (M) = 6
	- NotT        (V) = 7
	- NotU        (V) = 7
	- Thymine     (T) = 8
	- Uridine     (U) = 8
	- Keto        (K) = 9
	- Weak        (W) = 10
	- NotC        (D) = 11
	- Pyrimidine  (Y) = 12
	- NotA        (B) = 13
	- NotG        (H) = 14
	- Any         (N) = 15

	@note T and U have the same integer value.
	@note groups of bases will have several bit sets.
	@note changing the maximum value returned by this table implies changes
	 into AT_INT_TO_BASE and in some other tables.
	@note @c BaseCodeToIntCode is safer to use.

	@see BaseCodeToIntCode
	@see SNucleotideCodeToIntCode
	@see AT_INT_TO_BASE
	@see nucleic_acids::SNucleotideCode
	@see I_NUCLEOTIDES_COUNT

	@hideinitializer
	*/
	const int AI_BASE_TO_INT[I_ASCII_CHAR_COUNT] =
	{
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // 00h-0Fh
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // 10h-1Fh
		//                                                                            -
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // 20h-2Fh
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // 30h-3Fh
		//    A     B     C     D     E     F     G     H     I     J     K     L     M     N     O
		0x00, 0x02, 0x0D, 0x04, 0x0B, 0x00, 0x00, 0x01, 0x0E, 0x00, 0x00, 0x09, 0x00, 0x06, 0x0F, 0x00, // 40h-4Fh
		//P   Q     R     S     T     U     V     W     X     Y     Z
		0x00, 0x00, 0x03, 0x05, 0x08, 0x08, 0x07, 0x0A, 0x00, 0x0C, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // 50h-5Fh
		//    a     b     c     d     e     f     g     h     i     j     k     l     m     n     o
		0x00, 0x02, 0x0D, 0x04, 0x0B, 0x00, 0x00, 0x01, 0x0E, 0x00, 0x00, 0x09, 0x00, 0x06, 0x0F, 0x00, // 60h-6Fh
		//p   q     r     s     t     u     v     w     x     y     z
		0x00, 0x00, 0x03, 0x05, 0x08, 0x08, 0x07, 0x0A, 0x00, 0x0C, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // 70h-7Fh
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // 80h-8Fh
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // 90h-9Fh
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // A0h-AFh
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // B0h-BFh
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // C0h-CFh
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // D0h-DFh
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // E0h-EFh
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00  // F0h-FFh
	};


	//! Converts a base code into its associated integer code.
	const int& BaseCodeToIntCode(const char &cBaseCode);


	//! Converts a nucleotide code into its associated integer code.
	const int& SNucleotideCodeToIntCode(const nucleic_acids::SNucleotideCode &tNucleotideCode);


	/**
	@var AT_INT_TO_BASE
	*******************
	@brief Converts an integer code into its associated base code.

	This lookup table can convert a integer code into an SNucleotideCode
	(or a base letter code using SNUC_TO_CHAR macro). It is faster than any
	procedure for this purpose.

	@note as T and U have the same integer value, only T is returned, never U.
	@note integer codes used are the ones provided by AI_BASE_TO_INT.

	@see BaseCodeToIntCode
	@see SNucleotideCodeToIntCode
	@see AI_BASE_TO_INT
	@see nucleic_acids::SNucleotideCode
	@see I_NUCLEOTIDES_COUNT

	@hideinitializer
	*/
	const nucleic_acids::SNucleotideCode AT_INT_TO_BASE[I_NUCLEOTIDES_COUNT] =
	{
		//            0x00,                   0x01,                   0x02,                  0x03,
		nucleic_acids::Gap, nucleic_acids::Guanine, nucleic_acids::Adenine, nucleic_acids::Purine,
		//                 0x04,                  0x05,                 0x06,                0x07,
		nucleic_acids::Cytosine, nucleic_acids::Strong, nucleic_acids::Amino, nucleic_acids::NotT,
		//                0x08,                0x09,                0x0A,                0x0B,
		nucleic_acids::Thymine, nucleic_acids::Keto, nucleic_acids::Weak, nucleic_acids::NotC,
		//                   0x0C,                0x0D,                0x0E,                0x0F
		nucleic_acids::Pyrimidine, nucleic_acids::NotA, nucleic_acids::NotG,  nucleic_acids::Any
	};


	//! Converts an integer code into its associated base code.
	const nucleic_acids::SNucleotideCode& IntCodeToSNucleotideCode(const int &iCode);


	/**
	@var AAB_BASES_MATCH
	********************
	@brief AAB_BASES_MATCH tells if two bases codes (with groups) match.

	This lookup table tells if the first base code given matches the second
	base code. If the first base or base group is equal to the second base
	or is part of the second base group, the table returns true. If the second
	base or base group is equal to the first base or is part of the first base
	group, the table returns true.

	@note index order does not import. For asymmetric matching, see
	@c AAB_BASES_MATCH_ASYMMETRIC.

	@see DoesBaseIntCodeMatchOther
	@see DoesBaseCodeMatchOther
	@see DoesSNucleotideCodeMatchOther
	@see AI_BASE_TO_INT
	@see AAB_BASES_MATCH_ASYMMETRIC

	@hideinitializer
	*/
	const bool AAB_BASES_MATCH[I_NUCLEOTIDES_COUNT][I_NUCLEOTIDES_COUNT] =
	{
		// second:
		//   -      G      A      R      C      S      M      V      T      K      W      D      Y      B      H      N
	// first:
	// -
		{ true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false},
	// G
		{ false,  true, false,  true, false,  true, false,  true, false,  true, false,  true, false,  true, false,  true},
	// A
		{ false, false,  true,  true, false, false,  true,  true, false, false,  true,  true, false, false,  true,  true},
	// R
		{ false,  true,  true,  true, false, false, false,  true, false, false, false,  true, false, false, false,  true},
	// C
		{ false, false, false, false,  true,  true,  true,  true, false, false, false, false,  true,  true,  true,  true},
	// S
		{ false,  true, false, false,  true,  true, false,  true, false, false, false, false, false,  true, false,  true},
	// M
		{ false, false,  true, false,  true, false,  true,  true, false, false, false, false, false, false,  true,  true},
	// V
		{ false,  true,  true,  true,  true,  true,  true,  true, false, false, false, false, false, false, false,  true},
	// T
		{ false, false, false, false, false, false, false, false,  true,  true,  true,  true,  true,  true,  true,  true},
	// K
		{ false,  true, false, false, false, false, false, false,  true,  true, false,  true, false,  true, false,  true},
	// W
		{ false, false,  true, false, false, false, false, false,  true, false,  true,  true, false, false,  true,  true},
	// D
		{ false,  true,  true,  true, false, false, false, false,  true,  true,  true,  true, false, false, false,  true},
	// Y
		{ false, false, false, false,  true, false, false, false,  true, false, false, false,  true,  true,  true,  true},
	// B
		{ false,  true, false, false,  true,  true, false, false,  true,  true, false, false,  true,  true, false,  true},
	// H
		{ false, false,  true, false,  true, false,  true, false,  true, false,  true, false,  true, false,  true,  true},
	// N
		{ false,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true},
	};


	/**
	@var AAB_BASES_MATCH_ASYMMETRIC
	*******************************
	@brief AAB_BASES_MATCH_ASYMMETRIC tells if two bases codes (with groups) match (with asymmetry).

	This lookup table tells if the first base code given matches the second
	base code. It means that all the bases included in the first base code
	must be also included in the second base code, however all the bases in
	the second base code may not be in the first base code. For example,
	both matching Adenine and Purine or matching Keto and Any will return
	true while matching Any and Guanine will return false.
	If one needs to exact match groups, the use of "==" operator would be
	more appropriate.

	@note index order is important as matching 'A' with 'N' is different from
	 matching 'N' and 'A' ("true" in the first case, "false" in the second)!
	 It means that 'A' is included into 'N' but 'N' is not included into 'A'.
	 For symmetric matching, see @c AAB_BASES_MATCH.

	@see DoesBaseIntCodeMatchOtherA
	@see DoesBaseCodeMatchOtherA
	@see DoesSNucleotideCodeMatchOtherA
	@see AI_BASE_TO_INT
	@see AAB_BASES_MATCH

	@hideinitializer
	*/
	const bool AAB_BASES_MATCH_ASYMMETRIC[I_NUCLEOTIDES_COUNT][I_NUCLEOTIDES_COUNT] =
	{
		// second:
		//   -      G      A      R      C      S      M      V      T      K      W      D      Y      B      H      N
	// first:
	// -
		{ true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false},
	// G
		{false,  true, false,  true, false,  true, false,  true, false,  true, false,  true, false,  true, false,  true},
	// A
		{false, false,  true,  true, false, false,  true,  true, false, false,  true,  true, false, false,  true,  true},
	// R
		{false, false, false,  true, false, false, false,  true, false, false, false,  true, false, false, false,  true},
	// C
		{false, false, false, false,  true,  true,  true,  true, false, false, false, false,  true,  true,  true,  true},
	// S
		{false, false, false, false, false,  true, false,  true, false, false, false, false, false,  true, false,  true},
	// M
		{false, false, false, false, false, false,  true,  true, false, false, false, false, false, false,  true,  true},
	// V
		{false, false, false, false, false, false, false,  true, false, false, false, false, false, false, false,  true},
	// T
		{false, false, false, false, false, false, false, false,  true,  true,  true,  true,  true,  true,  true,  true},
	// K
		{false, false, false, false, false, false, false, false, false,  true, false,  true, false,  true, false,  true},
	// W
		{false, false, false, false, false, false, false, false, false, false,  true,  true, false, false,  true,  true},
	// D
		{false, false, false, false, false, false, false, false, false, false, false,  true, false, false, false,  true},
	// Y
		{false, false, false, false, false, false, false, false, false, false, false, false,  true,  true,  true,  true},
	// B
		{false, false, false, false, false, false, false, false, false, false, false, false, false,  true, false,  true},
	// H
		{false, false, false, false, false, false, false, false, false, false, false, false, false, false,  true,  true},
	// N
		{false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,  true},
	};


	//! Tells if the 1st base/group code matches the 2nd.
	const bool& DoesBaseIntCodeMatchOther(const int &iFirstCode, const int& iSecondCode);


	//! Tells if the 1st base/group code matches the 2nd.
	const bool& DoesBaseCodeMatchOther(const char &cFirstCode, const char& cSecondCode);


	//! Tells if the 1st base/group code matches the 2nd.
	const bool& DoesSNucleotideCodeMatchOther(const nucleic_acids::SNucleotideCode &tFirstCode, const nucleic_acids::SNucleotideCode& tSecondCode);


	//! Tells if the 1st base/group code matches the 2nd or is included in its group.
	const bool& DoesBaseIntCodeMatchOtherA(const int &iFirstCode, const int& iSecondCode);


	//! Tells if the 1st base/group code matches the 2nd or is included in its group.
	const bool& DoesBaseCodeMatchOtherA(const char &cFirstCode, const char& cSecondCode);


	//! Tells if the 1st base/group code matches the 2nd or is included in its group.
	const bool& DoesSNucleotideCodeMatchOtherA(const nucleic_acids::SNucleotideCode &tFirstCode, const nucleic_acids::SNucleotideCode& tSecondCode);


	/**
	@var AI_BASE_TO_COMPLEMENT
	**************************
	@brief AI_BASE_TO_COMPLEMENT returns the complement of the given base.

	This lookup table can convert very quickly and efficiently a base into
	its complement. The table index is the integer code of the base to
	complement. Bases are complemented using the following rules:
	 - gap=>gap
	 - G=>C
	 - A=>T
	 - R=>Y
	 - C=>G
	 - S=>S
	 - M=>K
	 - V=>B
	 - T=>A
	 - U=>A
	 - K=>M
	 - W=>W
	 - D=>H
	 - Y=>R
	 - B=>V
	 - H=>D
	 - N=>N

	Example:
	In this example, we check that 'T' (integer code 0x08) is the complement
	of 'A' (integer code 0x02).
	@code assert(0x08 == AI_BASE_TO_COMPLEMENT[0x02]); @endcode

	@see AT_BASE_TO_COMPLEMENT

	@hideinitializer
	*/
	const int AI_BASE_TO_COMPLEMENT[I_NUCLEOTIDES_COUNT]=
	{
	//  -=>-  G=>C  A=>T  R=>Y  C=>G  S=>S  M=>K  V=>B
		0x00, 0x04, 0x08, 0x0C, 0x01, 0x05, 0x09, 0x0D,
	//  T=>A  K=>M  W=>W  D=>H  Y=>R  B=>V  H=>D  N=>N
		0x02, 0x06, 0x0A, 0x0E, 0x03, 0x07, 0x0B, 0x0F
	};


	//! Returns the complement of the given integer code.
	const int& GetComplementBaseIntCode(const int &iCode);


	/**
	@var AT_BASE_TO_COMPLEMENT
	**************************
	@brief AT_BASE_TO_COMPLEMENT returns the complement of the given base

	This lookup table can convert very quickly and efficiently a base code into
	its complement. The table index is either a character/base code or the
	@c SNucleotideCode of the base to complement. The table is case unsensitive
	and returns an @c SNucleotideCode (capital letters). Bases are complemented
	using the following rules:
	 - gap=>gap
	 - G=>C
	 - A=>T
	 - R=>Y
	 - C=>G
	 - S=>S
	 - M=>K
	 - V=>B
	 - T=>A
	 - U=>A
	 - K=>M
	 - W=>W
	 - D=>H
	 - Y=>R
	 - B=>V
	 - H=>D
	 - N=>N

	Example:
	In this example, we check that Thymine is the complement of Adenine.
	@code
	assert(nucleic_acids::Thymine == AT_BASE_TO_COMPLEMENT[nucleic_acids::Adenine]); // ok
	@endcode

	@see AI_BASE_TO_COMPLEMENT
	@see SNucleotideCodeToTableIndex

	@hideinitializer
	*/
	const nucleic_acids::SNucleotideCode AT_BASE_TO_COMPLEMENT[I_ASCII_CHAR_COUNT]=
	{
		nucleic_acids::Null, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // 00h-0Fh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // 10h-1Fh
		//               Space                       !                       "                       #                       $                       %                       &                       '                       (                       )                       *                       +                       ,                   -                       .                       /
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Gap, nucleic_acids::Unknown, nucleic_acids::Unknown, // 20h-2Fh
		//                   0                       1                       2                       3                       4                       5                       6                       7                       8                       9                       :                       ;                       <                       =                       >                       ?
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // 30h-3Fh
		//                   @                       A                    B                       C                    D                       E                       F                        G                    H                       I                       J                     K                       L                    M                   N                       O
		nucleic_acids::Unknown, nucleic_acids::Thymine, nucleic_acids::NotT, nucleic_acids::Guanine, nucleic_acids::NotG, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Cytosine, nucleic_acids::NotC, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Amino, nucleic_acids::Unknown, nucleic_acids::Keto, nucleic_acids::Any, nucleic_acids::Unknown, // 40h-4Fh
		//                   P                       Q                          R                      S                       T                       U                    V                    W                       X                      Y                       Z                       [                       \                       ]                       ^                       _
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Pyrimidine, nucleic_acids::Strong, nucleic_acids::Adenine, nucleic_acids::Adenine, nucleic_acids::NotA, nucleic_acids::Weak, nucleic_acids::Unknown, nucleic_acids::Purine, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // 50h-5Fh
		//                   `                       a                    b                       c                    d                       e                       f                        g                    h                       i                       j                     k                       l                    m                   n                       o
		nucleic_acids::Unknown, nucleic_acids::Thymine, nucleic_acids::NotT, nucleic_acids::Guanine, nucleic_acids::NotG, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Cytosine, nucleic_acids::NotC, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Amino, nucleic_acids::Unknown, nucleic_acids::Keto, nucleic_acids::Any, nucleic_acids::Unknown, // 60h-6Fh
		//                   p                       q                          r                      s                       t                       u                    v                    w                       x                      y                       z                       {                       |                       }                       ~
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Pyrimidine, nucleic_acids::Strong, nucleic_acids::Adenine, nucleic_acids::Adenine, nucleic_acids::NotA, nucleic_acids::Weak, nucleic_acids::Unknown, nucleic_acids::Purine, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // 70h-7Fh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // 80h-8Fh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // 90h-9Fh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // A0h-AFh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // B0h-BFh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // C0h-CFh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // D0h-DFh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // E0h-EFh
		nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, nucleic_acids::Unknown, // F0h-FFh
	};


	//! Returns the complement code of the given base code.
	const nucleic_acids::SNucleotideCode& GetComplementSNucleotideCode(const char &cBaseCode);


	//! Returns the complement code of the given nucleotide code.
	const nucleic_acids::SNucleotideCode& GetComplementSNucleotideCode(const nucleic_acids::SNucleotideCode& tNucleotideCode);


    // placeholder


	//! Returns the minimum value between two floats
	const float& GetMin(const float& fFirstValue, const float& fSecondValue);

	//! Returns the maximum value between two floats
	const float& GetMax(const float& fFirstValue, const float& fSecondValue);

}; // namespace biocpp
}; // namespace val
#endif //ifndef __BIOCPP_H
