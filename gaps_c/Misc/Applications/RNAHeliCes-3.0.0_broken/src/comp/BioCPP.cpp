//*****************************************************************************
#include <string>
#include <stdexcept>
#include "BioCPP.hpp"
#include "BioExceptions.hpp"
#include "DebugTools.hpp"

namespace val
{
namespace biocpp
{
	/**
	@fn nucleic_acids::SNucleotideCode& CharToSNucleotideCode(char &cBaseCode)
	**************************************************************************
	Converts a single character representing a nucleotide one-letter-code into
	its associated base code as a structured type. @c SNucleotideCode type is
	safer to use than the char type. The way the convertion is done almost
	brings no cost but the converted letter-code will be wrong if the given
	@c cBaseCode wasn't a valid nucleotide base code so the convertion can
	bring some unwanted values if it is not used safely. To use a safer
	converting tool, see @c CheckAndUpcaseBaseCodes.

	Example:
	@code
	char cBaseCode = 'a'; // lower case
	nucleic_acids::SNucleotideCode tMyBase = CharToSNucleotideCode(cBaseCode);
	// warning: here, (nucleic_acids::Adenine != tMyBase) because
	// "nucleic_acids::Adenine" is an upper case 'A' and
	// not a lower case 'a' like tMyBase.
	std::cout << tMyBase.m_cCode << std::endl; // displays "a"
	@endcode

	@param char &cBaseCode:
	 the character to convert.

	@return nucleic_acids::SNucleotideCode&:
	 a reference of type @c SNucleotideCode on the same memory used by the
	 given char @c cBaseCode.

	@see nucleic_acids::SNucleotideCode
	@see CharArrayToSNucleotideCodeArray
	@see AT_CHAR_TO_SNUC
	@see CheckBaseCodes
	@see CheckAndUpcaseBaseCodes
	@see SNucleotideCodeToChar
	@see SNucleotideCodeArrayToCharArray
	*/
	nucleic_acids::SNucleotideCode& CharToSNucleotideCode(char &cBaseCode)
	{
		return (*(nucleic_acids::SNucleotideCode*) &cBaseCode);
	}


	/**
	@fn nucleic_acids::SNucleotideCode* CharArrayToSNucleotideCodeArray(char *acBaseCodes)
	**************************************************************************************
	Does the same job as @c CharToSNucleotideCode but for arrays. It has the same
	safety problems as @c CharToSNucleotideCode.
	
	@param char *acBaseCodes:
	 the array to convert.

	@return nucleic_acids::SNucleotideCode*:
	 a @c SNucleotideCode* type pointer to the same array @c acBaseCodes.

	@see nucleic_acids::SNucleotideCode
	@see CharToSNucleotideCode
	@see CharToSNucleotideCodeSafe
	@see AT_CHAR_TO_SNUC
	@see CheckBaseCodes
	@see CheckAndUpcaseBaseCodes
	@see SNucleotideCodeToChar
	@see SNucleotideCodeArrayToCharArray
	*/
	nucleic_acids::SNucleotideCode* CharArrayToSNucleotideCodeArray(char *acBaseCodes)
	{
		return (nucleic_acids::SNucleotideCode*)acBaseCodes;
	}


	/**
	@fn const nucleic_acids::SNucleotideCode& CharToSNucleotideCodeSafe(const char &cBaseCode)
	******************************************************************************************
	Converts a single character representing a base code into its associated
	nucleotide code. If the given base code is unknown, it will return
	@c nucleic_acids::Unknown. This version is a little bit slower but safer
	than @c CharToSNucleotide because it also converts the character case.
	However, the function uses a lookup table and is really fast.
	
	@note This function is just a wrapper for the table @c AT_CHAR_TO_SNUC.
	 Using the wrapper is safer than using the lookup table directly because
	 it checks the type of the given parameters and insures the char used for
	 indexing is unsigned.

	Example:
	@code
	nucleic_acids::SNucleotideCode tMyBase;
	tMyBase = CharToSNucleotideCodeSafe('a'); // lower case
	assert(nucleic_acids::Adenine == tMyBase); // ok
	std::cout << tMyBase.m_cBaseCode << std::endl; // displays "A", upper case
	@endcode

	@param const char &cBaseCode:
	 the base code to convert.

	@return const nucleic_acids::SNucleotideCode&:
	 the nucleotide code corresponding to the given character or
	 @c nucleic_acids::Unknown if the base code is unknown.

	@see AT_CHAR_TO_SNUC
	@see nucleic_acids::SNucleotideCode
	@see CharToSNucleotideCode
	@see CharArrayToSNucleotideCodeArray
	@see CheckBaseCodes
	@see CheckAndUpcaseBaseCodes
	@see SNucleotideCodeToChar
	@see SNucleotideCodeArrayToCharArray
	*/
	const nucleic_acids::SNucleotideCode& CharToSNucleotideCodeSafe(const char &cBaseCode)
	{
		return AT_CHAR_TO_SNUC[(unsigned)cBaseCode];
	}


	/**
	@fn void CheckBaseCodes(const char *acBaseArray, const int &iArraySize)
	***********************************************************************
	Checks the given array of char for invalid base codes and throws a
	@c bioexceptions::CInvalidNucleotide exception if an invalid base code is
	found. This is a case sensitive function (for example 'a' is an invalid
	code while 'A' (==nucleic_acids::Adenine) is valid).

	@note nucleotide codes check is made using AT_CHAR_TO_SNUC table.

	@param const char *acBaseArray:
	 the array to check.

	@param const int &iArraySize:
	 the array size or a negative value if it's a null-terminated array.
	 Default value: -1

	@throw bioexceptions::CInvalidNucleotide:
	 if an invalid nucleotide was met. Both the invalid nucleotide code and
	 position field of the exception object are filled.

	@see nucleic_acids::SNucleotideCode
	@see CharToSNucleotideCode
	@see CharToSNucleotideCodeSafe
	@see CharArrayToSNucleotideCodeArray
	@see AT_CHAR_TO_SNUC
	@see CheckAndUpcaseBaseCodes
	@see SNucleotideCodeToChar
	@see SNucleotideCodeArrayToCharArray
	*/
	void CheckBaseCodes(const char *acBaseArray, const int &iArraySize)
	{
		TRACED("void CheckBaseCodes(const char *acBaseArray, const int &iArraySize)");
		if (0 <= iArraySize)
		{ // array size was specified
			for (int i = 0; i < iArraySize; i++)
			{
				// check if the char has a non-unknown associated nucleotide
				// and if it's an upper case char
				if ((nucleic_acids::Unknown == AT_CHAR_TO_SNUC[(unsigned int)acBaseArray[i]])
					|| (('a' <= acBaseArray[i]) && (acBaseArray[i] <= 'z')))
				{
					throw bioexceptions::CInvalidNucleotide(STR_ERROR_INVALID_NUCLEOTIDE, acBaseArray[i], i);
				}
			}
		}
		else
		{ // null-terminated array
			int i = 0;
			while('\0' != acBaseArray[i])
			{
				// check if the char has a non-unknown associated nucleotide
				// and if it's an upper case char
				if ((nucleic_acids::Unknown == AT_CHAR_TO_SNUC[(unsigned int)acBaseArray[i]])
					|| (('a' <= acBaseArray[i]) && (acBaseArray[i] <= 'z')))
				{
					throw bioexceptions::CInvalidNucleotide(STR_ERROR_INVALID_NUCLEOTIDE, acBaseArray[i], i);
				}
				i++;
			}
		}
	} // void CheckBaseCodes


	/**
	@fn void CheckAndUpcaseBaseCodes(char *acBaseArray, const int &iArraySize)
	**************************************************************************
	Checks the given array of char for invalid nucleotides codes and throws a
	@c bioexceptions::CInvalidNucleotide exception if an invalid nucleotide is
	found. If a lower case character is found, it is changed to its uppercase
	version. This is a safe way to convert an array of char into an array of
	@c nucleic_acids::SNucleotideCode base codes.

	Example:
	@code
	 using namespace val::biocpp;
	 char acMySequence[9] = "cgtaATGC";
	 nucleic_acids::SNucleotideCode* atNucleotideSequence = ACHAR_TO_ASNUC(acMySequence);
	 // unsafe to use atNucleotideSequence here
	 try
	 {
	 	CheckAndUpcaseBaseCodes(acMySequence);
		// now safe to use atNucleotideSequence from here
		assert(nucleic_acids::Cytosine == atNucleotideSequence[0]); // ok here as atNucleotideSequence[0] == 'C'
	 }
	 catch (bioexceptions::CInvalidNucleotide tException)
	 {
	 	std::cout << "Invalid nucleotide code '" << tException.GetCode() << "' found at position " << tException.GetPosition()+1 << std::endl;
	 }
	@endcode

	@note nucleotide codes check is made using AT_CHAR_TO_SNUC table.

	@param char *acBaseArray:
	 the array to check and capitalize.

	@param const int &iArraySize:
	 the array size or a negative value if it's a null-terminated array.
	 Default value: -1

	@throw bioexceptions::CInvalidNucleotide:
	 if an invalid nucleotide was met. Both the invalid nucleotide code and
	 position field of the exception object are filled.

	@see nucleic_acids::SNucleotideCode
	@see CharToSNucleotideCode
	@see CharToSNucleotideCodeSafe
	@see CharArrayToSNucleotideCodeArray
	@see AT_CHAR_TO_SNUC
	@see CheckBaseCodes
	@see SNucleotideCodeToChar
	@see SNucleotideCodeArrayToCharArray
	*/
	void CheckAndUpcaseBaseCodes(char *acBaseArray, const int &iArraySize)
	{
		TRACED("void CheckAndUpcaseBaseCodes(char *acBaseArray, const int &iArraySize)");
		if (0 <= iArraySize)
		{ // array size was specified
			for (int i = 0; i < iArraySize; i++)
			{
				// check if the char has a non-unknown associated nucleotide
				// and if it's an upper case char
				if (nucleic_acids::Unknown == AT_CHAR_TO_SNUC[(unsigned int)acBaseArray[i]])
				{
					throw bioexceptions::CInvalidNucleotide(STR_ERROR_INVALID_NUCLEOTIDE, acBaseArray[i], i);
				}
				if (('a' <= acBaseArray[i]) && (acBaseArray[i] <= 'z'))
				{
					acBaseArray[i] += ('A' - 'a');
				}
			}
		}
		else
		{ // null-terminated array
			int i = 0;
			while('\0' != acBaseArray[i])
			{
				// check if the char has a non-unknown associated nucleotide
				// and if it's an upper case char
				if (nucleic_acids::Unknown == AT_CHAR_TO_SNUC[(unsigned int)acBaseArray[i]])
				{
					throw bioexceptions::CInvalidNucleotide(STR_ERROR_INVALID_NUCLEOTIDE, acBaseArray[i], i);
				}
				if (('a' <= acBaseArray[i]) && (acBaseArray[i] <= 'z'))
				{
					acBaseArray[i] += ('A' - 'a');
				}
				i++;
			}
		}
	} // void CheckAndUpcaseBaseCodes (char version)
	void CheckAndUpcaseBaseCodes(std::string strBaseArray, const int &iArraySize)
	{
		TRACED("void CheckAndUpcaseBaseCodes(std::string strBaseArray, const int &iArraySize)");
		int iRealArraySize = iArraySize;
		if (0 > iRealArraySize)
		{ // array size was specified
			iRealArraySize = (signed)strBaseArray.length();
		}
		for (int i = 0; i < iRealArraySize; i++)
		{
			// check if the char has a non-unknown associated nucleotide
			// and if it's an upper case char
			if (nucleic_acids::Unknown == AT_CHAR_TO_SNUC[(unsigned int)strBaseArray[i]])
			{
				throw bioexceptions::CInvalidNucleotide(STR_ERROR_INVALID_NUCLEOTIDE, strBaseArray[i], i);
			}
			if (('a' <= strBaseArray[i]) && (strBaseArray[i] <= 'z'))
			{
				strBaseArray[i] += ('A' - 'a');
			}
		}
	} // void CheckAndUpcaseBaseCodes (string version)


	/**
	@fn char& SNucleotideCodeToChar(nucleic_acids::SNucleotideCode& tNucleotideCode)
	****************************************************************************
	Converts a nucleotide code into its associated character representing its
	base code. Converted characters are usually upper case as their values
	initially come from one of the main defined bases constant. However, if
	custom nucleotides codes are used, the converted character may be lower
	case.
	
	@note @c SNucleotideCode type is safer to use than the char type.

	@param nucleic_acids::SNucleotideCode& tNucleotideCode:
	 the nucleotide code to convert (type cast).

	@return char&:
	 a char reference to the same memory space used by @c tNucleotideCode.

	@see nucleic_acids::SNucleotideCode
	@see CharToSNucleotideCode
	@see CharToSNucleotideCodeSafe
	@see CharArrayToSNucleotideCodeArray
	@see AT_CHAR_TO_SNUC
	@see CheckBaseCodes
	@see CheckAndUpcaseBaseCodes
	@see SNucleotideCodeArrayToCharArray
	*/
	char& SNucleotideCodeToChar(nucleic_acids::SNucleotideCode& tNucleotideCode)
	{
		return (*(char*)&tNucleotideCode);
	}


	/**
	@fn const char& SNucleotideCodeToChar(const nucleic_acids::SNucleotideCode& tNucleotideCode)
	****************************************************************************************
	Converts a nucleotide code into its associated character representing its
	base code. Converted characters are usually upper case as their values
	initially come from one of the main defined bases constant. However, if
	custom nucleotides codes are used, the converted character may be lower
	case.
	
	@note @c SNucleotideCode type is safer to use than the char type.

	@param const nucleic_acids::SNucleotideCode& tNucleotideCode:
	 the nucleotide code to convert (type cast).

	@return const char&:
	 a char reference to the same memory space used by @c tNucleotideCode.

	@see nucleic_acids::SNucleotideCode
	@see CharToSNucleotideCode
	@see CharToSNucleotideCodeSafe
	@see CharArrayToSNucleotideCodeArray
	@see AT_CHAR_TO_SNUC
	@see CheckBaseCodes
	@see CheckAndUpcaseBaseCodes
	@see SNucleotideCodeArrayToCharArray
	*/
	const char& SNucleotideCodeToChar(const nucleic_acids::SNucleotideCode &tNucleotideCode)
	{
		return (*(const char*)&tNucleotideCode);
	}


	/**
	@fn char* SNucleotideCodeArrayToCharArray(nucleic_acids::SNucleotideCode* atNucleotideCode)
	***************************************************************************************
	Does the same job as @c SNucleotideCodeToChar but for arrays.
	
	@param nucleic_acids::SNucleotideCode* atNucleotideCode:
	 the array to convert.

	@return char*:
	 a char pointer to the same array @c atNucleotideCode.

	@see nucleic_acids::SNucleotideCode
	@see CharToSNucleotideCode
	@see CharToSNucleotideCodeSafe
	@see CharArrayToSNucleotideCodeArray
	@see AT_CHAR_TO_SNUC
	@see CheckBaseCodes
	@see CheckAndUpcaseBaseCodes
	@see SNucleotideCodeToChar
	*/
	char* SNucleotideCodeArrayToCharArray(nucleic_acids::SNucleotideCode* atNucleotideCode)
	{
		return (char*)atNucleotideCode;
	}


	/**
	@fn const char* SNucleotideCodeArrayToCharArray(const nucleic_acids::SNucleotideCode* atNucleotideCode)
	*******************************************************************************************************
	Does the same job as @c SNucleotideCodeToChar but for arrays.
	
	@param const nucleic_acids::SNucleotideCode* atNucleotideCode:
	 the array to convert.

	@return const char*:
	 a char pointer to the same array @c atNucleotideCode.

	@see nucleic_acids::SNucleotideCode
	@see CharToSNucleotideCode
	@see CharToSNucleotideCodeSafe
	@see CharArrayToSNucleotideCodeArray
	@see AT_CHAR_TO_SNUC
	@see CheckBaseCodes
	@see CheckAndUpcaseBaseCodes
	@see SNucleotideCodeToChar
	*/
	const char* SNucleotideCodeArrayToCharArray(const nucleic_acids::SNucleotideCode* atNucleotideCode)
	{
		return (const char*)atNucleotideCode;
	}


	/**
	@fn const unsigned char& SNucleotideCodeToTableIndex(const nucleic_acids::SNucleotideCode& tNucleotideCode)
	*******************************************************************************************************
	Converts a single nucleotide code into a [0 to 255] index for lookup
	tables that contain @c I_ASCII_CHAR_COUNT elements. The index value is the
	ASCII code of the character representing the nucleotide base code.

	@param const nucleic_acids::SNucleotideCode& tNucleotideCode:
	 the nucleotide code to convert.

	@return const unsigned char&:
	 an index for lookup tables of size @c I_ASCII_CHAR_COUNT.

	@see nucleic_acids::SNucleotideCode
	*/
	const unsigned char& SNucleotideCodeToTableIndex(const nucleic_acids::SNucleotideCode& tNucleotideCode)
	{
		return (*(const unsigned char*)&tNucleotideCode);
	}


	/*
	@fn const int& BaseCodeToIntCode(const char &cBaseCode)
	*******************************************************
	Converts a base code into its associated integer code. This function is
	just a wrapper for the lookup table @c AI_BASE_TO_INT and is really fast.
	It makes sure the table is properly indexed (by un unsigned char).
	To have more details about how codes are associated, please refer to
	@c AI_BASE_TO_INT.

	@param const char &cBaseCode:
	 the base code to convert.

	@return const int&:
	 the integer code associated to the given base code.

	@see IntCodeToSNucleotideCode
	@see SNucleotideCodeToIntCode
	@see AI_BASE_TO_INT
	@see AT_INT_TO_BASE
	@see nucleic_acids::SNucleotideCode
	*/
	const int& BaseCodeToIntCode(const char &cBaseCode)
	{
		return AI_BASE_TO_INT[(unsigned)cBaseCode];
	}


	/*
	@fn const int& SNucleotideCodeToIntCode(const nucleic_acids::SNucleotideCode &tNucleotideCode)
	**********************************************************************************************
	Converts a nucleotide code into its associated integer code. This function
	is just a wrapper for the lookup table @c AI_BASE_TO_INT and is really
	fast. It makes sure the table is properly indexed. To have more details
	about how codes are associated, please refer to @c AI_BASE_TO_INT.

	@param const nucleic_acids::SNucleotideCode &tNucleotideCode:
	 the nucleotide code to convert.

	@return const int&:
	 the integer code associated to the given nucleotide code.

	@see IntCodeToSNucleotideCode
	@see AI_BASE_TO_INT
	@see AT_INT_TO_BASE
	@see nucleic_acids::SNucleotideCode
	*/
	const int& SNucleotideCodeToIntCode(const nucleic_acids::SNucleotideCode &tNucleotideCode)
	{
		return AI_BASE_TO_INT[SNucleotideCodeToTableIndex(tNucleotideCode)];
	}


	/**
	@fn const nucleic_acids::SNucleotideCode& IntCodeToSNucleotideCode(const int &iCode)
	************************************************************************************
	This function is just a wrapper for @c AT_INT_TO_BASE lookup table. It does
	a parameter check and is safer to use than the lookup table directly even
	if it is a very little bit slower.

	@note as T and U have the same integer value, only T is returned, never U.
	@note integer codes used are the ones provided by @c AI_BASE_TO_INT.

	@param const int &iCode:
	 the integer code associated to a nucleotide code.

	@return const nucleic_acids::SNucleotideCode&:
	 the associated nucleotide code for the given integer code.

	@throw std::range_error:
	 if the given integer code is not a valid code associated to a nucleotide.
	 It happens if @c iCode is negative or above or equal to
	 @c I_NUCLEOTIDES_COUNT.

	@see SNucleotideCodeToIntCode
	@see BaseCodeToIntCode
	@see AT_INT_TO_BASE
	@see AI_BASE_TO_INT
	@see nucleic_acids::SNucleotideCode
	*/
	const nucleic_acids::SNucleotideCode& IntCodeToSNucleotideCode(const int &iCode)
	{
		if ((0 > iCode) || (I_NUCLEOTIDES_COUNT <= iCode))
		{
			throw std::range_error("Integer code given to IntCodeToSNucleotideCode is not in the correct range!");
		}
		return AT_INT_TO_BASE[iCode];
	}


	/*
	@fn const bool& DoesBaseIntCodeMatchOther(const int &iFirstCode, const int& iSecondCode)
	****************************************************************************************
	Tells if the 1st integer code matches the 2nd or is included in its group.
	This function is just a wrapper for @c AAB_BASES_MATCH lookup table. It
	insures that the table is correctly indexed. See @c AAB_BASES_MATCH for
	matching details.

	@note index order is not important. See @c DoesBaseIntCodeMatchOtherA for
	asymmetric matching.
	
	@param const int &iFirstCode:
	 the first integer code (or group code).

	@param const int& iSecondCode:
	 the second integer code (or group code).

	@return const bool&:
	 true if the nucleotide or the group of nucleotides represented by the
	 first integer matches the nucleotide or the group of nucleotides
	 represented by the second integer. Returns false otherwise.

	@throw std::range_error:
	 if one of the given integer codes is not a valid code associated to a
	 nucleotide. It happens if an integer code is negative or above or equal to
	 @c I_NUCLEOTIDES_COUNT.

	@see AAB_BASES_MATCH
	@see I_NUCLEOTIDES_COUNT
	@see DoesBaseIntCodeMatchOtherA
	*/
	const bool& DoesBaseIntCodeMatchOther(const int &iFirstCode, const int& iSecondCode)
	{
		if ((0 > iFirstCode) || (I_NUCLEOTIDES_COUNT <= iFirstCode))
		{
			throw std::range_error("First integer code given to DoesBaseIntCodeMatchOther is not in the correct range!");
		}
		if ((0 > iSecondCode) || (I_NUCLEOTIDES_COUNT <= iSecondCode))
		{
			throw std::range_error("Second integer code given to DoesBaseIntCodeMatchOther is not in the correct range!");
		}
		return AAB_BASES_MATCH[iFirstCode][iSecondCode];
	}


	/*
	@fn const bool& DoesBaseCodeMatchOther(const char &cFirstCode, const char& cSecondCode)
	***************************************************************************************
	Tells if the 1st base code matches the 2nd or is included in its group.
	This function is just a wrapper for @c AAB_BASES_MATCH lookup table. It
	insures that the table is correctly indexed. See @c AAB_BASES_MATCH for
	matching details.

	@note index order is not important. See @c DoesBaseCodeMatchOtherA for
	asymmetric matching.
	
	@param const char &cFirstCode:
	 the first base code (or group code).

	@param const char& cSecondCode:
	 the second base code (or group code).

	@return const bool&:
	 true if the nucleotide or the group of nucleotides represented by the
	 first base code matches the nucleotide or the group of nucleotides
	 represented by the second base code. Returns false otherwise.

	@see AAB_BASES_MATCH
	@see DoesBaseCodeMatchOtherA
	*/
	const bool& DoesBaseCodeMatchOther(const char &cFirstCode, const char& cSecondCode)
	{
		return AAB_BASES_MATCH[BaseCodeToIntCode(cFirstCode)][BaseCodeToIntCode(cSecondCode)];
	}


	/*
	@fn const bool& DoesSNucleotideCodeMatchOther(const nucleic_acids::SNucleotideCode &tFirstCode, const nucleic_acids::SNucleotideCode &tSecondCode)
	**************************************************************************************************************************************************
	Tells if the 1st nucleotide/group code matches the 2nd. This function is
	just a wrapper for @c AAB_BASES_MATCH lookup table.	It insures that the
	table is correctly indexed. See @c AAB_BASES_MATCH for matching details.

	@note index order is not important. See @c DoesSNucleotideCodeMatchOther
	for asymmetric matching.
	
	@param const nucleic_acids::SNucleotideCode &tFirstCode:
	 the first nucleotide code (or group code).

	@param const nucleic_acids::SNucleotideCode &tSecondCode:
	 the second nucleotide code (or group code).

	@return const bool&:
	 true if the nucleotide or the group of nucleotides code is matches the
	 nucleotide or the group of nucleotides of second nucleotide code.
	 Returns false otherwise.
	
	@see AAB_BASES_MATCH
	@see DoesSNucleotideCodeMatchOtherA
	*/
	const bool& DoesSNucleotideCodeMatchOther(const nucleic_acids::SNucleotideCode &tFirstCode, const nucleic_acids::SNucleotideCode &tSecondCode)
	{
		return AAB_BASES_MATCH[SNucleotideCodeToIntCode(tFirstCode)][SNucleotideCodeToIntCode(tSecondCode)];
	}


	/*
	@fn const bool& DoesBaseIntCodeMatchOtherA(const int &iFirstCode, const int& iSecondCode)
	*****************************************************************************************
	Tells if the 1st integer code matches the 2nd or is included in its group.
	This function is just a wrapper for @c AAB_BASES_MATCH_ASYMMETRIC lookup
	table. It insures that the table is correctly indexed. See
	@c AAB_BASES_MATCH_ASYMMETRIC for matching details.

	@note index order is important as matching 'A' with 'N' is different from
	 matching 'N' and 'A' ("true" in the first case, "false" in the second)!
	 It means that 'A' is included into 'N' but 'N' is not included into 'A'.
	 See @c DoesBaseIntCodeMatchOther for symmetric matching.
	
	@param const int &iFirstCode:
	 the first integer code (or group code).

	@param const int& iSecondCode:
	 the second integer code (or group code).

	@return const bool&:
	 true if the nucleotide or the group of nucleotides represented by the
	 first integer is equal or included in the nucleotide or the group of
	 nucleotides represented by the second integer. Returns false otherwise
	 (false if at least one base code included in the first group is not
	 included in the second group).

	@throw std::range_error:
	 if one of the given integer codes is not a valid code associated to a
	 nucleotide. It happens if an integer code is negative or above or equal to
	 @c I_NUCLEOTIDES_COUNT.

	@see AAB_BASES_MATCH_ASYMMETRIC
	@see I_NUCLEOTIDES_COUNT
	@see DoesBaseIntCodeMatchOther
	*/
	const bool& DoesBaseIntCodeMatchOtherA(const int &iFirstCode, const int& iSecondCode)
	{
		if ((0 > iFirstCode) || (I_NUCLEOTIDES_COUNT <= iFirstCode))
		{
			throw std::range_error("First integer code given to DoesBaseIntCodeMatchOtherA is not in the correct range!");
		}
		if ((0 > iSecondCode) || (I_NUCLEOTIDES_COUNT <= iSecondCode))
		{
			throw std::range_error("Second integer code given to DoesBaseIntCodeMatchOtherA is not in the correct range!");
		}
		return AAB_BASES_MATCH_ASYMMETRIC[iFirstCode][iSecondCode];
	}


	/*
	@fn const bool& DoesBaseCodeMatchOtherA(const char &cFirstCode, const char& cSecondCode)
	****************************************************************************************
	Tells if the 1st base code matches the 2nd or is included in its group.
	This function is just a wrapper for @c AAB_BASES_MATCH_ASYMMETRIC lookup
	table. It insures that the table is correctly indexed. See
	@c AAB_BASES_MATCH_ASYMMETRIC for matching details.

	@note index order is important as matching 'A' with 'N' is different from
	 matching 'N' and 'A' ("true" in the first case, "false" in the second)!
	 It means that 'A' is included into 'N' but 'N' is not included into 'A'.
	 See @c DoesBaseCodeMatchOther for symmetric matching.
	
	@param const char &cFirstCode:
	 the first base code (or group code).

	@param const char& cSecondCode:
	 the second base code (or group code).

	@return const bool&:
	 true if the nucleotide or the group of nucleotides represented by the
	 first base code is equal or included in the nucleotide or the group of
	 nucleotides represented by the second base code. Returns false otherwise
	 (false if at least one base code included in the first group is not
	 included in the second group).

	@see AAB_BASES_MATCH_ASYMMETRIC
	@see DoesBaseCodeMatchOther
	*/
	const bool& DoesBaseCodeMatchOtherA(const char &cFirstCode, const char& cSecondCode)
	{
		return AAB_BASES_MATCH_ASYMMETRIC[BaseCodeToIntCode(cFirstCode)][BaseCodeToIntCode(cSecondCode)];
	}


	/*
	@fn const bool& DoesSNucleotideCodeMatchOtherA(const nucleic_acids::SNucleotideCode &tFirstCode, const nucleic_acids::SNucleotideCode &tSecondCode)
	***************************************************************************************************************************************************
	Tells if the 1st nucleotide code matches the 2nd or is included in its
	group. This function is just a wrapper for @c AAB_BASES_MATCH_ASYMMETRIC
	lookup table. It insures that the table is correctly indexed. See
	@c AAB_BASES_MATCH_ASYMMETRIC for matching details.

	@note index order is important as matching 'A' with 'N' is different from
	 matching 'N' and 'A' ("true" in the first case, "false" in the second)!
	 It means that 'A' is included into 'N' but 'N' is not included into 'A'.
	 See @c DoesSNucleotideCodeMatchOther for symmetric matching.
	
	@param const nucleic_acids::SNucleotideCode &tFirstCode:
	 the first nucleotide code (or group code).

	@param const nucleic_acids::SNucleotideCode &tSecondCode:
	 the second nucleotide code (or group code).

	@return const bool&:
	 true if the nucleotide or the group of nucleotides code is equal or
	 included in the nucleotide or the group of nucleotides of second
	 nucleotide code. Returns false otherwise (false if at least one nucleotide
	 code included in the first group is not included in the second group).
	
	@see AAB_BASES_MATCH_ASYMMETRIC
	@see DoesSNucleotideCodeMatchOther
	*/
	const bool& DoesSNucleotideCodeMatchOtherA(const nucleic_acids::SNucleotideCode &tFirstCode, const nucleic_acids::SNucleotideCode &tSecondCode)
	{
		return AAB_BASES_MATCH_ASYMMETRIC[SNucleotideCodeToIntCode(tFirstCode)][SNucleotideCodeToIntCode(tSecondCode)];
	}


	/**
	@fn const int& GetComplementBaseIntCode(const int &iCode)
	*********************************************************
	This function is just a wrapper for @c AI_BASE_TO_COMPLEMENT lookup table.
	It returns the complement (integer code) of the given integer code. This
	function is safer to use than the lookup table directly because it checks
	the parameter. Please refer to @c AI_BASE_TO_COMPLEMENT for
	complementation details.

	@param const int &iCode:
	 the integer code of the nucleotide (or group of nucleotides) to
	 complement.

	@return const int&:
	 the integer code of the complement nucleotide (or group of nucleotides).

	@throw std::range_error:
	 if the given integer code @c iCode is not a valid code associated to a
	 nucleotide. It happens if @c iCode is negative or above or equal to
	 @c I_NUCLEOTIDES_COUNT.

	@see AI_BASE_TO_COMPLEMENT
	*/
	const int& GetComplementBaseIntCode(const int &iCode)
	{
		if ((0 > iCode) || (I_NUCLEOTIDES_COUNT <= iCode))
		{
			throw std::range_error("Integer code given to GetComplement is not in the correct range!");
		}
		return AI_BASE_TO_COMPLEMENT[iCode];
	}


	/**
	@fn const nucleic_acids::SNucleotideCode& GetComplementSNucleotideCode(const char &cBaseCode)
	*********************************************************************************************
	This function is just a wrapper for @c AT_BASE_TO_COMPLEMENT lookup table.
	It returns the complement (nucleotide code) of the given base code. This
	function is safer to use than the lookup table directly because it checks
	the parameter. Please refer to @c AT_BASE_TO_COMPLEMENT for
	complementation details.

	@note: this function returns Thymine as the complement of Adenine. This
	 can be an issue in case of RNA as the user would expect Uridine. In that
	 case, the user should work with GetComplementBaseIntCode function instead.

	@param const char &cBaseCode:
	 the base code of the nucleotide (or group of nucleotides) to
	 complement.

	@return const nucleic_acids::SNucleotideCode&:
	 the nucleotide code of the complement nucleotide (or group of nucleotides).

	@see AT_BASE_TO_COMPLEMENT
	*/
	const nucleic_acids::SNucleotideCode& GetComplementSNucleotideCode(const char &cBaseCode)
	{
		return AT_BASE_TO_COMPLEMENT[(unsigned)cBaseCode];
	}


	/**
	@fn const nucleic_acids::SNucleotideCode& GetComplementSNucleotideCode(const nucleic_acids::SNucleotideCode& tNucleotideCode)
	*****************************************************************************************************************************
	This function is just a wrapper for @c AT_BASE_TO_COMPLEMENT lookup table.
	It returns the complement (nucleotide code) of the given nucleotide code.
	This function is safer to use than the lookup table directly because it
	checks the parameter. Please refer to @c AT_BASE_TO_COMPLEMENT for
	complementation details.

	@note: this function returns Thymine as the complement of Adenine. This
	 can be an issue in case of RNA as the user would expect Uridine. In that
	 case, the user should work with GetComplementBaseIntCode function instead.

	@param const nucleic_acids::SNucleotideCode& tNucleotideCode:
	 the nucleotide code of the nucleotide (or group of nucleotides) to
	 complement.

	@return const nucleic_acids::SNucleotideCode&:
	 the nucleotide code of the complement nucleotide (or group of nucleotides).

	@see AT_BASE_TO_COMPLEMENT
	*/
	const nucleic_acids::SNucleotideCode& GetComplementSNucleotideCode(const nucleic_acids::SNucleotideCode& tNucleotideCode)
	{
		return AT_BASE_TO_COMPLEMENT[SNucleotideCodeToTableIndex(tNucleotideCode)];
	}


        // placeholder 


	/**
	@fn const float& GetMin(const float& fVal1, const float& fVal2)
	***************************************************************
	Returns the minimum value between the 2 given floats.
	
	@param const float& fFirstValue:
	 the first float value.
	
	@param const float& fSecondValue:
	 the second float value.

	@return const float&:
	 fFirstValue if fFirstValue is less than fSecondValue, fSecondValue otherwise.

	*/
	const float& GetMin(const float& fFirstValue, const float& fSecondValue)
	{
		return (fFirstValue < fSecondValue)? fFirstValue:fSecondValue;
	}


	/**
	@fn const float& GetMax(const float& fVal1, const float& fVal2)
	***************************************************************
	Returns the maximum value between the 2 given floats.
	
	@param const float& fFirstValue:
	 the first float value.
	
	@param const float& fSecondValue:
	 the second float value.

	@return const float&:
	 fFirstValue if fFirstValue is greater than fSecondValue, fSecondValue otherwise.

	*/
	const float& GetMax(const float& fFirstValue, const float& fSecondValue)
	{
		return (fFirstValue > fSecondValue)? fFirstValue:fSecondValue;
	}


}; // namespace biocpp
}; // namespace val
