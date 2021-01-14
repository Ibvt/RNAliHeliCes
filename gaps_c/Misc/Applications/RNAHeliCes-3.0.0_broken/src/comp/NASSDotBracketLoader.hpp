//*****************************************************************************
#ifndef __NASSDOTBRACKETLOADER_H
#define __NASSDOTBRACKETLOADER_H 1


/**
@file NASSDotBracketLoader.h
*****************************
@brief NASSDotBracketLoader.h is the header file for
@c val::biocpp::loaders::SNASSDotBracketLoader functor.

NASSDotBracketLoader.h is the header for
@c val::biocpp::loaders::SNASSDotBracketLoader. @c SNASSDotBracketLoader
loads a RNA secondary structure from a string which contains the secondary
structure in dot-bracket notation.

@see val::biocpp::loaders
@see val::biocpp::loaders::SLoader

@author Valentin GUIGNON
@version 1.0
@date 06/09/2004
*/
#include "Loader.hpp"
#include "NASecondaryStructure.hpp"
#include "BioCPP.hpp"
#include "DebugTools.hpp"
#include <stdexcept>
#include <string>
#include <sstream>
#include <stack>


namespace val
{
namespace biocpp
{
namespace loaders
{
/**
@struct SNASSDotBracketLoader   NASSDotBracketLoader.h "NASSDotBracketLoader.hpp"
*****************************
@brief Loads a secondary structure from a dot-bracket format string.

@c SNASSDotBracketLoader loads a RNA secondary structure from a string
which contains the secondary structure in dot-bracket notation.
*/
struct SNASSDotBracketLoader: public SLoader<CNASecondaryStructure, std::string>
{
	/**
	@fn CNASecondaryStructure operator()(const std::string &strDotBracketSequence)
	******************************************************************************
	@brief Loads a secondary structure from a string using dot-bracket notation.

	Loads a RNA secondary structure from a string which contains two ti three
	lines:
	the first line is optional and can contain a name in FASTA format
	the next line is the base codes sequence and the last line is the
	dot-bracket notation for that sequence.\n
	Example of string:
	@code
	>A Sequence
	ACCGCCCT
	(..(..))
	@endcode

	@param const std::string &strDotBracketSequence:\n
	 the string has two to three lines; the first one is optional and contains
	 a name for the sequence; the next line is the base codes sequence and the
	 last one is the dot-bracket annotation. The two last lines must have the
	 same length. The sequence line must only contain valid base codes; the
	 dot-bracket line must only contains dot-bracket and dots (no spaces!).

	@return CNASecondaryStructure:
	 a new instance of a nucleic acid secondary structure.

	@see CNASecondaryStructure
	*/
	CNASecondaryStructure operator()(const std::string &strDotBracketSequence)
	{
		std::stack<int> tIndexStack;
		CNASecondaryStructure tNewSecondaryStructure;

		std::istringstream tSequenceStream(strDotBracketSequence);
		std::string strLine;
		// skip empty lines or comments
		while (((0 == strLine.size()) || ('#' == strLine[0])) && (false == tSequenceStream.eof()))
		{
			std::getline(tSequenceStream, strLine);
		}
		
		// get the sequence name if one
		if ((false == tSequenceStream.eof()) && ('>' == strLine[0]))
		{
			// extract structure name
			tNewSecondaryStructure.SetName(strLine.substr(1));
			// get next line
			std::getline(tSequenceStream, strLine);
		}

		// init primary structure
		bool bInNucleicAcidChain = true;
		std::string strNAChain;
		int iSplitPosition = 0;
		while ((false == tSequenceStream.eof()) && (true == bInNucleicAcidChain))
		{
			int iInLinePosition = 0;
			if (0 < strLine.size())
			{
				while (((signed)strLine.size() > iInLinePosition)
				       && ((val::biocpp::nucleic_acids::Unknown != val::biocpp::CharToSNucleotideCodeSafe(strLine[iInLinePosition])
				       	   || (val::biocpp::na_secondary_structure::C_STEM_SEPARATOR == strLine[iInLinePosition]))))
				{
					// check for stem split position
					if (val::biocpp::na_secondary_structure::C_STEM_SEPARATOR == strLine[iInLinePosition])
					{
						// check if we already got a split position
						if (0 != iSplitPosition)
						{
							std::ostringstream tOutputStream;
							tOutputStream << "SNASSDotBracketLoader: The given sequence contains more than one split position (second split position detected at position " << strNAChain.size()+1 << ")!";
							throw std::invalid_argument(tOutputStream.str());
						}
						else
						{
							iSplitPosition = strNAChain.size() + 1;
						}
					}
					else
					{
						// we got a valid nucleotide
						strNAChain += strLine[iInLinePosition];
					}
					++iInLinePosition;
				}
				// check if we got invisible characters to skip
				while (((signed)strLine.size() > iInLinePosition)
				       && ((' ' == strLine[iInLinePosition])
				           || ('\t' == strLine[iInLinePosition])
				           || ('\n' == strLine[iInLinePosition])
				           || ('\r' == strLine[iInLinePosition])))
				{
					++iInLinePosition;
				}
				// check if we reached the end of line
				if ((signed)strLine.size() > iInLinePosition)
				{
					// end of line not reached, check for a comment
					if ('#' != strLine[iInLinePosition])
					{
						// encountered an unwanted character!
						std::ostringstream tOutputStream;
						tOutputStream << "SNASSDotBracketLoader: The given sequence contains invalid character at position " << strNAChain.size()+1 << " (character: \"" << strLine[iInLinePosition] << "\")!";
						throw std::invalid_argument(tOutputStream.str());
					}
				}
			}
			// get next line
			std::getline(tSequenceStream, strLine);
			// check for dot-bracket annotations
			if ((0 < strLine.size())
				&& (('(' == strLine[0])
				    || ('.' == strLine[0])))
			{
				bInNucleicAcidChain = false;
			}
		}

		// set nucleotides chain
		tNewSecondaryStructure.AppendSequence(strNAChain);

		// init secondary structure
		// get sequence length
		int iSequenceLength = (signed)strNAChain.length();
		int iStructureLength = 0;
		while ((false == tSequenceStream.eof()) && (iStructureLength < iSequenceLength))
		{
			int iInLinePosition = 0;
			if (0 < strLine.size())
			{
				while (((signed)strLine.size() > iInLinePosition)
				       && (('(' == strLine[iInLinePosition])
				           || (')' == strLine[iInLinePosition])
				           || ('.' == strLine[iInLinePosition])
				       	   || (val::biocpp::na_secondary_structure::C_STEM_SEPARATOR == strLine[iInLinePosition])))
				{
					// check structure content
					if ('(' == strLine[iInLinePosition])
					{
						tIndexStack.push(++iStructureLength);
					}
					else if (')' == strLine[iInLinePosition])
					{
						if (true == tIndexStack.empty())
						{
							throw std::invalid_argument("SNASSDotBracketLoader: The given sequence does not contain matching brackets!");
						}
						tNewSecondaryStructure.Pair(tIndexStack.top(), ++iStructureLength);
						tIndexStack.pop();
					}
					else if ('.' == strLine[iInLinePosition])
					{
						++iStructureLength;
					}
					else if (val::biocpp::na_secondary_structure::C_STEM_SEPARATOR == strLine[iInLinePosition])
					{
						// check if split position is correct
						if (iStructureLength+1 != iSplitPosition)
						{
							std::ostringstream tOutputStream;
							tOutputStream << "SNASSDotBracketLoader: The given split position in dot-bracket annotation does not correspond to the position given in the sequence!";
							throw std::invalid_argument(tOutputStream.str());
						}
					}
					else
					{
						// encountered an unwanted character!
						std::ostringstream tOutputStream;
						tOutputStream << "SNASSDotBracketLoader: The given dot-bracket annotated sequence contains invalid characters (others than '(', ')', '.' or '_' (at a good position)) at position " << iInLinePosition << " (\"" << strLine[iInLinePosition] << "\")!";
						throw std::invalid_argument(tOutputStream.str());
					}
					++iInLinePosition;
				}
				// check if we got invisible characters to skip
				while (((signed)strLine.size() > iInLinePosition)
				       && ((' ' == strLine[iInLinePosition])
				           || ('\t' == strLine[iInLinePosition])
				           || ('\n' == strLine[iInLinePosition])
				           || ('\r' == strLine[iInLinePosition])))
				{
					++iInLinePosition;
				}
				// check if we reached the end of line
				if ((signed)strLine.size() > iInLinePosition)
				{
					// end of line not reached, check for a comment
					if ('#' != strLine[iInLinePosition])
					{
						// encountered an unwanted character!
						std::ostringstream tOutputStream;
						tOutputStream << "SNASSDotBracketLoader: The given dot-bracket annotation contains invalid character at position " << iStructureLength << " (character: \"" << strLine[iInLinePosition] << "\")!";
						throw std::invalid_argument(tOutputStream.str());
					}
				}
			}
			// get next line
			std::getline(tSequenceStream, strLine);
		}
		// check for a simple stem
		if (0 != iSplitPosition)
		{
			try
			{
				tNewSecondaryStructure.SetStemSplitPosition(iSplitPosition);
			}
			catch (std::invalid_argument &ex)
			{
				// encountered an unwanted character!
				std::ostringstream tOutputStream;
				tOutputStream << "SNASSDotBracketLoader: Invalid split position! Make sure the structure contains only one stem and the split position (character \"" << val::biocpp::na_secondary_structure::C_STEM_SEPARATOR << "\") is not before an opening pair.";
				throw std::invalid_argument(tOutputStream.str());
			}
		}
		return tNewSecondaryStructure;
	}
};
}; // namespace loaders
}; // namespace biocpp
}; // namespace val
#endif //ifndef __NASSDOTBRACKETLOADER_H
