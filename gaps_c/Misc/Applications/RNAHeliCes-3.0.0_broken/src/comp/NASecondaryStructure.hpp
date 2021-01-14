//*****************************************************************************
#ifndef __NASECONDARYSTRUCUTRE_H
#define __NASECONDARYSTRUCUTRE_H 1

/**
@file NASecondaryStructure.h
****************************
@brief NASecondaryStructure.h is the header file for
@c val::biocpp::CNASecondaryStructure class.

SecondaryStructure.h contains @c val::biocpp::CNASecondaryStructure class
declaration and methods prototypes.

@see val::biocpp::CNASecondaryStructure

@author Valentin GUIGNON
@version 1.0
@date 18/07/2004
*/

#include <map>
#include <list>
#include "BioCPP.hpp"
#include "NAPrimaryStructure.hpp"

#include "TreeNode.hpp"
#include "HelixIndex.hpp"

namespace val
{
namespace biocpp
{

//! Error message when the user specified an invalid stem or stem-loop number.
const std::string STR_ERROR_INVALID_STEM_OR_STEMLOOP = std::string("Error: The stem or stem-loop number specified is not valid!\n");

//! Error message when the user specified an invalid stem or stem-loop number.
const std::string STR_ERROR_CORRUPTED_STRUCTURE_DATA = std::string("Error: Stem or stem-loop data is not valid! Maybe the structure was not parsed correctly.\n");

//! Error message when a function requiers the structure to be parsed while the structure was not parsed.
const std::string STR_ERROR_NOT_PARSED_STRUCTURE = std::string("Error: Secondary structre has not been parsed before this call!\n");

namespace na_secondary_structure
{
	/**
	@name String format constants
	*****************************
	The following constants are used to select a string format of the secondary
	structure when calling @c ToString.

	@see ToString
	*/
	//@{
	//! No secondary structure information, only a simple sequence.
	const int I_ONE_LINE_NO_BRACKETS	        = 0;
	//! Brackets are integrated in the sequence line (non-standard).
	const int I_ONE_LINE_WITH_BRACKETS		    = 1;
	//! Brackets are on a separate line (standard notation).
	const int I_TWO_LINES						= 2;
	//! Display stems and stem-loops in a semi-graphical way.
	const int I_TEXT_STEMLOOPS					= 3;
	//! Display all stems and stem-loops in dot-bracket format.
	const int I_TWO_LINES_STEMLOOPS				= 4;
	//@}

	/**
	@name Nucleotides configuration constants
	*****************************************
	The following constants are used to describe a nucleotide configuration.

	@see GetNucleotideConfig
	*/
	//@{
	//! The nucleotide is in an unknown configuration.
	const int I_UNKNOWN_CONFIG           = 0;
	//! The nucleotide is an external base.
	const int I_EXTERNAL_BASE            = 1;
	//! The nucleotide is in a stacked pair.
	const int I_STACKED_PAIR             = 2;
	//! The nucleotide is in a bulge.
	const int I_BULGE                    = 3;
	//! The nucleotide is in an internal loop.
	const int I_INTERNAL_LOOP            = 4;
	//! The nucleotide is in a hairpin loop.
	const int I_HAIRPIN_LOOP             = 5;
	//! The nucleotide is in a multibranched loop.
	const int I_MULTIBRANCHED_LOOP       = 6;
	//! The nucleotide is in a pseudoknot.
	const int I_PSEUDOKNOT               = 7;
	//! The nucleotide is in a kissing hairpin.
	const int I_KISSING_HAIRPIN	        = 8;
	//! The nucleotide is in a kissing loop.
	const int I_KISSING_LOOP             = 9;
	//! The nucleotide is in a triple helix.
	const int I_TRIPLE_HELIX             = 10;
	//@}


	/**
	@name Structural component constants
	************************************
	The following constants are used to describe structural components nature.

	*/
	//@{
	//! Unknown component nature.
	const int I_UNKNOWN_NATURE = 0;
	//! Unknown component nature.
	const int I_STEM_TYPE      = 1;
	//! Unknown component nature.
	const int I_STEMLOOP_TYPE  = 2;
	//@}


	/**
	@name StemLoop extraction mode constants
	****************************************
	The following constants are used to select a way to extract a stem-loop
	from the structure.

	@see GetStemLoop
	*/
	//@{
	//! Extract stem or stem-loop.
	const int I_STEM_OR_STEM_LOOP                 = 0;
	//! Extract stems only.
	const int I_STEM_ONLY                         = 1;
	//! Extract stem-loop only.
	const int I_STEM_LOOP_ONLY                    = 2;
	//! The stem or stem-loop will not include any free bases at its ends.
	const int I_CUT_AT_PAIR                       = 4;
	//! The stem or stem-loop will include free bases at 5' end.
	const int I_CUT_INCLUDING_5PRIME_FREE_BASES   = 5;
	//! The stem or stem-loop will include free bases at 3' end.
	const int I_CUT_INCLUDING_3PRIME_FREE_BASES   = 6;
	//! The stem or stem-loop will include free bases at both ends.
	const int I_CUT_INCLUDING_NEIGHBOR_FREE_BASES = 7;
	//! The stem or stem-loop will include half of the free bases at both ends.
	const int I_CUT_INCLUDING_HALF_FREE_BASES     = 8;
	//! The stem or stem-loop will include free bases at 5' end.
	const int I_CUT_INCLUDING_5PRIME_FREE_BASES_WITH_LAST_3PRIME_FREE_BASES    = 9;
	//! The stem or stem-loop will include free bases at 3' end.
	const int I_CUT_INCLUDING_3PRIME_FREE_BASES_WITH_FIRST_5PRIME_FREE_BASES   = 10;
	//@}


	//! Char used to separate the two segments of sequences in a stem in dot-bracket representation.
	const char C_STEM_SEPARATOR	= '_';

}; // end namespace CNASecondaryStructure

/**
@class CNASecondaryStructure   NASecondaryStructure.h "NASecondaryStructure.hpp"
****************************
@brief CNASecondaryStructure implements a nucleic acids chain secondary
structure (DNA, RNA) and its basic operations.

This class works with @c C3LinksNucleotide nucleotides.

@author Valentin GUIGNON
@version 1.0
@date 13/03/2007
*/
class CNASecondaryStructure
: public CNAPrimaryStructure
{
public:
	/**
	@struct SSubstructureData
	*****************
	@brief SSubstructureData contains stem or stem-loop data.
	
	@c SSubstructureData contains stem or stem-loop data.

	*/
	struct SSubstructureData
	{
		//! Contains I_STEM_TYPE or I_STEMLOOP_TYPE (see CNASecondaryStructure constants).
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

	//! Default constructor.
	CNASecondaryStructure(void);
	//! Copy constructor.
	CNASecondaryStructure(const CNASecondaryStructure& tSecondaryStructure);
	//! Destructor.
	virtual ~CNASecondaryStructure(void);
	//! Assignation operator.
	virtual CNASecondaryStructure& operator=(const CNASecondaryStructure &tSecondaryStructure);
	//! Inherited assignation operator.
	virtual CNAPrimaryStructure& operator=(const val::biocpp::CNAPrimaryStructure &tSequence);
	//! Concatenates given sequence to the structure
	virtual CNAPrimaryStructure& operator+=(const val::biocpp::CNAPrimaryStructure &tRightSequence);
//	//! Concatenates 2 sequences
//	virtual CNASecondaryStructure operator+(const val::biocpp::CNASecondaryStructure &tRightStructure) const;

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
	//! Inserts iCopyCount copies of tBaseCode at the given position.
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

	//! Pairs two nucleotides together.
	virtual void Pair(int iFirstPosition, int iSecondPosition);
	//! Pairs two intervals of nucleotides together.
	virtual void PairInterval(int iFirstStartPosition, int iSecondStartPosition, int iPairsCount, bool bIncreaseSecondPosition=false);
	//! Unpairs the specified nucleotide (and his pair).
	virtual void BreakPair(int iAtPosition);
	//! Unpairs the specified interval of nucleotides.
	virtual void BreakPairInterval(int iStartPosition, int iPairsCount, bool bIncreasePosition=true);
	//! Returns the pair index if a nucleotide is paired.
	virtual int GetPairedBaseIndex(int iAtPosition) const;
	//! Initializes nodes information (call it before using bases configuration related methods)
	virtual void ParseStructure();
	//! Returns the configuration of a nucleotide.
	virtual int GetNucleotideConfig(int iAtPosition) const;
	//! Returns the tertiary structure configuration of a nucleotide.
	virtual int GetNucleotideTertiaryConfig(int iAtPosition) const;
	//! Tells if a nucleotide is an extern.
	virtual bool IsExternalBase(int iAtPosition) const;
	//! Tells if a nucleotide is paired.
	virtual bool IsPaired(int iAtPosition) const;
	//! Tells if a nucleotide is in a bulge.
	virtual bool IsInBulge(int iAtPosition) const;
	//! Tells if a nucleotide is in an internal loop.
	virtual bool IsInInternalLoop(int iAtPosition) const;
	//! Tells if a nucleotide is in a hairpin loop.
	virtual bool IsInHairpinLoop(int iAtPosition) const;
	//! Tells if a nucleotide is in a multibranched loop.
	virtual bool IsInMultibranchedLoop(int iAtPosition) const;
	//! Tells if a nucleotide is in a pseudoknot loop.
	virtual bool IsInPseudoknot(int iAtPosition) const;
	//! Tells if a nucleotide is in a kissing hairpin loop.
	virtual bool IsInKissingHairpin(int iAtPosition) const;
	//! Tells if a nucleotide is in a kissing loop.
	virtual bool IsInKissingLoop(int iAtPosition) const;
	//! Tells if a nucleotide is in a triple helix.
	virtual bool IsInTripleHelix(int iAtPosition) const;
	//! Returns the number of external bases.
	virtual int GetExternalBasesCount() const;
	//! Returns the number of bases in stacked pairs (a pair is 2 nucleotides).
	virtual int GetStackedPairBasesCount() const;
	//! Returns the number of bulges.
	virtual int GetBulgeBasesCount() const;
	//! Returns the number of bases in internal loops.
	virtual int GetInternalLoopBasesCount() const;
	//! Returns the number of bases in hairpin loops.
	virtual int GetHairpinLoopBasesCount() const;
	//! Returns the number of bases in multibranched loops.
	virtual int GetMultibranchedLoopBasesCount() const;
	//! Returns the number of base in pseudoknots.
	virtual int GetPseudoknotBasesCount() const;
	//! Returns the number of base in kissing hairpin loops.
	virtual int GetKissingHairpinBasesCount() const;
	//! Returns the number of base in kissing loops.
	virtual int GetKissingLoopBasesCount() const;
	//! Returns the number of bases in triple helix.
	virtual int GetTripleHelixBasesCount() const;
	//! Returns the number of stem-loops.
	virtual int GetStemLoopsCount() const;
	//! Returns the number of stem (without hairpin loops).
	virtual int GetStemsCount() const;
	//! Returns the number of stem and stem-loops.
	virtual int GetStemsAndStemLoopsCount() const;
	//! Returns the secondary structure as an arc-annotated sequence.
	virtual std::string ToString(int iDisplayMode=na_secondary_structure::I_ONE_LINE_NO_BRACKETS) const;
	//! Tells if given stem/stem-loop index correspond to a stem.
	virtual bool IsStem(int iStemLoopNumber) const;
	//! Tells if given stem/stem-loop index correspond to a stem-loop.
	virtual bool IsStemLoop(int iStemLoopNumber) const;
	//! Returns the secondary structure of a stem or a stem-loop.
	virtual CNASecondaryStructure GetStemOrStemLoop(int iStemLoopNumber, int iCutMethod=val::biocpp::na_secondary_structure::I_CUT_AT_PAIR) const;
	//! Returns the position of the end of the first segment of a stem or 0 for non-stem.
	virtual int GetStemSplitPosition() const;
	//! Sets the position of the end of the first segment of a stem.
	virtual void SetStemSplitPosition(int iPosition);
	//! Returns data about a stem or a stem-loop.
	virtual SSubstructureData GetStemOrStemLoopData(int iStemLoopNumber, int iCutMethod=val::biocpp::na_secondary_structure::I_CUT_AT_PAIR) const;
	
	typedef val::biocpp::tools::CTreeNode<hishape::comp::CHelixIndex> CHelixIndexTreeNode;
	virtual CHelixIndexTreeNode::CTreeNodePointer ParseStructureWithTree();

protected:
	//! Tells if the structure has been parsed or not.
	bool m_bIsStructureParsed;
	int m_iExternalBasesCount;
	int m_iStackedPairBasesCount;
	int m_iBulgeBasesCount;
	int m_iInternalLoopBasesCount;
	int m_iHairpinLoopBasesCount;
	int m_iMultibrachedLoopBasesCount;
	int m_iPseudoknotBasesCount;
	int m_iKissingHairpinBasesCount;
	int m_iKissingLoopBasesCount;
	int m_iTripleHelixBasesCount;
	int m_iStemSplitPosition;
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
	/**
	m_tNucleotidesTertiaryInteractions
	**********************************
	m_tNucleotidesTertiaryInteractions follows the same scheme as
	m_tNucleotidesPairs. m_tNucleotidesTertiaryInteractions does not store
	pair index but tertiary interaction. For example, if base 4, 17 and 38
	are linked together, the strongest interaction, let's say 4-17 will be
	stored in m_tNucleotidesPairs and the interaction with 38 will be stored
	in m_tNucleotidesTertiaryInteractions as following:
	m_tNucleotidesPairs[4] == 17
	m_tNucleotidesPairs[17] == 4
	m_tNucleotidesPairs[38] == 0
	m_tNucleotidesTertiaryInteractions[4] == 38
	m_tNucleotidesTertiaryInteractions[17] == 38
	m_tNucleotidesTertiaryInteractions[38] == 4 (or 17, one of the 2 other bases, doesn't matter wich one)

	Note: m_tNucleotidesTertiaryInteractions must always be GetLength() + 1
	      even if the structure has not been parsed yet.
	*/
	std::deque<int> m_tNucleotidesTertiaryInteractions;

	/**
	m_tNucleotidesConfig
	********************
	1-based array. First element is not used (wasted). The length of this
	array is set in ParseStructure() so it can only be used when
	m_bIsStructureParsed is true.
	*/
	std::deque<int> m_tNucleotidesConfig;

	/**
	m_tNucleotidesConfigTertiaryStructure
	*************************************
	1-based array. First element is not used (wasted). The length of this
	array is set in ParseStructure() so it can only be used when
	m_bIsStructureParsed is true.
	*/
	std::deque<int> m_tNucleotidesConfigTertiaryStructure;
	int m_iStemsCount;
	int m_iStemLoopsCount;

	/**
	m_tStemsAndStemLoopsNature
	**************************
	Associate to a stem/stem-loop number a nature (ie. 'stem' or 'stem-loop' or 'unknown').
	Can only be used when m_bIsStructureParsed is true.
	*/
	std::deque<int> m_tStemsAndStemLoopsNature;

	/**
	m_tStemsAndStemLoopsStartingBases
	*********************************
	m_tStemsAndStemLoopsStartingBases has 2 way of use:
	-for stem, 'pair.first' is the first base of the first base-pair and
	 'pair.second' is the last base of the last base-pair;
	-for stem-loop, 'pair.first' is the first base of the first base-pair and
	 'pair.second' is the second base of the first base-pair.
	m_tStemsAndStemLoopsStartingBases.size() is the number of stem and
	stem-loop found in the parsed structure.
	Can only be used when m_bIsStructureParsed is true.
	Important note: stored base index always belong to paired bases.
	*/
	std::deque<std::pair<int, int> > m_tStemsAndStemLoopsStartingBases;

	//! type of the iterator on the array of nucleotides.
	typedef std::deque<int>::iterator CDequeIntIterator;
	//! returns an iterator at the given position.
	CDequeIntIterator GetPairsDataIteratorAt(int iAtPosition);
	//! returns an iterator at the given position.
	CDequeIntIterator GetTertiaryInteractionsDataIteratorAt(int iAtPosition);
	
	


    
}; // class CNASecondaryStructure
}; // namespace biocpp
}; // namespace val
#endif //ifndef __NASECONDARYSTRUCUTRE_H
