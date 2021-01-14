//*****************************************************************************
#ifndef __ALIGNER_H
#define __ALIGNER_H 1


/**
@file Aligner.h
***************
@brief Aligner.h is the header file for @c val::biocpp::CAligner class.

Aligner.h contains @c val::biocpp::algorithms::CAligner class declaration and
methods prototypes.

@see val::biocpp::algorithms::CAligner

@author Valentin GUIGNON
@version 1.0
@date 06/09/2004
*/
#include <string>
#include <ctime>
#include "Algorithm.hpp"
#include "ScoreScheme.hpp"


namespace val
{
namespace biocpp
{
namespace algorithms
{
/**
@class CAligner   Aligner.h "Aligner.hpp"
***************
@brief This class is the base class for alignement algorithm.

@c CAligner implements the base class for alignment algorithm. It contains the
main requiered methods to align sequences.

@author Valentin GUIGNON
@version 1.1
@date 29/08/2007
*/
class CAligner
: public CAlgorithm
{
public:
	//! Default constructor.
	CAligner(void);
	//! Copy constructor.
	CAligner(const CAligner &tAligner);
	// Constructor with a score scheme.
	CAligner(const hishape::comp::CScoreScheme &tScoreScheme);
	//! Destructor
	virtual ~CAligner(void);
	//! Alignment function that returns the alignment score.
	virtual float Align(bool bStoreDistMatrix=true, bool bComputeBacktracking=true) = 0;
	//! Sets the score scheme to use.
	virtual void SetScoreScheme(const hishape::comp::CScoreScheme &tScoreScheme);
	//! Returns current score scheme.
	virtual const hishape::comp::CScoreScheme &GetScoreScheme() const;
	//! Returns a string containing algorithm output.
	virtual const std::string &GetOutput() const;
	//! Returns a string containing the distance matrix.
	virtual const std::string &GetDistanceMatrix() const;
	//! Returns a string containing the backtracking matrix.
	virtual const std::string &GetBacktrackingMatrix() const;
	//! Returns a string containing the alignment.
	virtual const std::string &GetAlignment() const;
	//! Returns the last score computed by @c Align.
	virtual float GetScore() const;
	//! Returns the cycles count of the last call of @c Align.
	virtual int GetCycles() const;
	//! Returns the estimated number of cycles needed.
	virtual int EstimateCycles() const;
	//! Assignation operator.
	virtual CAligner& operator=(const CAligner &tAligner);

protected:
	//! Member score scheme.
	hishape::comp::CScoreScheme* m_ptScoreScheme;
	//! Last score.
	float m_fLastScore;
	//! Algorithm string output.
	std::string m_strOutput;
	//! Distance matrix
	std::string m_strDistanceMatrix;
	//! backtracking matrix
	std::string m_strBacktrackingMatrix;
	//! Alignment.
	std::string m_strAlignment;
	//! Cycles.
	int m_iCyclesCount;
}; // class CAligner
}; // namespace algorithms
}; // namespace biocpp
}; // namespace val
#endif //ifndef __ALIGNER_H
