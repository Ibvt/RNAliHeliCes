//*****************************************************************************
#ifndef __SCORESCHEME_H
#define __SCORESCHEME_H 1

/**
@file ScoreScheme.h
*******************
@brief ScoreScheme.h is the header file for @c hishape::comp::CScoreScheme class.

ScoreScheme.h contains @c hishape::comp::CScoreScheme class declaration and
methods prototypes.

@see hishape::comp::CScoreScheme
*/

#include <vector>
#include <string>
#include <functional>


namespace hishape
{
namespace comp
{
/**
@class CScoreScheme   ScoreScheme.h "ScoreScheme.hpp"
*******************
@brief This class is the base class for alignement score computation.

@c CScoreScheme is the base class for score scheme handling objects. It
contains the main functions needed to work with scoring.
*/
class CScoreScheme
{
protected:
	//! Score matrix.
	std::vector<std::vector<float> > m_tScoreMatrix;
	//! Score matrix width.
	int m_iScoreMatrixWidth;
	//! Score matrix height.
	int m_iScoreMatrixHeight;
	

public:
	//! Default constructor.
	CScoreScheme(void);
	//! Copy constructor.
	CScoreScheme(const CScoreScheme &tScoreScheme);
	//! Constructor with parameters.
	CScoreScheme(const int &iWidth, const int &iHeight);
	//! Destructor.
	virtual ~CScoreScheme(void);

	//! Sets the score for matching a reference element with another.
	virtual void SetScore(const int &iFirstElement, const int &iSecondElement, const float &fScore);
	//! Gets the score of matching a reference element with another.
	virtual const float& GetScore(const int &iFirstElement, const int &iSecondElement) const;
	//! Returns the score matrix as a multi-lines string.
	virtual std::string ToString() const;
	//! Assignation operator.
	virtual CScoreScheme& operator=(const CScoreScheme &tRightScoreScheme);
	

}; // class CScoreScheme
}; // namespace comp
}; // namespace hishape
#endif //ifndef __SCORESCHEME_H
