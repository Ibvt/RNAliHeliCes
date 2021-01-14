//*****************************************************************************
#include "Aligner.hpp"
#include "DebugTools.hpp"


namespace val
{
namespace biocpp
{
namespace algorithms
{
/**
Default constructor
*******************
Returns a new instance of @c CAligner child.

@return CAligner:\n
 a new instance of @c CAligner.
*/
CAligner::CAligner(void)
: m_ptScoreScheme(NULL)
, m_fLastScore(0)
, m_iCyclesCount(0)
{
	TRACEI("CAligner::CAligner(void)");
}


/**
Destructor
**********
Standard virtual destructor.
*/
CAligner::~CAligner(void)
{
	TRACEI("CAligner::~CAligner(void)");
	delete m_ptScoreScheme;
}


/**
Copy constructor
****************
*/
CAligner::CAligner(const CAligner &tAligner)
: m_ptScoreScheme(new hishape::comp::CScoreScheme(*(tAligner.m_ptScoreScheme)))
, m_fLastScore(tAligner.m_fLastScore)
, m_strOutput(tAligner.m_strOutput)
, m_strDistanceMatrix(tAligner.m_strDistanceMatrix)
, m_strBacktrackingMatrix(tAligner.m_strBacktrackingMatrix)
, m_strAlignment(tAligner.m_strAlignment)
{
	TRACEI("CAligner::CAligner(const CAligner &tAligner)");
}


/**
Constructor with score scheme
*****************************
Constructor with a score scheme provided.
*/
CAligner::CAligner(const hishape::comp::CScoreScheme &tScoreScheme)
: m_ptScoreScheme(new hishape::comp::CScoreScheme(tScoreScheme))
, m_fLastScore(0)
{
	TRACEI("CAligner::CAligner(const CScoreScheme &tScoreScheme)");
}


/**
@fn void CAligner::SetScoreScheme(const CScoreScheme &tScoreScheme)
*******************************************************************
Sets the score scheme to use for alignements. This methode is the only one
always called each time the score scheme needs to be changed.

@param const CScoreScheme &tScoreScheme:
 a new score scheme to use.
*/
void CAligner::SetScoreScheme(const hishape::comp::CScoreScheme &tScoreScheme)
{
	TRACEI("void CAligner::SetScoreScheme(const CScoreScheme &tScoreScheme)");
	delete m_ptScoreScheme;
	m_ptScoreScheme = new hishape::comp::CScoreScheme(tScoreScheme);
}


/**
@fn const CScoreScheme &tScoreScheme GetScoreScheme() const
***********************************************************
Returns the score scheme used by current alignment algorithm.

@return const CScoreScheme &:
 a reference to current score scheme.
*/
const hishape::comp::CScoreScheme &CAligner::GetScoreScheme() const
{
	TRACEI("const CScoreScheme &CAligner::GetScoreScheme() const");
	return *m_ptScoreScheme;
}


/**
@fn const std::string &CAligner::GetOutput()
********************************************
Returns a string containing algorithm output.

@return std::string:
 a string containing stuff the algorithm did output.
*/
const std::string &CAligner::GetOutput() const
{
	TRACEI("const std::string &CAligner::GetOutput() const");
	return m_strOutput;
}


/**
@fn const std::string &CAligner::GetDistanceMatrix()
****************************************************
Returns a string containing the distance matrix.

@return std::string:
 a string containing the distance matrix.
*/
const std::string &CAligner::GetDistanceMatrix() const
{
	TRACEI("const std::string &CAligner::GetDistanceMatrix() const");
	return m_strDistanceMatrix;
}


/**
@fn const std::string &CAligner::GetBacktrackingMatrix()
********************************************************
Returns a string containing the backtracking matrix.

@return std::string:
 a string containing the backtracking matrix.
*/
const std::string &CAligner::GetBacktrackingMatrix() const
{
	TRACEI("const std::string &CAligner::GetBacktrackingMatrix() const");
	return m_strBacktrackingMatrix;
}


/**
@fn const std::string &CAligner::GetAlignment()
***********************************************
Returns a string containing the alignment.

@return std::string:
 a string containing the alignment.
*/
const std::string &CAligner::GetAlignment() const
{
	TRACEI("const std::string &CAligner::GetAlignment() const");
	return m_strAlignment;
}


/**
@fn float CAligner::GetScore()
******************************
Returns the last score computed by @c Align.

@return float:
 last alignment score computed by the last @c Align call.
*/
float CAligner::GetScore() const
{
	TRACEI("float CAligner::GetScore() const");
	return m_fLastScore;
}


/**
@fn int CAligner::GetCycles()
*********************************
Returns the cycles of the last call of @c Align.

@return int:
 the cycles count (algorithm loops) of the last call of @c Align.
*/
int CAligner::GetCycles() const
{
	TRACEI("int CAligner::GetCycles() const");
	return m_iCyclesCount;
}


/**
@fn int CAligner::EstimateCycles() const
********************************************
Returns the estimated number of cycles needed if align is run with current
parameters.

@return int:
 Estimated cycles count or 0 if no estimation is available.
*/
int CAligner::EstimateCycles() const
{
	TRACEI("int CAligner::EstimateCycles() const");
	return 0;
}


/**
@fn CAligner& CAligner::operator=(const CAligner &tAligner)
***********************************************************
*/
CAligner& CAligner::operator=(const CAligner &tAligner)
{
	TRACEI("CAligner& CAligner::operator=(const CAligner &tAligner)");
	if (this != &tAligner)
	{
		m_fLastScore = tAligner.m_fLastScore;
		m_strOutput = tAligner.m_strOutput;
		m_strDistanceMatrix = tAligner.m_strDistanceMatrix;
		m_strBacktrackingMatrix = tAligner.m_strBacktrackingMatrix;
		m_strAlignment = tAligner.m_strAlignment;
		SetScoreScheme(tAligner.GetScoreScheme());
	}
	return *this;
}


}; // namespace algorithms
}; // namespace val
}; // namespace biocpp
