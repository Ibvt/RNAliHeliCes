//*****************************************************************************
#include "Algorithm.hpp"
#include "DebugTools.hpp"


namespace val
{
namespace biocpp
{
namespace algorithms
{
/**
Default Constructor
*******************
Constructs a new @c CAlgorithm object that does nothing.

@return CAlgorithm:\n
 a new instance of CAlgorithm.
*/
CAlgorithm::CAlgorithm(const bool &bThreaded)
{
	TRACEI("CAlgorithm::CAlgorithm(const bool &bThreaded)");
}

/**
Destructor
**********
A standard virtual destructor.
*/
CAlgorithm::~CAlgorithm(void)
{
	TRACEI("CAlgorithm::~CAlgorithm(void)");
}

}; // namespace algorithms
}; // namespace biocpp
}; // namespace val
