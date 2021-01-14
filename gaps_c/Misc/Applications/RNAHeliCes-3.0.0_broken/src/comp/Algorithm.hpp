//*****************************************************************************
#ifndef __ALGORITHM_H
#define __ALGORITHM_H 1

/**
@file Algorithm.h
*****************
@brief Algorithm.h is the base class header for computational algorithms.

Algorithm.h is the base class header for computational algorithms. It defines
which functions are needed to handle algorithms.
It defines the namespace val::biocpp::algorithms.

@see val::biocpp::algorithms

@author Valentin GUIGNON
@version 1.0
@date 06/09/2004
*/
#include <time.h>


namespace val
{
namespace biocpp
{
/**
@namespace val::biocpp::algorithms
**********************************
@brief algorithms contains computational algorithms.

@c algorithms is a namespace that contains various bioinformatics algorithms
to align sequences, stuctures and more.
*/
namespace algorithms
{
/**
@class CAlgorithm   Algorithm.h "Algorithm.hpp"
*****************
@brief CAlgorithm implements the basic requiered behaviour of an algorithm.

@c CAlgorithm implements the basic requiered behaviour of an algorithm.
*/
class CAlgorithm
{
protected:
	bool m_bStop;
	//+FIXME: thread member
	//+FIXME: time member
	/*
        clock_t t_start;
        string proc_name = "new module";
        t_start = clock();
        clock_t e_time;
        e_time = clock();
        double duration = (double)(e_time - t_start) /CLOCKS_PER_SEC;
	*/

public:
	CAlgorithm(const bool &bThreaded=false);
	virtual ~CAlgorithm(void);
//	virtual void Run(const int &iCyclesCount=-1);
//	virtual void Run(const clock_t &tMaxTimeRun);
//	virtual void Pause();
//	virtual void Resume();
//	virtual void Stop();
//	virtual void Continue(const int &iCyclesCount=-1);
//	virtual void Continue(const clock_t &tMaxTimeRun);
//	virtual void SetThreaded(const bool &bUseThread);
//	virtual bool IsThreaded() const;
//	virtual bool IsDone() const;
//	virtual bool IsComputing() const;
//	virtual double GetCyclesDone() const;
//	virtual double EstimateTotalCyclesCount() const;
//	virtual double EstimateCyclesRemaining() const;
//	virtual const clock_t &GetComputationTime();
}; // class CAlgorithm
}; // namespace algorithms
}; // namespace biocpp
}; // namespace val
#endif //ifndef __ALGORITHM_H
