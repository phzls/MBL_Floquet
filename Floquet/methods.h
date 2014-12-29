#ifndef METHODS_H
#define METHODS_H

#include "evol_class.h"
#include "results.h"
#include "parameters.h"

using namespace std;

/**
 ** This file contains methods which are used in evolution/level statistics.
 **/

/*
 * Compute level statistics where the results are redirected to a data of type T2.
 */
template <class T1, class T2>
void level_cal(const AllPara&, vector<EvolMatrix<T1>*>&, ResultsOutput<EvolMatrix<T1>*, T2>*,
	T2&);

/*
 * Compute level statistics and print out the results.
 */
template <class T1, class T2>
void level_cal(const AllPara&, vector<EvolMatrix<T1>*>&, 
	ResultsOutput<EvolMatrix<T1>*, T2>*);

#include "level_cal.tpp"

#endif