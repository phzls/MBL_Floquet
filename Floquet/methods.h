#ifndef METHODS_H
#define METHODS_H

#include <utility>
#include <vector>
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

/*
 * Construct the matrix representation of a sigma z operator at rightmost site of the chain.
 * T1 gives the type of the matrix represation. T2 is the type of basis vectors which the
 * representation is about to be written in. The basis vectors must be written in binary basis,
 * i.e., the basis of eigenstates of sigma_z matrix at each site, and with spin down = 0, 
 * spin_up = 1, and the amplitudes are listed in ascending order of the corresponding binary
 * numbers. The string can take values as: "entry" and "norm", which specifies what the matrix
 * store, i.e., the exact entry or the norm. If the T1 object passed in has dimension 0, then 
 * it will be resized to the appropriate size, otherwise a check of size will be performed, and
 * new values will be added to the existing values.
 */

template <class T1, class T2>
void rightmost_sigma_z_sum(T1&, const vector< vector<T2> >&, const string&);

#include "level_cal.tpp"
#include "rightmost_sigma_z_sum.tpp"

#endif