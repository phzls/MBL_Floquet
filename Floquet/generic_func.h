#ifndef GENERIC_H
#define GENERIC_H

#include <complex>

using namespace std;

/**
 ** This file contains some small functions which generalize some functions from c++ so that 
 ** they can take a larger range of arguments.
 **/

/* 
 * Generic conjugate function which also handles double and integer
 */
template<class T>
complex<T> generic_conj(complex<T> val){ return conj(val);}

double generic_conj(double);
int generic_conj(int);

/*
 * Generic norm function which also handles double and integer
 */
template<class T>
T generic_norm(complex<T> val){ return norm(val);}

double generic_norm(double);
int generic_norm(int);
#endif