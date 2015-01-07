#include "generic_func.h"

using namespace std;

/**
 ** Implementation of functions in generic_func.h without templates to avoid duplicates when
 ** included in different files.
 **/

/* 
 * Generic conjugate function
 */
double generic_conj(double val){return val;}
int generic_conj(int val){return val;}

/*
 * Generic norm function
 */
double generic_norm(double val){return val*val;}
int generic_norm(int val){return val*val;}