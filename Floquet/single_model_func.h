#ifndef SINGLE_MODEL_FUNC
#define SINGLE_MODEL_FUNC

#include <iostream>
#include "parameters.h"
#include "evol_class.h"

using namespace std;

/**
 ** This file contains functions which can be used under single_model calculation. They accept
 ** a set of parameters and a pointer to a model, and return void.
 **/

/*
 * This function computes entry and norm representations of the left and right sigma matrix at
 * the end of chain for a floquet system.
 */
void flo_chain_end_sigma_z_under_one_model(const AllPara&,
EvolMatrix<ComplexEigenSolver<MatrixXcd> >*);

/*
 * This function computes the time evolution of a floquet system under a simple markov scheme.
 * It can have multiple initial conditions. For each set of parameters, there are possible
 * multiple realizations to average with. It directly deals with density matrix is written in
 * basic binary basis. For now, the different markov states (i.e, different evol_op from the
 * floquet) are assumed of equal probability.
 */

void flo_evolution_simple_markov_under_one_model(const AllPara&, 
EvolMatrix<ComplexEigenSolver<MatrixXcd> >*);



#endif