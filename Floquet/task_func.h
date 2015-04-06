#ifndef TASK_FUNC_H
#define TASKS_FUNC_H

#include "parameters.h"

/**
 ** This file includes all the tasks that can be called from the main function.
 **/

/*
 * Compute level statistics for floquet system.
 */
void flo_level(const AllPara&);

/*
 * Compute entry and norm representations of the right most sigma matrix at the end of chain
 * for floquet systems.
 */
void flo_rightmost_sigma_z(const AllPara&);

/*
 * Compute the time evolution for a floquet system of pure states
 */
void flo_evolution(const AllPara&);

/*
 * Compute the simple Markov time evolution for a floquet system
 */
void flo_evolution_simple_markov(const AllPara&);

/*
 * Compute entry and norm representations of the leftmost sigma matrix at the end of chain
 * for floquet systems.
 */
void flo_leftmost_sigma_z(const AllPara&);

/*
 * Compute entry and norm representations of the leftmost and rightmostsigma matrix at the 
 * end of chain for floquet systems.
 */
void flo_chain_end_sigma_z(const AllPara&);

/*
 * Compute the simple Markov time evolution for a floquet system with only one model and 
 * possibly multiple sets of initial parameters
 */
void flo_evolution_simple_markov_one_model(const AllPara&);

/*
 * Call other functions under a single model
 */
void single_model(const AllPara&);

/*
 * Compute the time evolution for a floquet system of density matrices
 */
void flo_evolution_density(const AllPara&);

/*
 * Study transition properties of a floquet model
 */
void flo_transition(const AllPara&);

/*
 * Eigenstates related properties of a floquet model
 */
void flo_eigen(const AllPara&);

#endif