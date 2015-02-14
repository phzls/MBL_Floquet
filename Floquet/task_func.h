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
 * Compute the time evolution for a floquet system
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

#endif