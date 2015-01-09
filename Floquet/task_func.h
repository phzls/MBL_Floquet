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
void flo_rightmost_sigma_z(const AllPara& parameters);

#endif