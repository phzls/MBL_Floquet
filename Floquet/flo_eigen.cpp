//
// Created by Liangsheng Zhang on 4/3/15.
//

#include <ctime>
#include <string>
#include <cstdlib>
#include <iostream>
#include <Eigen/Core>
#include <utility>
#include <sstream>
#include "flo_evol_model.h"
#include "flo_eigen_func.h"
#include "output_func.h"
#include "methods.h"
#include "parameters.h"
#include "tasks_models.h"

using namespace std;

/*
 * Compute quantities related to eigenstates for a Floquet operator. Level statistics haven't been moved into this
 * implementation.
 */

extern TasksModels tasks_models; // Record all the tasks and methods. Defined in main.

void flo_eigen(const AllPara& parameters){

    // System Size
    const int size = parameters.generic.size;
    // Number of realizations.
    const int num_realization = parameters.generic.num_realizations;
    const string model = parameters.generic.model;

    const int width = parameters.output.width; // Output spacing

    AllPara local_parameters(parameters); // Local parameters which can be changed

    stringstream base_filename; // Filename

    ofstream level_out; // Output for level spacings
    ofstream mean_out; // Output for mean of level spacings
    ofstream square_mean_out; // Output for square mean of level spacings




}

