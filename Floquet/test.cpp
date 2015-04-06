//
// Created by Liangsheng Zhang on 3/30/15.
//

#include <iostream>
#include <ctime>
#include <string>
#include "constants.h"
#include "parameters.h"
#include "tasks_models.h"
#include "randomc.h"

using namespace std;

CRandomMersenne RanGen_mersenne(time(NULL));
TasksModels tasks_models; // Record all the tasks and methods

int main(){

    AllPara parameters;

    parameters.generic.task = "Flo Eigen";
    parameters.generic.model = "XXZ Random Simple Shift Flo";

    parameters.generic.size = 2; // System size
    parameters.generic.evec = true; // Whether compute eigenvectors, so far only called in
    // level statistics calculation
    parameters.generic.erase = true; // Whether erase matrix after diagonization
    parameters.generic.debug = true; // Whether output debug information
    parameters.generic.iso_keep = true; // Whether isolated part is kept
    parameters.generic.version = 1; // Version of the output

    parameters.output.width = 30; // Width for spacing in output files
    parameters.output.filename_output = true; // Whether print out file names

    parameters.floquet.J_N = 30; // Number of points of coupling strength
    parameters.floquet.J_min = 0.1; // Minimum J
    parameters.floquet.J_max = 0.9; // Maximum J
    parameters.floquet.tau = 0.8; // Time step size
    parameters.floquet.J = 0.7;

    parameters.markov.K = 0.8; // Coupling strength to the bath in markov models

    string model = parameters.generic.model;
    EvolMatrix<ComplexEigenSolver<MatrixXcd> >* floquet;
    tasks_models.Model(model, parameters, floquet);
    floquet -> Evol_Para_Init();
    floquet -> Evol_Construct();
    floquet -> Evol_Diag();

    tasks_models.Print_Model();
    tasks_models.Print_Task();

    return 0;
}

