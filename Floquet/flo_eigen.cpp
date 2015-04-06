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
#include "eigen_output.h"
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

    // Number of realizations.
    const int num_realization = parameters.generic.num_realizations;
    const string model = parameters.generic.model;
    const bool debug = parameters.generic.debug;
    const int threads_N = parameters.generic.threads_N;
    string name;

    FloEigenFunc flo_eigen_func(parameters);

    #pragma omp parallel num_threads(threads_N)
    {
        #pragma omp for
        for (int i=0; i < num_realization; i++){
            EvolMatrix<ComplexEigenSolver<MatrixXcd> >* floquet;
            tasks_models.Model(model, parameters, floquet);
            floquet -> Evol_Para_Init();
            floquet -> Evol_Construct();
            floquet -> Evol_Diag();
            floquet -> Evol_Erase(); // Time evolution operator is not kept

            if (debug){
                cout << "Realization " << i << endl;
                cout << "Eigenvalues:" << endl;
                for (int i=0; i<floquet -> eigen.size(); i++){
                    cout << "Sector " << i << endl;
                    complex_matrix_write(floquet -> eigen[i].eigenvalues());
                    cout << endl;
                }

                cout << "Eigenvectors:" << endl;
                for (int i=0; i<floquet -> eigen.size(); i++){
                    cout << "Sector " << i << endl;
                    complex_matrix_write(floquet -> eigen[i].eigenvectors());
                    cout << endl;
                }
            }

            if (i == 0) name = floquet -> Type();

            LocalInfo local_info;
            local_info.dim = floquet -> Get_Dim();
            local_info.realization_index = i;

            flo_eigen_func.Compute(parameters, floquet, local_info);

            delete floquet;
            floquet = NULL;
        }
    }

    cout << "Output data." << endl;
    flo_eigen_func.Output(parameters, name);
}

