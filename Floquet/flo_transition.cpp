//
// Created by Liangsheng Zhang on 3/24/15.
//

#include <iostream>
#include "parameters.h"
#include "tasks_models.h"
#include "flo_model_transition.h"
#include "eigen_output.h"

using namespace std;

/**
 ** This file implements functions related to studying transition of floquet models from thermal to localization.
 **/

extern TasksModels tasks_models; // Record all the tasks and methods. Defined in main.

void flo_transition(const AllPara& parameters){
    const int J_N = parameters.floquet.J_N; // Number of points for J
    const double J_max = parameters.floquet.J_max; // Maximum J
    const double J_min = parameters.floquet.J_min; // Minimum J
    const int num_realization = parameters.generic.num_realizations;
    const int threads_N = parameters.generic.threads_N;
    const string model = parameters.generic.model;
    const bool debug = parameters.generic.debug;
    string name;

    FloModelTransition flo_model_transition(parameters);

    AllPara local_parameters(parameters); // Local parameters which can be changed


    Eigen::initParallel();

    cout << "Floquet model transition:" << endl;
    for (int i=0; i<J_N; i++){
        double J;
        if (J_N > 1) J = J_min + i * (J_max - J_min)/(J_N-1);
        else J = J_min;

        local_parameters.floquet.J = J;

        cout << "J: " << local_parameters.floquet.J << endl;
        cout << "Start computation." << endl;

        #pragma omp parallel num_threads(threads_N)
        {
            #pragma omp for
            for (int k=0; k < num_realization; k++){
                EvolMatrix<ComplexEigenSolver<MatrixXcd> >* floquet;
                tasks_models.Model(model, local_parameters, floquet);
                floquet -> Evol_Para_Init();
                floquet -> Evol_Construct();
                floquet -> Evol_Diag();
                floquet -> Evol_Erase(); // Time evolution operator is not kept

                if (debug){
                    cout << "Realization " << k << endl;
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

                if (i == 0 && k == 0) name = floquet -> Type();

                LocalInfo local_info;
                local_info.J_index = i;
                local_info.realization_index = k;

                flo_model_transition.Compute(parameters, floquet, local_info);

                delete floquet;
                floquet = NULL;
            }
        }
    }

    cout << "Output data." << endl;
    flo_model_transition.Output(parameters, name);
}

