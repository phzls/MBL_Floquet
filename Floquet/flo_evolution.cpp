#include <iostream>
#include <omp.h>
#include <cmath>
#include "evol_class.h"
#include "parameters.h"
#include "initial_obj.h"
#include "evol_data.h"
#include "tasks_models.h"

using namespace std;

/**
 ** This file implements the time evolution of a floquet system. There will be only one initial
 ** state in any case. If the operator has multi-sector, then the states of different sectors
 ** stored in the initial state in an order consisten with the order of sectors in eigen vector
 ** of evol_class.
 **/

extern TasksModels tasks_models; // Record all the tasks and methods. Defined in main.

void flo_evolution(const AllPara& parameters){ 
	// Number of realizations.
	const int num_realization = parameters.generic.num_realizations;
	const string model = parameters.generic.model;
	const bool erase = parameters.generic.erase; // Whether erase the evolution matrix
	const bool debug = parameters.generic.debug; // Whether print debug information	
	const int threads_N = parameters.generic.threads_N; // Number of threads for parallelization

	const int width = parameters.output.width; // Width in output file

	const int time_step = parameters.evolution.time_step; // Number of time steps
	const int model_num = parameters.evolution.model_num; // Number of models
	const int jump = parameters.evolution.jump; // jump of time points
	const string init_func_name = parameters.evolution.init_func_name; 

	EvolMatrix<ComplexEigenSolver<MatrixXcd> >* floquet;
	TransitionMatrix transition;
	VectorXcd init_state; // Initial state

	InitObj init_obj;
	InitInfo init_info; // Information used for initial state

	init_info.size = parameters.generic.size;
	init_info.norm_delta = 1.0e-15;

	EvolData evol_data(parameters);

	for (int i=0; i<model_num;i++){
		cout << i << "th model:" << endl;

		cout << "Initialize Model." << endl;
		tasks_models.Model(model, parameters, floquet);
		init_info.dim = floquet -> Get_Dim();

		cout << "Diagonalize Evolution Operator." << endl;
		floquet -> Evol_Para_Init();
		floquet -> Evol_Construct();
		floquet -> Evol_Diag();
		if (erase) floquet -> Evol_Erase();

		cout << "Construct Transition Matrix." << endl;
		floquet -> Transition_Compute(transition, "Basic_Full");

		for (int n=0; n<num_realization; n++){
			cout << endl;
			cout << n << "th realization:" << endl;

			cout << "Construct Initial State." << endl;
			init_state = VectorXcd::Zero(floquet -> Get_Dim());
			init_obj.Init_Func(init_func_name)(init_info, transition, init_state);

			cout << "Time evolution starts." << endl;

			#pragma omp parallel num_threads(threads_N)
			{
				#pragma omp for
				for (int t=0; t < time_step; t++){
					int power = t*jump;

					VectorXcd state_evec(init_state.size()); // Current state in evec basis
					VectorXcd state_basic(init_state.size()); // Current state in binary basis

					int index = 0;
					for (int j=0; j<floquet -> eigen.size(); j++){
						for (int k=0; k<floquet -> eigen[j].eigenvalues().rows(); k++){
							complex<double> eval = floquet -> eigen[j].eigenvalues()(k);
							state_evec(index) = pow(eval, power) * init_state(index);
							index ++;
						}
					}

					state_basic = transition.Matrix("Basic_Full") * state_evec;

					if (state_basic.size() != init_state.size()){
						cout << "Size of current state in binary basis is wrong." << endl;
						cout << "Expected size: " << init_state.size() << endl;
						cout << "Obtained size: " << state_basic.size() << endl;
						abort();
					}

					StepInfo info;
					info.model = i;
					info.realization = n;
					info.time = t;

					evol_data.Data_Compute(state_basic, info);
				}
			}

			cout << "Time evolution ends." << endl;

		}

		cout << "Output data." << endl;
		evol_data.Data_Output(parameters, floquet -> Repr());

		cout << endl;
		cout << endl;

		delete floquet;
		floquet = NULL;
		transition.Erase_All();
	}

}