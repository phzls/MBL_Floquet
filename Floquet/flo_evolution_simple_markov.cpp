#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include "evol_class.h"
#include "parameters.h"
#include "initial_obj.h"
#include "evol_data.h"
#include "tasks_models.h"
#include "eigen_output.h"

using namespace std;

/**
 ** This file implements the time evolution of a floquet system under a simple markov scheme. 
 ** It directly deals with density states at any step. The density matrix is written in basic
 ** binary basis. For now, the different markov states (i.e, different evol_op from the 
 ** floquet) are assumed of equal probability.
 **/

extern TasksModels tasks_models; // Record all the tasks and methods. Defined in main.

void flo_evolution_simple_markov(const AllPara& parameters){ 
	// Number of realizations.
	const int num_realization = parameters.generic.num_realizations;
	const string model = parameters.generic.model;
	const bool debug = parameters.generic.debug; // Whether print debug information	
	const int threads_N = parameters.generic.threads_N; // Number of threads for parallelization

	const int width = parameters.output.width; // Width in output file

	const int time_step = parameters.evolution.time_step; // Number of time steps
	const int model_num = parameters.evolution.model_num; // Number of models
	const string init_func_name = parameters.evolution.init_func_name; 

	EvolMatrix<ComplexEigenSolver<MatrixXcd> >* floquet;

	InitObj init_obj;
	InitInfo init_info; // Information used for initial state

	init_info.size = parameters.generic.size;
	init_info.norm_delta = 1.0e-15;
	init_info.debug = debug;

	EvolData evol_data(parameters);

	for (int i=0; i<model_num;i++){
		cout << i << "th model:" << endl;

		cout << "Initialize Model." << endl;
		tasks_models.Model(model, parameters, floquet);
		init_info.dim = floquet -> Get_Dim();

		cout << "Diagonalize Evolution Operator." << endl;
		floquet -> Evol_Para_Init();
		floquet -> Evol_Construct();

		const int sec_num = floquet -> Get_Sector_Dim().size();

		// Parallel the initial states, assuming time evolution is still serial
		#pragma omp parallel for num_threads(threads_N)
		for (int n=0; n<num_realization; n++){
			cout << endl;
			cout << n << "th realization:" << endl;

			MatrixXcd state_density; // Density matrix of the state

			cout << "Construct Initial State." << endl;
			state_density = MatrixXcd::Zero(floquet -> Get_Dim(), floquet -> Get_Dim());
			init_obj.Init_Func_C(init_func_name)(init_info, state_density);

			if (state_density.rows() != floquet -> Get_Dim() || 
				state_density.cols() != floquet -> Get_Dim()){
				cout << "The size of initial density matrix is wrong." << endl;
				cout << "Expected size: " << floquet -> Get_Dim() << endl;
				cout << "Rows: " << state_density.rows() << endl;
				cout << "Cols: " << state_density.cols() << endl;
				abort();
			}

			// A temporary holder for density matrix
			MatrixXcd temp_density;

			cout << "Time evolution starts." << endl;

			for (int t=0; t < time_step; t++){

				StepInfo info;
				info.model = i;
				info.realization = n;
				info.time = t;
				info.debug = debug;

				// This comes first because we start with t=0	
				evol_data.Data_Compute(state_density, info);

				// Evol the state to t+1
				temp_density = MatrixXcd::Zero(init_info.dim, init_info.dim);

				for (int i=0; i<sec_num; i++){
					temp_density += floquet -> Get_U(i) * state_density * 
						floquet -> Get_U(i).adjoint();
				}

				if (temp_density.rows() != state_density.rows() || 
					temp_density.cols() != state_density.rows()){
					cout << "Size of current state density matrix is wrong." << endl;
					cout << "Expected size: " << state_density.rows() << endl;
					cout << "Obtained rows: " << temp_density.rows() << endl;
					cout << "Obtained cols: " << temp_density.cols() << endl;
					abort();
				}

				state_density = 1.0/double(sec_num) * temp_density;

				if (debug){
					cout << "Time step:" << t+1 << endl;
					cout << "State density matrix in binary basis:"<<endl;
					complex_matrix_write(state_density);
					cout << endl;
				}

			}

			cout << "Time evolution ends." << endl;

		}

		cout << "Output data." << endl;
		string init_string = init_func_name;
		replace(init_string.begin(), init_string.end(),' ','_');
		evol_data.Data_Output(parameters, floquet -> Repr() + ",Init_" + init_string);

		cout << endl;
		cout << endl;

		delete floquet;
		floquet = NULL;
	}

}