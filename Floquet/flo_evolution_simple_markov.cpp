#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <omp.h>
#include "evol_class.h"
#include "parameters.h"
#include "initial_obj.h"
#include "evol_data.h"
#include "evol_data_total.h"
#include "tasks_models.h"
#include "eigen_output.h"
#include "flo_evol_model.h"

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

	// The jump time in markov time evolution
	const int markov_time_jump = parameters.evolution.markov_time_jump; 
	const bool markov_jump = parameters.evolution.markov_jump;

	InitObj init_obj;

	EvolData evol_data(parameters);
	EvolDataTotal evol_data_total(parameters);

	// Parallelize model part instead of realization part when model_num > realization
	bool model_parallel = false;
	if (model_num > num_realization) model_parallel = true;

	// A string used to record model for last output
	string model_type;

	// Parallel the models, assuming time evolution is still serial
	#pragma omp parallel for num_threads(threads_N) if (model_parallel)
	for (int i=0; i<model_num;i++){
		cout << i << "th model:" << endl;

		InitInfo init_info; // Information used for initial state

		init_info.size = parameters.generic.size;
		init_info.norm_delta = 1.0e-15;
		init_info.debug = debug;
		init_info.leftmost_spin_z_index = parameters.evolution.leftmost_spin_z_index;

		EvolMatrix<ComplexEigenSolver<MatrixXcd> >* floquet;

		tasks_models.Model(model, parameters, floquet);
		init_info.dim = floquet -> Get_Dim();

		cout << "Diagonalize Evolution Operator." << endl;
		floquet -> Evol_Para_Init();
		floquet -> Evol_Construct();
		floquet -> Evol_Diag();

		// Record eigensystems of the corresponding isolated system
		init_info.complex_eigen.resize(0);
		for (int j=0; j<floquet -> eigen.size(); j++){
			if (debug){
				cout << "Sector " << j << endl;
				cout << "Matrix:" << endl;
				complex_matrix_write(floquet -> Get_U(j));
				cout << "Eigenvalues:"<< endl;
				complex_matrix_write(floquet -> eigen[j].eigenvalues());
				cout << endl;
			}

			if (floquet -> eigen_name[j] == "Isolated"){
				init_info.complex_eigen.push_back(& (floquet -> eigen[j]) );
				if (debug) cout << "Attached eigen: " << j << endl;
			}
		}

		const int sec_num = floquet -> Get_Sector_Dim().size();

		// Parallel the initial states when time evolution is still serial
		#pragma omp parallel for num_threads(threads_N) if (!model_parallel)
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

			// A matrix used to advance density matrix if markov_time_jump > 1
			vector<MatrixXcd> markov_adv(sec_num);

			for (int j=0; j<sec_num; j++){
				int row = floquet -> Get_U(j).rows();
				int col = floquet -> Get_U(j).cols();

				if (row != col){
					cout << "Evolution Matrix in " << j <<"th"
						 << " sector is not square." << endl;
					abort(); 
				}

				if (floquet -> eigen_name[j] != "Isolated"){
					markov_adv[j] = MatrixXcd::Zero(row,col);

					for (int l=0; l<row; l++){
					markov_adv[j](l,l) = pow(floquet -> eigen[j].eigenvalues()(l),
						markov_time_jump);
					}

					markov_adv[j] = floquet -> eigen[j].eigenvectors() * markov_adv[j] *
						floquet -> eigen[j].eigenvectors().adjoint();
				}
				else markov_adv[j] = MatrixXcd::Zero(0,0);
			}

			StepInfo info;
			info.model = i;
			info.realization = n;
			info.debug = debug;
			info.delta = init_info.norm_delta;

			cout << "Time evolution starts." << endl;

			for (int t=0; t < time_step; t++){

				info.time = t;

				// This comes first because we start with t=0	
				evol_data.Data_Compute(state_density, info);
				evol_data_total.Data_Compute(state_density, info);

				// Evol the state to t+1
				temp_density = MatrixXcd::Zero(init_info.dim, init_info.dim);

				int true_sec_num = 0; // Sector number without isolated sectorss
				for (int j=0; j<sec_num; j++){
					if (floquet -> eigen_name[j] != "Isolated"){
						if (!markov_jump || markov_time_jump == 1){
							temp_density += floquet -> Get_U(j) * state_density * 
								floquet -> Get_U(j).adjoint();
						}
						else{
							// Not flip the spin at every step in the bath
							temp_density += markov_adv[j] * state_density * 
								markov_adv[j].adjoint();
						}

						true_sec_num ++;
					}
				}

				if (temp_density.rows() != state_density.rows() || 
					temp_density.cols() != state_density.rows()){
					cout << "Size of current state density matrix is wrong." << endl;
					cout << "Expected size: " << state_density.rows() << endl;
					cout << "Obtained rows: " << temp_density.rows() << endl;
					cout << "Obtained cols: " << temp_density.cols() << endl;
					abort();
				}

				state_density = 1.0/double(true_sec_num) * temp_density;

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

		string task_string = parameters.generic.task;
		replace(task_string.begin(), task_string.end(),' ','_');

		evol_data.Data_Output(parameters, floquet -> Repr() + ",Task_" + task_string + ",Init_"
		 + init_string + init_obj.Init_Para_String(init_func_name, init_info));

		cout << endl;
		cout << endl;

		if (i == 0) model_type = floquet -> Repr();

		delete floquet;
		floquet = NULL;
	}

	string init_string = init_func_name;
	replace(init_string.begin(), init_string.end(),' ','_');

	string task_string = parameters.generic.task;
	replace(task_string.begin(), task_string.end(),' ','_');

	evol_data_total.Data_Total_Output(parameters, model_type + ",Task_" + task_string + ",Init_"
	+ init_string + init_obj.Init_Para_String(init_func_name, parameters));

}