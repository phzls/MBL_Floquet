//
// Created by Liangsheng Zhang on 3/18/15.
//

#include <iostream>
#include <omp.h>
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
 ** This file implements the time evolution of a floquet system in density matrix. There will be
 ** only one initial density matrix in any case. If the operator has multi-sector, then the density
 ** matricex of different sectors are stored in the initial system in an order consistent with the
 ** order of sectors in eigen vector of evol_class.
 **/

extern TasksModels tasks_models; // Record all the tasks and methods. Defined in main.

void flo_evolution_density(const AllPara& parameters){
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

	// Whether time changes logarithmically
	const bool log_time = parameters.evolution.log_time;
	// The base under which time changes logarithmically
	const int log_time_jump = parameters.evolution.log_time_jump;

	EvolMatrix<ComplexEigenSolver<MatrixXcd> >* floquet;
	TransitionMatrix transition;
	MatrixXcd init_density; // Initial density matrix of the system

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
		floquet -> Evol_Diag();
		if (erase) floquet -> Evol_Erase();

		if (debug){
			cout << "Eigenvectors and eigenvalues:" << endl;

			for (int i=0; i<floquet -> eigen.size(); i++){
				cout << "Sector " << i <<" :" << endl;
				cout << "Eigenvectors:" << endl;
				complex_matrix_write(floquet -> eigen[i].eigenvectors());
				cout << endl;
				cout << "Eigenvalues:" << endl;
				complex_matrix_write(floquet -> eigen[i].eigenvalues());
				cout << endl;
			}
		}

		cout << "Construct Transition Matrix." << endl;
		floquet -> Transition_Compute(transition, "Basic_Full");

		if (debug){
			cout << "Full to Basic transtion matrix:" << endl;
			complex_matrix_write(transition.Matrix("Basic_Full"));
			cout << endl;
		}

		for (int n=0; n<num_realization; n++){
			cout << endl;
			cout << n << "th realization:" << endl;

			cout << "Construct Initial State." << endl;
			init_density = MatrixXcd::Zero(floquet -> Get_Dim(), floquet -> Get_Dim());
			init_obj.Init_Func_C(init_func_name)(init_info, init_density);

			// Transform initial density matrix into the basis of evec states
			init_density = transition.Matrix("Basic_Full").adjoint() * init_density
					* transition.Matrix("Basic_Full");

			cout << "Time evolution starts." << endl;

			#pragma omp parallel num_threads(threads_N)
			{
				#pragma omp for
				for (int t=0; t < time_step; t++){
					long long int power = t*jump;

					// If time changes logarithmically
					if (log_time){
						power = pow(log_time_jump,t);
					}

					if (power < 0){
						cout << "Overflow happens in time evolution." << endl;
						abort();
					}

					// Current density matrix in evec basis
					MatrixXcd density_evec(init_density.rows(), init_density.cols());

					// Current density matrix in binary basis
					MatrixXcd density_basic(init_density.rows(), init_density.cols());

					int row_index = 0;
					int col_index = 0;
					for (int j=0; j<floquet -> eigen.size(); j++){
						for (int k=0; k<floquet -> eigen[j].eigenvalues().rows(); k++){
							complex<double> row_eval = conj( floquet -> eigen[j].eigenvalues()[k] );

							for (int l=0; l<floquet -> eigen.size(); l++){
								for (int m = 0; m < floquet -> eigen[l].eigenvalues().rows(); l++){
									complex<double> col_eval = floquet -> eigen[l].eigenvalues()[m];

									density_evec(row_index, col_index) = pow(row_eval, power)
											* init_density(row_index, col_index) * pow(col_eval, power);
									col_index ++;
 								}
							}
							row_index ++;
						}
					}

					density_basic = transition.Matrix("Basic_Full") * density_evec
							* transition.Matrix("Basic_Full").adjoint();

					if (density_basic.rows() != init_density.rows()){
						cout << "# rows of current density matrixx in binary basis is wrong." << endl;
						cout << "Expected # of rows: " << init_density.rows() << endl;
						cout << "Obtained # rows: " << density_basic.rows() << endl;
						abort();
					}

					if (density_basic.cols() != init_density.cols()){
						cout << "# cols of current density matrixx in binary basis is wrong." << endl;
						cout << "Expected # of cols: " << init_density.cols() << endl;
						cout << "Obtained # cols: " << density_basic.cols() << endl;
						abort();
					}

					for (int l=0; l<density_basic.rows();l++){
						for (int m=l; m<density_basic.cols();m++){
							if(abs( conj(density_basic(l,m)) - density_basic(m,l) )<init_info.norm_delta){
								cout << "At row " << l <<" and col " << m << " density matrix in basic"
									<<	" basis is not Hermitian" << endl;
								cout << "(l,m): " << density_basic(l,m) << endl;
								cout << "(m,l): " << density_basic(m,l) << endl;
								abort();
							}
						}
					}

					if (debug){
						cout << "Time step:" << t << endl;
						cout << "Density matrix in binary basis:"<<endl;
						complex_matrix_write(density_basic);
						cout << endl;

						cout << "Density matrix in evec basis:" << endl;
						complex_matrix_write(density_evec);
						cout << endl;
					}

					StepInfo info;
					info.model = i;
					info.realization = n;
					info.time = t;
					info.debug = debug;
					info.left_size = parameters.evolution.left_size;

					evol_data.Data_Compute(density_basic, info);
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

		delete floquet;
		floquet = NULL;
		transition.Erase_All();
	}

}
