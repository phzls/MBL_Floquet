#include <cmath>
#include <Eigen/Eigenvalues>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "evol_data.h"
#include "methods.h"
#include "generic_func.h"
#include "eigen_output.h"

using namespace std;
using namespace Eigen;

/**
 ** This file contains functions related to leftmost_spin_z_per_model. It is the average of
 ** leftmost spin z at a spin chain model whose each site has only spin down and spin up. 
 ** Spin up is labeled as 1 and spin down is labeled as -1.
 **/

/*
 * Initialize the leftmost_spin_z. The outer index is time; the inner index is realization
 */
 void EvolData::Leftmost_Spin_Z_Per_Model_Init_(const AllPara& parameters){
	const int num_realizations = parameters.generic.num_realizations;
	const int time_step = parameters.evolution.time_step;

	leftmost_spin_z_per_model_.resize(time_step);

	for (int i=0; i<time_step; i++){
		leftmost_spin_z_per_model_[i].resize(num_realizations);

		for (int j=0; j<num_realizations; j++){
			leftmost_spin_z_per_model_[i][j] = 0;
		}
	}
}

/*
 * Compute the leftmost_spin_z at a given time step and realization, given a complex density
 * matrix in basic binary basis, where 0 is the rightmost position
 */
 void EvolData::Leftmost_Spin_Z_Per_Model_Cal_C_(const MatrixXcd& density_matrix, 
 const StepInfo& info){
	const int realization = info.realization;
	const int time = info.time;

	// Check whether density_matrix is Hermitian
	if (density_matrix.rows() != density_matrix.cols()){
		cout << "Density matrix passed in entropy_per_model_cal_C is not square." << endl;
		cout << "Rows: " << density_matrix.rows() << endl;
		cout << "Cols: " << density_matrix.cols() << endl;
		abort();
	}

	const int row_num = density_matrix.rows();
	const int size = __builtin_popcount(row_num - 1); // Size of the chain
	const int down = 1 << (size-1); // Number below this will have leftmost spin down

	leftmost_spin_z_per_model_[time][realization] = 0;

	for (int i=0; i<row_num; i++){
		if (i<down) leftmost_spin_z_per_model_[time][realization] -= real(density_matrix(i,i));
		else leftmost_spin_z_per_model_[time][realization] += real(density_matrix(i,i)); 
	}

	if (info.debug){
		cout << "Average leftmost spin z per model:" << endl;
		cout << leftmost_spin_z_per_model_[time][realization] << endl;
		cout << endl;
	}
}

void EvolData::Leftmost_Spin_Z_Per_Model_Out_(const AllPara& parameters, const string& name){
	const int time_step = parameters.evolution.time_step;
	const int num_realizations = parameters.generic.num_realizations;
	const int version = parameters.generic.version;
	const int jump = parameters.evolution.jump;
	const bool output = parameters.output.filename_output;
	const int width = parameters.output.width;
	const double step_size = parameters.evolution.step_size;
	const bool log_time = parameters.evolution.log_time;
	const int log_time_jump = parameters.evolution.log_time_jump;
	const bool markov_jump = parameters.evolution.markov_jump;
	const int markov_time_jump = parameters.evolution.markov_time_jump;

	if (leftmost_spin_z_per_model_.size() != time_step){
		cout << "Not enough time steps are computed for leftmost_spin_z." << endl;
		cout << "Expected: " << time_step << endl;
		cout << "Computed: " << leftmost_spin_z_per_model_.size() << endl;
	}

	for(int i=0; i<time_step; i++){
		if (num_realizations != leftmost_spin_z_per_model_[i].size()){
			cout << "leftmost_spin_z is not computed with enough realizations for time " << i 
				 << endl;
			cout << "Expect: " << num_realizations << endl;
			cout << "Computed: " << leftmost_spin_z_per_model_[i].size() << endl;
			abort();
		}
	}

	vector<double> mean(time_step);
	vector<double> sd(time_step);

	stringstream filename;
	filename << name <<",Run=" << num_realizations << ",Total=" << time_step;

	if (log_time) filename <<",log_time";

	if (markov_jump) filename <<",markov_jump=" << markov_time_jump;

	filename << ",jump=" << jump << ",left_spin_z_per_model";
	if (version > 0) filename <<",v" << version;
	filename << ".txt";

	if (output) cout << filename.str() <<endl;

	ofstream fout( filename.str().c_str() );

	for (int t=0; t<time_step; t++){
		double time = t * step_size * jump; 

		if (markov_jump) time *= markov_time_jump;

		if (log_time){
			long long int power = pow(log_time_jump,t);
			time = step_size * power;
		}
		
		generic_mean_sd(leftmost_spin_z_per_model_[t], mean[t], sd[t]);
		fout << setw(10) << time << setw(width) << mean[t] << setw(width) << sd[t] << endl;
	}
}

/*
 * Compute the leftmost_spin_z at a given time step and realization, given a state vector. 
 * Not implemented yet.
 */
 void EvolData::Leftmost_Spin_Z_Per_Model_Cal_(const VectorXcd& init_state, 
 const StepInfo& info){
 	cout << "Leftmost_Spin_Z_Per_Model_Cal_ Not Implemented Yet." << endl;
 	abort();
 }