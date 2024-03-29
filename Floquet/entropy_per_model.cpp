#include <cmath>
#include <Eigen/Eigenvalues>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "evol_data.h"
#include "methods.h"
#include "generic_func.h"
#include "eigen_output.h"

using namespace std;
using namespace Eigen;

/**
 ** This file contains functions related to entropy_per_model. It is the entropy of the 
 ** left part of a spin chain model whose local dimension is 2. The entropy is computed in bits.
 **/

/*
 * Initialize the entropy. The outer index is time; the inner index is realization
 */
void EvolData::Entropy_Per_Model_Init_(const AllPara& parameters){
	const int num_realizations = parameters.generic.num_realizations;
	const int time_step = parameters.evolution.time_step;

	entropy_per_model_.resize(time_step);

	for (int i=0; i<time_step; i++){
		entropy_per_model_[i].resize(num_realizations);

		for (int j=0; j<num_realizations; j++){
			entropy_per_model_[i][j] = 0;
		}
	}
}

/*
 * Compute the entropy at a given time step and realization, given a state vector
 */
void EvolData::Entropy_Per_Model_Cal_(const VectorXcd& state_basic, const StepInfo& info){
	const int realization = info.realization;
	const int time = info.time;
	const int left_size = info.left_size;

	MatrixXcd reduced_density; // Reduced density matrix for the left part
	reduced_density_left_2(state_basic, size_, left_size, reduced_density);

	if (info.debug){
		cout << "Reduced density matrix:" << endl;
		complex_matrix_write(reduced_density);
		cout << endl;
	}

	SelfAdjointEigenSolver<MatrixXcd> density_eigen; // Eigen for reduced density matrix

	density_eigen.compute(reduced_density, false); // Eigenvectors not computed

	entropy_per_model_[time][realization] = 0;

	for (int i=0; i<density_eigen.eigenvalues().rows();i++){
		double eval = density_eigen.eigenvalues()(i);

		if (abs(eval)>1.0e-15)
		{
			if (eval<0){
				cout << "Density matrix has significant negative eigenvalues." << endl;
				cout << eval << endl;
				abort();
			}
			entropy_per_model_[time][realization] += -eval*log2(eval);
		}
	}

	if (info.debug){
		cout << "Entropy per model:" << endl;
		cout << entropy_per_model_[time][realization] << endl;
		cout << endl;
	}
}

/*
 * Compute the entropy at a given time step and realization, given a complex density matrix
 */
void EvolData::Entropy_Per_Model_Cal_C_(const MatrixXcd& density_matrix, const StepInfo& info){
	const int realization = info.realization;
	const int time = info.time;

	// Check whether density_matrix is Hermitian
	if (density_matrix.rows() != density_matrix.cols()){
		cout << "Density matrix passed in entropy_per_model_cal_C is not square." << endl;
		cout << "Rows: " << density_matrix.rows() << endl;
		cout << "Cols: " << density_matrix.cols() << endl;
		abort();
	}

	for (int i=0; i< density_matrix.rows(); i++){
		for (int j=i; j<density_matrix.rows();j++){
			if ( norm( density_matrix(i,j) - conj(density_matrix(j,i)) ) > info.delta ){
				cout << "Density matrix is not Hermitian at (" << i << "," << j << ")." << endl;
				cout << "At (" << i << "," << j << "): " << density_matrix(i,j) << endl;
				cout << "At (" << j << "," << i << "): " <<  density_matrix(j,i) << endl;
				abort();
			}
		}
	}

	SelfAdjointEigenSolver<MatrixXcd> density_eigen; // Eigen for reduced density matrix

	density_eigen.compute(density_matrix, false); // Eigenvectors not computed

	entropy_per_model_[time][realization] = 0;

	for (int i=0; i<density_eigen.eigenvalues().rows();i++){
		double eval = density_eigen.eigenvalues()(i);

		if (abs(eval)>1.0e-12)
		{
			if ( eval*log2(eval) != eval*log2(eval) ){
				cout << "Significant negative eigenvalues of density matrix." << endl;
				cout << "eval: " << eval << endl;
				cout << "Time: " << time << "  Realization: " << realization << endl;
				abort();
			}
			entropy_per_model_[time][realization] += -eval*log2(eval);
		}
	}

	if (info.debug){
		cout << "Entropy per model:" << endl;
		cout << entropy_per_model_[time][realization] << endl;
		cout << endl;
	}
}

/*
 * Output entropy per model. The entroy will be averaged over different realizations, and the
 * sample variance will be computed, where the sample is taken at a fixed time step with different
 * realizations. The time step, average and the standard deviation computed from the sample 
 * variance will be computed. Note this standard deviation is not the standard deviation of 
 * average yet. The name taken in will become part of the output file name.
 */

void EvolData::Entropy_Per_Model_Out_(const AllPara& parameters, const string& name){
	const int time_step = parameters.evolution.time_step;
	const int num_realizations = parameters.generic.num_realizations;
	const int jump = parameters.evolution.jump;
	const bool output = parameters.output.filename_output;
	const int width = parameters.output.width;
	const double step_size = parameters.evolution.step_size;
	const bool log_time = parameters.evolution.log_time;
	const int log_time_jump = parameters.evolution.log_time_jump;
	const bool markov_jump = parameters.evolution.markov_jump;
	const int markov_time_jump = parameters.evolution.markov_time_jump;

	if (entropy_per_model_.size() != time_step){
		cout << "Not enough time steps are computed for entropy." << endl;
		cout << "Expected: " << time_step << endl;
		cout << "Computed: " << entropy_per_model_.size() << endl;
	}

	for(int i=0; i<time_step; i++){
		if (num_realizations != entropy_per_model_[i].size()){
			cout << "Entropy is not computed with enough realizations for time " << i << endl;
			cout << "Expect: " << num_realizations << endl;
			cout << "Computed: " << entropy_per_model_[i].size() << endl;
			abort();
		}
	}

	vector<double> mean(time_step);
	vector<double> sd(time_step);

	stringstream filename;
	filename << name <<",Run=" << num_realizations << ",Total=" << time_step;

	if (log_time) filename <<",log_time";

	if (markov_jump) filename <<",markov_jump=" << markov_time_jump;

	filename << ",jump=" << jump << ",entropy_per_model.txt";

	if (output) cout << filename.str() <<endl;

	ofstream fout( filename.str().c_str() );

	for (int t=0; t<time_step; t++){
		double time = t * step_size * jump; // An overflow still happens for entropy

		if (markov_jump) time *= markov_time_jump;

		if (log_time){
			long long int power = pow(log_time_jump,t);
			time = step_size * power;
		}
		
		generic_mean_sd(entropy_per_model_[t], mean[t], sd[t]);
		fout << setw(10) << time << setw(width) << mean[t] << setw(width) << sd[t] << endl;
	}
}