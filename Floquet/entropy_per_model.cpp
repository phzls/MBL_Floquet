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
 * Compute the entropy at a given time step and realization
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

	if (entropy_per_model_.size() != time_step){
		cout << "Not enough time steps are computed for entropy." << endl;
		cout << "Expected: " << num_realizations << endl;
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
	filename << name <<",Realizations=" << num_realizations << ",Total_time_step=" << time_step
		     << ",jump=" << jump << ",entropy_per_model.txt";

	if (output) cout << filename.str() <<endl;

	ofstream fout( filename.str().c_str() );

	for (int t=0; t<time_step; t++){
		double time = t * step_size * jump;
		generic_mean_sd(entropy_per_model_[t], mean[t], sd[t]);
		fout << setw(10) << time << setw(width) << mean[t] << setw(width) << sd[t] << endl;
	}
}