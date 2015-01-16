#include <cmath>
#include <Eigen/Eigenvalues>
#include "evol_data.h"
#include "methods.h"

using namespace std;
using namespace Eigen;

/**
 ** This file contains functions related to entropy_per_model. It is the entropy of the 
 ** left part of a spin chain model whose local dimension is 2. The entropy is computed in bits.
 **/

/*
 * Initialize the entropy
 */
void EvolData::Entropy_Per_Model_Init_(const AllPara& parameters){
	const int num_realizations = parameters.generic.num_realizations;
	const int time_step = parameters.evolution.time_step;

	entropy_per_model_.resize(num_realizations);

	for (int i=0; i<num_realizations; i++){
		entropy_per_model_[i].resize(time_step);

		for (int j=0; j<time_step; j++){
			entropy_per_model_[i][j] = 0;
		}
	}
}

/*
 * Compute the entropy at a given time step and realization
 */
void EvolData::Entropy_Per_Model_Cal_(const VectorXcd& state_basic, const StepInfo& info){
	const int realization = info.realization;
	const int time_step = info.time_step;
	const int left_size = info.left_size;

	MatrixXcd reduced_density; // Reduced density matrix for the left part
	reduced_density_left_2(state_basic, size_, left_size, reduced_density);

	SelfAdjointEigenSolver<MatrixXcd> density_eigen; // Eigen for reduced density matrix

	density_eigen.compute(reduced_density, false); // Eigenvectors not computed

	entropy_per_model_[realization][left_size] = 0;

	for (int i=0; i<density_eigen.eigenvalues().rows();i++){
		double eval = density_eigen.eigenvalues()(i);

		if (abs(eval)>1.0e-15)
		{
			entropy_per_model_[realization][left_size] += -eval*log2(eval);
		}
	}
}