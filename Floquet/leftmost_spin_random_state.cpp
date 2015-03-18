//
// Created by Liangsheng Zhang on 3/18/15.
//

#include <iostream>
#include "randomc.h"
#include "initial_obj.h"
#include "eigen_output.h"
#include "constants.h"

using namespace std;
using namespace Eigen;

extern CRandomMersenne RanGen_mersenne; // points in [0,1)

/**
 ** This file creates initial density matrix of a mixed state, which is a random pure state on
 ** bloch sphere on the leftmost site, and an identity density matrix for the rest spins. The
 ** leftmost spin and the rest spins have no entanglement, so a tensor product of these two density
 ** matrices are taken. We treat 0 in the binary basis system as the rightmost position.
 **/


/*
 * This function gives the initial state in density matrix, and the evolution model must
 * be written in a complex matrix.
 */
void leftmost_spin_random_state(const InitInfo& init_info, MatrixXcd& init_state_density){
	init_state_density = MatrixXcd::Zero(init_info.dim, init_info.dim);

	const int rest_dim = init_info.dim / 2; // The dimension of the rest spins
	MatrixXcd rest_id = MatrixXcd::Identity(rest_dim, rest_dim);

	double theta = RanGen_mersenne.Random() * Pi;
	double phi = RanGen_mersenne.Random() * 2*Pi;

	for (int i=0; i<rest_dim; i++){
		for (int j=0; j<init_info.dim; j++){
			if ( j<rest_dim ){
				init_state_density(i,j) = (1+cos(theta)) / double(init_info.dim) * rest_id(i,j);
			}
			else{
				init_state_density(i,j) = (1+sin(theta)) / double(init_info.dim) *
						exp(-Complex_I * phi) * rest_id(i, j-rest_dim);
			}
		}
	}

	for (int i=rest_dim; i<init_info.dim; i++){
		for (int j=0; j<init_info.dim; j++){
			if ( j<rest_dim ){
				init_state_density(i,j) = conj(init_state_density(j,i));
			}
			else{
				init_state_density(i,j) = (1-cos(theta)) / double(init_info.dim) *
						rest_id(i-rest_dim, j-rest_dim);
			}
		}
	}

	if ( abs(init_state_density.trace() - 1.0) < 1.0e-8 ){
		cout << "Initial density matrix has trace not equal to 1." << endl;
		cout << "Obtained trace: " << init_state_density.trace() << endl;
		abort();
	}

	if (init_info.debug){
		cout << "Initial density matrix random angle: " << endl;
		cout << "theta: " << theta << endl;
		cout << "phi: " << phi << endl;
		cout << endl;
		cout << "Initial state in density matrix:" << endl;
		complex_matrix_write(init_state_density);
		cout << endl;
	}
}

