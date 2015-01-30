#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

/**
 ** This function converts a complex state vector to a complex density matrix
 **/

void state_to_density(const VectorXcd& state, MatrixXcd& matrix){
	matrix = MatrixXcd::Zero(state.size(), state.size());

	for (int i=0; i<state.size(); i++){
		for (int j=i; j< state.size(); j++){
			matrix(i,j) = state(i) * conj(state(j));
			if (i!=j) matrix(j,i) = conj(matrix(i,j));
		}
	}
}