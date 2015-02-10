#include <utility>
#include <vector>
#include <iostream>
#include <complex>
#include <Eigen/Core>
#include "generic_func.h"

using namespace std;
using namespace Eigen;

/*
 * This file computes average
 */

template <class T1, class T2>
void leftmost_sigma_z_sum(Matrix<T1, Dynamic, Dynamic>& m, const vector< vector<T2> >& basis, 
const string& type){
	// dimension of space defined by the bases
	const int dim = basis.size();

	// Check the size of basis vectors
	for (int i=0; i<dim; i++){
		if ( basis[i].size() != dim ){
			cout << "The size of basis vectors is not correct." << endl;
			abort();
		}
	}

	// Check the size of m
	if (m.size() == 0){
		m.resize(dim, dim);

		for (int i=0; i<dim; i++){
			for (int j=0; j<dim; j++){
				m(i,j) = 0;
			}
		}
	}
	else{
		if (m.rows() != dim || m.cols() != dim){
			cout << "In rightmost sigma z construction, the size of matrix does not match." 
				 << endl;
			abort();
		}
	}

	const int down = dim / 2; // Below this number the binary representation of spins has
							  // leftmost spin down

	// Compute the entry
	for (int i=0; i<dim; i++){
		for (int j=0; j<dim; j++){
			T2 sum = 0;
			for (int k=0; k<dim; k++){
				if (k < down){ // Leftmost bit is down
					sum -= generic_conj(basis[i][k]) * basis[j][k];
				}
				else{ // Leftmost bit is up
					sum += generic_conj(basis[i][k]) * basis[j][k];
				}
			}

			T1 element;

			if (type == "entry" || type == "Entry"){
				update(element, sum);
				m(i,j) += element;
			}
			else if (type == "norm" || type == "Norm"){
				m(i,j) += generic_norm(sum);
			}
			else{
				cout << "Type of representation does not exist." << endl;
				abort();
			}
		}
	}
}