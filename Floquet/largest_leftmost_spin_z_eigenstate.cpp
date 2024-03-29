#include <iostream>
#include <utility>
#include "initial_obj.h"
#include "eigen_output.h"
#include "sort.h"

using namespace std;

/**
 ** This file creates initial state, which is an eigenstate of a particular evolution model.
 ** This eigenstate has the largest value for the leftmost spin z, where spin up is 1
 ** and spin down is -1. We treat 0 in the binary basis system as the rightmost position.
 **/


/*
 * This function gives the initial state in density matrix, and the evolution model must
 * be written in a complex matrix.
 */
void largest_leftmost_spin_z_complex_eigenstate(const InitInfo& init_info, 
MatrixXcd& init_state_density){

	// Any number smaller than this has leftmost spin down
	const int down = 1 << (init_info.size - 1);

	const vector<const ComplexEigenSolver<MatrixXcd>* > eigen = init_info.complex_eigen;

	if (init_info.debug){
		cout << "Eigenvectors and eigenvalues of init_model:" << endl;
			
		for (int i=0; i<eigen.size(); i++){
			cout << "Sector " << i <<" :" << endl;
			cout << "Eigenvectors:" << endl;
			complex_matrix_write(eigen[i] -> eigenvectors());
			cout << endl;
			cout << "Eigenvalues:" << endl;
			complex_matrix_write(eigen[i] -> eigenvalues());
			cout << endl;
		}
	}

	// Record magnitude of leftmost spin z. The first integer in the pair is for sector,
	// the second is for the position in that sector
	vector<pair<double, pair<int, int> > > leftmost_spin_z_pos(init_info.dim); 

	int index = 0; // A continuous index of all eigenstates
	for (int i=0; i<eigen.size();i++){
		for (int j=0; j<eigen[i] -> eigenvectors().cols(); j++){
			leftmost_spin_z_pos[index].first = 0;
			leftmost_spin_z_pos[index].second.first = i;
			leftmost_spin_z_pos[index].second.second = j;

			double temp = 0;

			for (int k=0; k<eigen[i] -> eigenvectors().rows(); k++){
				if (k<down) temp -= norm( eigen[i] -> eigenvectors()(k,j) );
				else temp += norm( eigen[i] -> eigenvectors()(k,j) );
			}

			leftmost_spin_z_pos[index].first = temp;
			index ++;
		}
	}

	if (index != leftmost_spin_z_pos.size()){
		cout << "Number of eigenstates is not correct." << endl;
		cout << "Expected Number: " << leftmost_spin_z_pos.size() << endl;
		cout << "Obtained Number: " << index << endl;
		abort();
	}

	// Sort according to leftmost spin z value
	sort(leftmost_spin_z_pos.begin(), leftmost_spin_z_pos.end(), 
		Vec_Pair_Double_First_Sort<pair<int,int> >);

	// Sector of the eigenstate with the largest leftmost spin z value
	const int sec = leftmost_spin_z_pos[index-1].second.first;
	// Relative position in the sector for the eigenstate with the largest leftmost 
	// spin z value
	const int rel_pos = leftmost_spin_z_pos[index-1].second.second;

	// The eigenstate with the largest leftmost spin z value
	VectorXcd state(index);
	for (int i=0; i<index;i++){
		state(i) = eigen[sec] -> eigenvectors()(i,rel_pos);
	}

	if (init_info.debug){
		cout << "Largest leftmost spin z value:  " << leftmost_spin_z_pos[index-1].first << endl;
		cout << endl;
		cout << "Initial state in binary basis:" << endl;
		complex_matrix_write(state);
		cout << endl;
	}

	state_to_density(state, init_state_density);

	if (init_info.debug){
		cout << "Initial state in density matrix:" << endl;
		complex_matrix_write(init_state_density);
		cout << endl;
	}
}