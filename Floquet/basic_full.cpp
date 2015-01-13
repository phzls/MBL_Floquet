#include <iostream>
#include <complex>
#include "transition.h"

using namespace std;
using namespace Eigen;

/**
 ** This file implements construction of transition matrix from full chain eigenstates to
 ** binary basis
 **/

/*
 * This function uses Basic_Parity and Parity_Full for construction
 */
void TransitionMatrix::Basic_Parity_Full(){
	if (Check_Matrix("Basic_Even_Full")){
		cout << "Transition matrix basic_even_full has been constructed." << endl;
		abort();
	}

	if (Check_Matrix("Basic_Odd_Full")){
		cout << "Transition matrix basic_odd_full has been constructed." << endl;
		abort();
	}

	if (!Check_Matrix("Basic_Even")){
		cout << "For basic_even_full, transition matrix basic_even has not been constructed." 
		     << endl;
		abort();
	}

	if (!Check_Matrix("Basic_Odd")){
		cout << "For basic_odd_full, transition matrix basic_odd has not been constructed." 
		     << endl;
		abort();
	}

	if (!Check_Matrix("Even_Full")){
		cout << "For basic_even_full, transition matrix even_full has not been constructed." 
			 << endl;
		abort();
	}

	if (!Check_Matrix("Odd_Full")){
		cout << "For basic_odd_full, transition matrix odd_full has not been constructed." 
			 << endl;
		abort();
	}

	basic_even_full_ = basic_even_ * even_full_;
	basic_odd_full_ = basic_odd_ * odd_full_;

	constructed_map_["Basic_Even_Full"] = basic_even_full_;
	constructed_map_["Basic_Odd_Full"] = basic_odd_full_;

}

/*
 * This function directly constructs the transition matrix from full to basic. The full chain
 * eigenstates have no parity, and are stored column-wise.
 */
void Basic_Full(const MatrixXd& evec){
 	const int total_rank = evec.cols(); // Total dimension

 	if (Check_Matrix("Basic_Full")){
		cout << "Transition matrix basic_full has been constructed." << endl;
		abort();
	}

	// Initialize the matrix
	basic_full_.resize(total_rank, total_rank);

	for (int i=0; i< total_rank; i++){
		for (int j=0; j<total_rank; j++){
			basic_full_(i,j) = complex<double>(evec(i,j), 0);
		}
	}

	constructed_map_["Basic_Full"] = basic_full_;
}

void Basic_Full(const MatrixXcd& evec){
 	const int total_rank = evec.cols(); // Total dimension

 	if (Check_Matrix("Basic_Full")){
		cout << "Transition matrix basic_full has been constructed." << endl;
		abort();
	}

	// Initialize the matrix
	basic_full_.resize(total_rank, total_rank);

	for (int i=0; i< total_rank; i++){
		for (int j=0; j<total_rank; j++){
			basic_full_(i,j) = evec(i,j);
		}
	}

	constructed_map_["Basic_Full"] = basic_full_;
}