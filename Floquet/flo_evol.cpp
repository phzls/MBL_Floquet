#include <iostream>
#include "flo_evol.h"

using namespace std;

/**
 ** This file implements some simple functions in flo_evol.h
 **/


void FloEvolVanilla::Transition_Compute(TransitionMatrix& transition, 
const string& matrix_name) const{
	if (matrix_name == "Basic_Full"){
		if (constructed_) transition.Basic_Full(eigen[0].eigenvectors());
		else{
			cout << "Evolution Operator " << Repr() <<" has not been constructed." << endl;
			abort();
		}
	}
	else{
		cout << "Transition matrix " << matrix_name << " cannot be established for "
			 << Type() << endl;
		abort();
	}
}

void FloEvolParity::Transition_Compute(TransitionMatrix& transition, 
const string& matrix_name) const{
	if (matrix_name == "Parity_Full"){
		if (eigen_info_)
			// Use only one thread now
			transition.Parity_Full(eigen[0].eigenvectors(), eigen[1].eigenvectors(),1);
		else{
			cout << "Evolution Operator " << Repr() <<" has not been diagonalized." << endl;
			abort();
		}
	}
	else if (matrix_name == "Basic_Parity"){
		if (eigen_info_)
			// For now use 1 thread for parallelization
			transition.Basic_Parity(even_parity_, odd_parity_, 1);
		else{
			cout << "Evolution Operator " << Repr() <<" has not been diagonalized." << endl;
			abort();
		}
	}
	else if (matrix_name == "Basic_Parity_Full"){
		if (transition.Check_Matrix("Even_Full") != transition.Check_Matrix("Odd_Full")){
			cout << "Parity_Full construction is not consistent." << endl;
			abort();
		}
		if (!transition.Check_Matrix("Even_Full")){
			if (eigen_info_)
				// Use only one thread now
				transition.Parity_Full(eigen[0].eigenvectors(), eigen[1].eigenvectors(), 1);
			else{
				cout << "Evolution Operator " << Repr() <<" has not been diagonalized." << endl;
				abort();
			}
		}

		if (transition.Check_Matrix("Basic_Even") != transition.Check_Matrix("Basic_Odd")){
			cout << "Parity_Full construction is not consistent." << endl;
			abort();
		}
		if (!transition.Check_Matrix("Basic_Even")){
			if (eigen_info_)
				// For now use 1 thread for parallelization
				transition.Basic_Parity(even_parity_, odd_parity_, 1);
			else{
				cout << "Evolution Operator " << Repr() <<" has not been diagonalized." << endl;
				abort();
			}
		}
		transition.Basic_Parity_Full();
	}
	else if (matrix_name == "Basic_Full"){
		if (!transition.Check_Matrix("Basic_Parity_Full")) 
			Transition_Compute(transition, "Basic_Parity_Full");
		transition.Basic_Full_From_Basic_Parity_Full();
	}
	else{
		cout << "Transition matrix " << matrix_name << " cannot be established for "
			 << Type() << endl;
		abort();
	}
}