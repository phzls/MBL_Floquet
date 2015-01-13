#include <iostream>
#include "transition.h"

using namespace std;

/**
 ** Implement simple functions in transition.h
 **/

const MatrixXcd& TransitionMatrix::Matrix(const string& matrix_name) const{
	map<string, MatrixXcd&>::iterator it;
	it = constructed_type_.find(matrix_name);

	if (it == constructed_type_.end()){
		cout << "Requested transition matrix has not been constructed." << endl;
		abort();
	}
	else return it -> second;
}

bool TransitionMatrix::Check_Matrix(const string& matrix_name) const{
	if (constructed_type_.find(matrix_name) == constructed_type_.end()) return false;
	else return true;
}

void TransitionMatrix::Erase_Matrix(const string& matrix_type){
	map<string, MatrixXcd&>::iterator it;
	it = constructed_type_.find(matrix_name);

	if (it != constructed_type_.end()){
		it -> second.resize(0,0);
		constructed_type_.erase(it);
	}
}

void TransitionMatrix::Erase_All(){
	map<string, MatrixXcd&>::iterator it;
	for (it = constructed_type_.begin(); it != constructed_type_.end(); it ++){
		it -> second.resize(0,0);
	}
	constructed_type_.clear();
}
