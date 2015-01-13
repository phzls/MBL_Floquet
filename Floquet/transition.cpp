#include <iostream>
#include <cmath>
#include <complex>
#include "transition.h"

using namespace std;

/**
 ** Implement simple functions in transition.h
 **/

const MatrixXcd& TransitionMatrix::Matrix(const string& matrix_name) const{
	map<string, MatrixXcd*>::const_iterator it;
	it = constructed_type_.find(matrix_name);

	if (it == constructed_type_.end()){
		cout << "Requested transition matrix has not been constructed." << endl;
		abort();
	}
	else return (*(it -> second));
}

bool TransitionMatrix::Check_Matrix(const string& matrix_name) const{
	if (constructed_type_.find(matrix_name) == constructed_type_.end()) return false;
	else return true;
}

void TransitionMatrix::Erase_Matrix(const string& matrix_name){
	map<string, MatrixXcd*>::iterator it;
	it = constructed_type_.find(matrix_name);

	if (it != constructed_type_.end()){
		(*(it -> second)).resize(0,0);
		constructed_type_.erase(it);
	}
}

void TransitionMatrix::Erase_All(){
	map<string, MatrixXcd*>::iterator it;
	for (it = constructed_type_.begin(); it != constructed_type_.end(); it ++){
		(*(it -> second)).resize(0,0);
	}
	constructed_type_.clear();
}

void TransitionMatrix::Print(const string& matrix_name) const {
	map<string, MatrixXcd*>::const_iterator it = constructed_type_.find(matrix_name);
	if (it == constructed_type_.end()) 
		cout << "Requested transition matrix " << matrix_name << " not constructed." <<endl;
	else{
		cout << "Transition Matrix: " << it -> first << endl;
		for (int i=0; i< (*(it -> second)).rows(); i++){
			for (int j=0; j< (*(it -> second)).cols(); j++){
				cout << real((*(it -> second))(i,j));
				double im = imag((*(it -> second))(i,j));
				if (abs(im) > 1.0e-15){
					if (im >= 0) cout <<"+";
					cout << im << "j";
				}
				cout << "    ";
			}
			cout << endl;
		}
	}
}

void TransitionMatrix::Print_All() const {
	map<string, MatrixXcd*>::const_iterator it;
	for (it = constructed_type_.begin(); it != constructed_type_.end(); it ++){
		Print(it -> first);
	}
}
