#include <complex>
#include "rightmost_sigma_z_sum.tpp"

using namespace std;

/**
 ** This file gives template specializations of function update used in rightmost_sigma_z_sum.tpp
 ** which is used to circumvent the type check at compile time in if control. It is used when
 ** entry update is needed.
 **/

template <>
void update(double& val1, complex<double> val2) {
	cout << "Type error in updating matrix." <<endl;
	abort();
}

template <>
void update(double& val1, complex<int> val2) { 
	cout << "Type error in updating matrix." <<endl;
	abort();
}

template <>
void update(int& val1, complex<double> val2) {
	cout << "Type error in updating matrix." <<endl;
	abort();
}

template <>
void update(int& val1, complex<int> val2) {
	cout << "Type error in updating matrix." <<endl;
	abort();
}