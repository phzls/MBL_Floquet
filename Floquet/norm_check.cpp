#include <complex>
#include <Eigen/Dense>
#include <iostream>

using namespace std;
using namespace Eigen;

/**
 ** This function checks whether total norm of a complex vector is close 1. delta gives the small
 ** number for individual component error tolerance. state_type gives the name of the state.
 **/

void norm_check(const VectorXcd& state, double delta, const string& state_type){
	double state_norm = 0;
	for (int i=0; i< state.size(); i++){
		state_norm += norm(state(i));
	}

	cout << state_type << " norm: " << state_norm << endl;

	// Consider delta to be the error tolerance for the individual component
	if (abs(state_norm - 1) > (delta * state.size())){
		cout << state_type << " norm is not 1." << endl;
		abort();
	}
}