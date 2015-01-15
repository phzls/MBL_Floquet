#include <complex>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;s

/**
 ** This function checks whether total norm of a complex vector is close 1. delta gives the small
 ** number. state_type gives the name of the state.
 **/

void norm_check(const VectorXcd& state, double delta, const string& state_type){
	double state_norm = 0;
	for (int i=0; i< state.size(); i++){
		state_norm += norm(state(i));
	}

	cout << state_type << " norm: " << state_norm << endl;
	if (abs(state_norm - 1) > delta){
		cout << state_type << " norm is not 1." << endl;
		abort();
	}s
}