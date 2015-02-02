#include <vector>
#include <cmath>
#include <complex>
#include "constants.h"
#include "initial_obj.h"
#include "eigen_output.h"

using namespace std;
using namespace Eigen;

void random_pure(const InitInfo& init_info, MatrixXcd& init_state_density){
	const int size = init_info.size; // System size
	const int total_rank = (1<<size); // Total dimension of Hilbert space
	const double delta = init_info.norm_delta; // A small quantity

	// The initial state written in basic binary basis
	VectorXcd init_basic = VectorXcd::Zero(total_rank);

	// Record the amplitude for each state
	vector<complex<double> > amp(total_rank); 

	// Construct amplitudes
	random_pure_amplitude(amp);

	for (int i=0; i<total_rank;i++) init_basic(i) = amp[i];

	norm_check(init_basic, delta, "Binary state");

	if (init_info.debug){
		cout << "Initial state in binary basis:" << endl;
		complex_matrix_write(init_basic);
		cout << endl;
	}

	state_to_density(init_basic, init_state_density);

	if (init_info.debug){
		cout << "Initial state in density matrix:" << endl;
		complex_matrix_write(init_state_density);
		cout << endl;
	}
}

