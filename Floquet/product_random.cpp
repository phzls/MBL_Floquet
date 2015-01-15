#include <vector>
#include <cmath>
#include <complex>
#include "constants.h"
#include "initial_obj.h"

using namespace std;
using namespace Eigen;

MTRand u1rand(time(NULL));

/**
 ** Construct a state which is a product of two random pure states of half chains 
 **/

void product_random(const InitInfo& init_info, const TransitionMatrix transition,
VectorXcd& init_state){
	const int size = init_info.size; // System size
	const int total_rank = init_state.size(); // Total dimension of Hilbert space
	const double delta = init_info.norm_delta; // A small quantity

	if (total_rank != (1<<size)){
		cout << "Total rank for initial state is wrong." << endl;
		cout << "Expected total rank: " << (1<<size) << endl;
		cout << "Obtained total rank: " << total_rank << endl;
		abort();
	

	if (size%2 != 0){
		cout << "product_random states are only constructed for chain of even number of spins."
			 << endl;
		cout << "Spin chain size: " << size << endl;
		abort();
	}

	const int half_rank = 1 << (size/2); // Dimension of the half-chain space

	// Store amplitudes for left and right random pure state
	vector<complex<double> > amplitude_left(half_rank); 
	vector<complex<double> > amplitude_right(half_rank);

	// Construct amplitudes
	random_pure_amplitude(amplitude_left);
	random_pure_amplitude(amplitude_right);

	VectorXcd& init_basic(total_rank); // The initial state written in basic binary basis

	for (int n=0;n<half_rank;n++)
	{
		for (int m=0;m<half_rank;m++)
		{
			int pos = n*half_rank + m;

			init_basic(pos) = amplitude_left[n] * amplitude_right[m];
		}
	}

	norm_check(init_basic, delta, "Binary state");

	init_state = transition.Matrix("Basic_Full").adjoint() * init_basic;

	if (init_state.size() != total_rank){
		cout << "init_state size is wrong." << endl;
		cout << "Expected size: " << total_rank << endl;
		cout << "Obtained size: " << init_state.size() << endl;
		abort();
	}

	norm_check(init_state, delta, "Initial state"); 
}