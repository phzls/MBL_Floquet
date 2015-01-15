#include <vector>
#include <cmath>
#include <complex>
#include "constants.h"
#include "initial_obj.h"

using namespace std;
using namespace Eigen;

MTRand u1rand(time(NULL));

/**
 ** This function creates an initial state which is the random product state. Specifically,
 ** the 1/2 spin pointing random direction on Bloch sphere at each site.
 **/

void random_product(const InitInfo& init_info, const TransitionMatrix transition,
VectorXcd& init_state){
	const int size = init_info.size; // System size
	const int total_rank = init_state.size(); // Total dimension of Hilbert space
	const double delta = init_info.norm_delta; // A small quantity

	if (total_rank != (1<<size)){
		cout << "Total rank for initial state is wrong." << endl;
		cout << "Expected total rank: " << (1<<size) << endl;
		cout << "Obtained total rank: " << total_rank << endl;
		abort();
	}

	VectorXcd& init_basic(total_rank); // The initial state written in basic binary basis

	// Record the amplitude of up and down spin at each site
	vector<vector<complex<double> > > spin_amp(size); 

	for (int i=0;i<size;i++)
	{
		double theta = u1rand() * Pi;
		double phi = u1rand() * 2*Pi;

		spin_amp.resize(2); // 0 for spin down, 1 for spin up

		spin_amp[i][1] = complex<double>(cos(theta/2),0);
		spin_amp[i][0] = exp( I*complex<double>(phi,0) ) * complex<double>(sin(theta/2),0);
	}

	// Calculate the amplitude of the state in binary basis
	for (int i=0;i<total_rank;i++)
	{
		complex<double> product = complex<double>(1,0);
		int site = i;

		for (int j=0;j<size;j++)
		{
			int spin = site & 1;
			site = site >> 1;

			product = product * spin_amp[j][spin];
		}

		init_basic(i) = product;
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