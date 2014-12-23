#include <fstream>
#include <iostream>
#include <complex>
#include "level_stats.h"

using namespace std;

void VanillaFloLevel::Data_Process(const EvolMatrix< ComplexEigenSolver<MatrixXcd> >& U){
	// The total dimension of Hilbert space
	int dim = U.GetDim();

	if (dim != U.eigen.eigenvalues().rows()){
		cout << "The dimension of eigenvalues of time evolution operator is not correct."<<endl;
		abort();
	}

	level_.resize(dim);

	vector<double> phases(dim); // Phases of all eigenvalues in ascending order
}