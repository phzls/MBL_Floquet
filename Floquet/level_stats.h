#ifndef LEVEL_STATS_H
#define LEVEL_STATS_H

#include <sstream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "results.h"
#include "evol_class.h"

using namespace std;
using namespace Eigen;

/**
This file includes output classes which output level statistics from Hamiltonian/unitary
time evolution operator
**/

/*
This class outputs level statistics for vanilla Floquet time evolution operators
*/
class VanillaFloLevel: public ResultsOutput< EvolMatrix< ComplexEigenSolver<MatrixXcd> > >
{
	private:
		vector<double> level_; // The vector that holds all the level spacings
		stringstream filename_; // filename
		const bool spacing_; // If it is true, level spacing will be outputted
		const bool mean_; // If it is true, the mean of level spacings will be outputted
		const bool mean_square_; // If it is true, the mean of square will be outputted

	public:
		VanillaFloLevel(bool spacing, bool mean, bool mean_square):
			spacing_(spacing), mean_(mean), mean_square_(mean_square) {};
			
		void Data_Process(const EvolMatrix< ComplexEigenSolver<MatrixXcd> >& );
		void Data_Output();

		~VanillaFloLevel() {};
};

#endif