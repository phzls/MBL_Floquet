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
This class outputs 
*/
class VanillaFloLevel: public ResultsOutput< EvolMatrix< ComplexEigenSolver<MatrixXcd> > >
{
	private:
		vector<double> level_;
		stringstream filename_;

	public:
		void Data_Process(const EvolMatrix< ComplexEigenSolver<MatrixXcd> >& );
		void Data_Output();

		~VanillaFloLevel() {};
};

#endif