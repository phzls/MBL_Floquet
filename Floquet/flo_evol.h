#ifndef FLO_EVOL_H
#define FLO_EVOL_H

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include "evol_class.h"

using namespace std;
using namespace Eigen;

/**
This file contains base classes for Floquet time evolution operators
**/

/*
The evolution operator for Floquet system with no apparent symmetry to reduce
time evolution operator. It uses ComplexEigenSolver<MatrixXcd> for the public
member eigen inherited from EvolMatrix.
*/

class FloEvolVanilla : public EvolMatrix< ComplexEigenSolver<MatrixXcd> >
{		
	protected:
		MatrixXcd evol_op_; // Time evolution operator

	public:
		// When local dimension is not given
		FloEvolVanilla(int size): EvolMatrix(size){
			evol_op_ = MatrixXcd::Zero(dim_, dim_);
		}

		// When local dimension is given
		FloEvolVanilla(int size, int local_dim): EvolMatrix(size, local_dim){
			evol_op_ = MatrixXcd::Zero(dim_, dim_);
		}

		// Diagnolize time evolution matrix with eigenvectors kept
		void Evol_Diag() {eigen.compute(evol_op_);} 

		// Diagnolize time evolution matrix, user can determine whether eigenvectors are kept
		// False is not kept; True is kept
		void Evol_Diag(bool keep) {eigen.compute(evol_op_, keep);} 

		virtual ~FloEvolVanilla() {};
};

#endif
