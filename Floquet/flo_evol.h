#ifndef FLO_EVOL_H
#define FLO_EVOL_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include "evol_class.h"

using namespace std;
using namespace Eigen;

/**
 ** This file contains base classes for Floquet time evolution operators
 **/

/*
 * The evolution operator for Floquet system with no apparent symmetry to reduce
 * time evolution operator. It uses ComplexEigenSolver<MatrixXcd> for the public
 * member eigen inherited from EvolMatrix.
 */

class FloEvolVanilla : public EvolMatrix< ComplexEigenSolver<MatrixXcd> >
{		
	protected:
		MatrixXcd evol_op_; // Time evolution operator
		bool constructed_; // Whether the matrix has been constructed and not erased

	public:
		// When local dimension is not given
		FloEvolVanilla(int size): EvolMatrix< ComplexEigenSolver<MatrixXcd> >(size),
			constructed_(false){}

		// When local dimension is given
		FloEvolVanilla(int size, int local_dim): 
			EvolMatrix< ComplexEigenSolver<MatrixXcd> >(size, local_dim), 
			constructed_(false){}

		// Diagnolize time evolution matrix with eigenvectors kept
		void Evol_Diag() {
			if (constructed_) eigen.compute(evol_op_);
			else{
				cout << "The matrix for diagonalization does not exist." <<endl;
				abort();
			}
		} 

		// Diagnolize time evolution matrix, user can determine whether eigenvectors are kept
		// False is not kept; True is kept
		void Evol_Diag(bool keep) {
			if (constructed_) eigen.compute(evol_op_, keep);
			else{
				cout << "The matrix for diagonalization does not exist." <<endl;
				abort();
			}
		} 

		virtual ~FloEvolVanilla() {};
};

#endif
