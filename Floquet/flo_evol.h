#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Eigenvalues> 
#include "evol_class.h"

using namespace std;
using namespace Eigen;

/**
This file contains base classes for Floquet time evolution operators
**/

/*
The evolution operator for Floquet system with no apparent symmetry to reduce
time evolution operator.
*/

class FloEvolVanilla : public EvolMatrix
{		
	protected:
		MatrixXcd evol_op_; // Time evolution operator

	public:
		ComplexEigenSolver<MatrixXcd> eigen; // Eigenvalues and Eigenvectors

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

