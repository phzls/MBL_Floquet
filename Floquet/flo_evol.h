#ifndef FLO_EVOL_H
#define FLO_EVOL_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "evol_class.h"

using namespace std;
using namespace Eigen;

/**
 ** This file contains base classes for Floquet time evolution operators
 **/

/*
 * The evolution operator for Floquet system with no apparent symmetry to reduce
 * time evolution operator. It uses ComplexEigenSolver<MatrixXcd> for the public
 * member eigen inherited from EvolMatrix. It only uses eigen.
 */

class FloEvolVanilla : public EvolMatrix< ComplexEigenSolver<MatrixXcd> >
{		
	protected:
		MatrixXcd evol_op_; // Time evolution operator
		bool constructed_; // Whether the matrix has been constructed and not erased

		stringstream repr_; // Representation string stream of the model
		string type_; // Type string of the model

	public:
		// When local dimension is not given
		FloEvolVanilla(int size): EvolMatrix< ComplexEigenSolver<MatrixXcd> >(size),
			constructed_(false){eigen.resize(1);}

		// When local dimension is given
		FloEvolVanilla(int size, int local_dim): 
			EvolMatrix< ComplexEigenSolver<MatrixXcd> >(size, local_dim), 
			constructed_(false){eigen.resize(1);}

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
			if (constructed_) eigen[0].compute(evol_op_, keep);
			else{
				cout << "The matrix for diagonalization does not exist." <<endl;
				abort();
			}
		}

		// Return the string format of representation string stream 
		string Repr() const {return repr_.str();} 

		// Return the type of the model
		string Type() const {return type_;} 

		// Erase the evolutionary operator
		void Evol_Erase() {evol_op_.resize(0,0); constructed_ = false;}

		virtual ~FloEvolVanilla() {};
};

/*
 * The evolution operator for Floquet system with parity symmetry to reduce
 * time evolution operator. It uses ComplexEigenSolver<MatrixXcd> for the public
 * member eigen inherited from EvolMatrix. It uses eigen_set of size 2, where
 * 0 is for even and 1 is for odd.
 */

class FloEvolParity : public EvolMatrix< ComplexEigenSolver<MatrixXcd> >
{		
	protected:
		MatrixXcd evol_op_even_; // Even time evolution operator
		MatrixXcd evol_op_odd_; // Odd time evolution operator
		bool constructed_; // Whether the matrix has been constructed and not erased

		stringstream repr_; // Representation string stream of the model
		string type_; // Type string of the model

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

		// Return the string format of representation string stream 
		string Repr() const {return repr_.str();} 

		// Return the type of the model
		string Type() const {return type_;} 

		// Erase the evolutionary operator
		void Evol_Erase() {evol_op_.resize(0,0); constructed_ = false;}

		virtual ~FloEvolVanilla() {};
};

#endif
