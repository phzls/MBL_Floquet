#ifndef TRANSITION_H
#define TRANSITION_H

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <map>

using namespace std;
using namespace Eigen;

/**
 ** This file defines the class for basis transition matrix among different bases. All 
 ** eigenvectors used are assumed to be stored column-wise, i.e., each column is an eigenvector.
 ** Note for general diagonalization, if the matrix has almost-degenerate eigenvalues, then there
 ** is no guarantee that the eigenvectors of orthogornal to each other. This will make the inverse
 ** transition matrix not the hermitian conjugate of the forward transtion matrix.
 **/

class TransitionMatrix
{
	private:
		// Transform an even vector written in even parity states to  binary basis. 
		// The inverse transform for an even vector is given by its adjoint
		MatrixXcd basic_even_;
		// Transform an even vector written in odd parity states to  binary basis. 
		// The inverse transform for an odd vector is given by its adjoint
		MatrixXcd basic_odd_;

		// Transform an even vector written in even full chain eigenstates to even parity states
		MatrixXcd even_full_;
		// Transform an odd vector written in odd full chain eigenstates to odd parity states
		MatrixXcd odd_full_;

		// Transform the even part of a vector written in half chain eigenstates basis to 
		// even parity states
		MatrixXcd even_half_;
		// Transform the odd part of a vector written in half chain eigenstates basis to 
		// odd parity states
		MatrixXcd odd_half_;

		/*
		 * Currently the above two matrices cannot be initialized.
		 */


		// Transform the full chain eigenstates directly to binary basis states
		MatrixXcd basic_full_;

		// Transform an even vector written in even full chain eigenstates to binary basis
		MatrixXcd basic_even_full_;
		// Transform an even vector written in even full chain eigenstates to binary basis
		MatrixXcd basic_odd_full_;

		// Store constructed matrices by its names
		map<string, MatrixXcd*> constructed_type_;

	public:
		// Return constant reference to corresponding matrix. If not constructed, abort
		// The string specifies the type of matrix
		const MatrixXcd& Matrix(const string&) const;

		// Check wether corresponding matrix has been constructed. The string specifies the
		// matrix type
		bool Check_Matrix(const string&) const; 

		// Erase the matrix content if it has been constructed.
		void Erase_Matrix(const string&);

		// Erase all matrix content
		void Erase_All();

		// Print out specific matrix
		void Print(const string&) const;

		// Print out all matrices
		void Print_All() const;

		// Construct above basci_parity matrices. The first vector is even_parity; the second 
		// is odd_parity. The last integer is the number of threads used for parallelization
		void Basic_Parity(const vector<vector<int> >&, const vector<vector<int> >&, int);

		// Construct above parity_full matrices. The first matrix is even_evec; the second 
		// is odd_evec. The last integer is the number of threads used for parallelization
		void Parity_Full(const MatrixXd&, const MatrixXd&, int);
		void Parity_Full(const MatrixXcd&, const MatrixXcd&, int);

		// Use basic_parity and parity_full to construct basic_parity_full
		void Basic_Parity_Full();

		// Directly use full chain eigenstates to construct basic_full. The eigenstates have
		// no parity
		void Basic_Full(const MatrixXd& evec);
		void Basic_Full(const MatrixXcd& evec);
};

#endif