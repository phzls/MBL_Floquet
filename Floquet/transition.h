#ifndef TRANSITION_H
#define TRANSITION_H

#include <iostream>
#include <Eigen/Dense>
#include <map>

using namespace std;
using namespace Eigen;

/**
 ** This file defines the class for basis transition matrix among different bases.
 **/

class TransitionMatrix
{
	private:
		// Transform an even vector written in even parity states to  binary basis. 
		// The inverse transform for an even vector is given by its adjoint
		MatrixXcd basic_even_(0,0);
		// Transform an even vector written in odd parity states to  binary basis. 
		// The inverse transform for an odd vector is given by its adjoint
		MatrixXcd basic_odd_(0,0);

		// Construct above two matrices given...
		void Basic_Parity_();

		// Transform an even vector written in even full chain eigenstates to even parity states
		MatrixXcd even_full_(0,0);
		// Transform an odd vector written in odd full chain eigenstates to odd parity states
		MatrixXcd odd_full_(0,0);


		// Transform the even part of a vector written in half chain eigenstates basis to 
		// even parity states
		MatrixXcd even_half_(0,0);
		// Transform the odd part of a vector written in half chain eigenstates basis to 
		// odd parity states
		MatrixXcd odd_half_(0,0);

		// Transform the full chain eigenstates directly to binary basis states
		MatrixXcd basic_full_(0,0);

		// Transform an even vector written in even full chain eigenstates to binary basis
		MatrixXcd basic_full_even_(0,0);
		// Transform an even vector written in even full chain eigenstates to binary basis
		MatrixXcd basic_full_odd_(0,0);

		// Store constructed matrices by its names
		map<string, MatrixXcd&> constructed_type_;

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
}

#endif