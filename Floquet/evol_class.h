#ifndef EVOL_CLASS_H
#define EVOL_CLASS_H

#include <cmath>
#include <string>
#include <vector>
#include "transition.h"

using namespace std;

/**
 ** A base class defines the evolution of a quantum system. The size of the system needs to
 ** be passed in. In this case the local dimension at each site is assumed to be 2, and the 
 ** total dimension is calculated. The local dimension can also be specified when 
 ** constructing the class. The template type of the class determines the type of public
 ** member eigen used for eigenvalues and eigenvectors.
 **/

template <class T>
class EvolMatrix
{
	protected:
		const int size_; // Size of the system
		const int dim_; // Dimension of the space given size
		const int local_dim_; // Local dimension of the Hilber space

	public:
		// When local dimension is not given
		EvolMatrix(int size): 
			size_(size), local_dim_(2), dim_(1 << size){}

		// When local dimension is explicitly given
		EvolMatrix(int size, int local_dim):
			size_(size), local_dim_(local_dim), dim_(int(pow(double(local_dim), double(size)))){}

		//T eigen; // Eigenvalues and Eigenvectors

		vector<T> eigen; // Eigenvalues and Eigenvectors for possible different sectors if 
						 // there is some symmetry to simplity the operator

		vector<string> eigen_name; // Name of each eigensectors

		// Initialize parameters which will be used in constructing time evolution operator.
		// These parameters will only exist for concrete model classes.
		virtual void Evol_Para_Init() = 0;

		// Constructing time evolution matrix
		virtual void Evol_Construct() = 0; 
		
		// Diagnolize time evolution matrix with eigenvectors kept
		virtual void Evol_Diag() = 0;

		// Diagnolize time evolution matrix, user can determine whether eigenvectors are kept
		// False is not kept; True is kept
		virtual void Evol_Diag(bool keep) = 0;

		// Compute various transition matrix. The string specifies what transition to compute
		virtual void Transition_Compute(TransitionMatrix&, const string&) const = 0;

		// Erase the matrix to free some memroy
		virtual void Evol_Erase() = 0;

		// Return the string format of representation string stream 
		virtual string Repr() const = 0;

		// Return the type of the model as a string, i.e., representation without
		// concerete parameters
		virtual string Type() const = 0;

		// Return the type of basis which eigenstates are written in
		virtual string Eigen_Type() const = 0;

		// Return the size of the system
		int Get_Size() const {return size_;}

		// Return the dimension of total Hilbert space
		int Get_Dim() const {return dim_;} 

		// Return the dimension of each sector of symmetry
		virtual vector<int> Get_Sector_Dim() const = 0;

		// Return the complex unitary operator according to a string
		virtual const MatrixXcd& Get_U(string) const = 0;

		// Return the complex unitary operator according to an integer
		virtual const MatrixXcd& Get_U(int) const = 0;

		// Return the real Hamiltonian operator according to a string
		virtual const MatrixXd& Get_real_H(string) const = 0;

		// Return the real Hamiltonian operator according to an integer
		virtual const MatrixXd& Get_real_H(int) const = 0;

		virtual ~EvolMatrix(){};
};

#endif