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
		int model_id_; // unique ID of this model
		const bool iso_keep_; // Whether keep another EvolMatrix pointer which points
							  // to the isolated part of a model coupled to the bath
		bool eigen_computed_; // Whether eigensystems of this model have been solved
		
		// Pointers pointing to isolated part of a model coupling to the bath
		EvolMatrix<ComplexEigenSolver<MatrixXcd> >* model_iso_complex_;
		EvolMatrix<EigenSolver<MatrixXd> >* model_iso_real_;

	public:
		// When local dimension is not given
		EvolMatrix(int size, bool iso_keep = false): 
			size_(size), local_dim_(2), dim_(1 << size), iso_keep_(iso_keep), 
			model_iso_complex_(NULL), model_iso_real_(NULL),eigen_computed_(false)
			{model_num++; model_id_ = model_num;}

		// When local dimension is explicitly given
		EvolMatrix(int size, int local_dim, bool iso_keep = false):
			size_(size), local_dim_(local_dim), dim_(int(pow(double(local_dim), double(size)))),
			iso_keep_(iso_keep), model_iso_complex_(NULL), model_iso_real_(NULL),
			eigen_computed_(false) {model_num++; model_id_ = model_num;}

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

		// Get the pointer to the complex isolated matrix
		void Get_Iso(EvolMatrix<ComplexEigenSolver<MatrixXcd> >*& model)
		{model = model_iso_complex_;}

		// Get the pointer to the real isolated matrix
		void Get_Iso(EvolMatrix<EigenSolver<MatrixXd> >*& model) const {model = model_iso_real_;}

		// Check whether eigensystems have been solved
		bool Get_Eigen_Computed() const {return eigen_computed_;}

		static int model_num;

		int Model_ID() const {return model_id_;}

		virtual ~EvolMatrix(){
			if (model_iso_complex_ != NULL) delete model_iso_complex_;
			if (model_iso_real_ != NULL) delete model_iso_real_;
		};
};

template <class T>
int EvolMatrix<T>::model_num = 0;

#endif