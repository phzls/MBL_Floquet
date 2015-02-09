#ifndef FLO_EVOL_H
#define FLO_EVOL_H

#include <iostream>
#include <vector>
#include <map>
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

		bool eigen_info_; // Whether eigenvectors have been computed

	public:
		// When local dimension is not given
		FloEvolVanilla(int size): EvolMatrix< ComplexEigenSolver<MatrixXcd> >(size),
			constructed_(false), eigen_info_(false){eigen.resize(1); eigen_name.resize(1,"");}

		// When local dimension is given
		FloEvolVanilla(int size, int local_dim): 
			EvolMatrix< ComplexEigenSolver<MatrixXcd> >(size, local_dim), 
			constructed_(false), eigen_info_(false){eigen.resize(1); eigen_name.resize(1,"");}

		// Diagnolize time evolution matrix with eigenvectors kept
		void Evol_Diag() {
			if (constructed_){
				eigen[0].compute(evol_op_);
				eigen_info_ = true;
				eigen_name[0] = "Full";
			}
			else{
				cout << "The matrix for diagonalization does not exist." <<endl;
				abort();
			}
		} 

		// Diagnolize time evolution matrix, user can determine whether eigenvectors are kept
		// False is not kept; True is kept
		void Evol_Diag(bool keep) {
			if (constructed_){
				eigen[0].compute(evol_op_, keep);
				eigen_info_ = keep;
				eigen_name[0] = "Full";
			}
			else{
				cout << "The matrix for diagonalization does not exist." <<endl;
				abort();
			}
		}

		// Return the string format of representation string stream 
		string Repr() const {return repr_.str();} 

		// Return the type of the model
		string Type() const {return type_;} 

		// Return the type of the basis that eigenstates are written in
		string Eigen_Type() const {return "Basic";}

		// Return dimension of each sector. Here only 1 sector exists, so total dimension
		// is returned.
		vector<int> Get_Sector_Dim() const{
			vector<int> dim(1);
			dim[0] = dim_;
			return dim;
		}

		// Erase the evolutionary operator
		void Evol_Erase() {evol_op_.resize(0,0); constructed_ = false;}

		// Construct Transition Matrix
		void Transition_Compute(TransitionMatrix&, const string&) const;

		// Return evol_op no matter what the input is
		const MatrixXcd& Get_U(string a) const {
			if (constructed_) return evol_op_;
			else{
				cout << Repr() << " has not been constructed yet." << endl;
				abort();
			}
		}

		const MatrixXcd& Get_U(int a=0) const {
			if (constructed_) return evol_op_;
			else{
				cout << Repr() << " has not been constructed yet." << endl;
				abort();
			}
		}

		// No real Hamiltonian to return
		const MatrixXd& Get_real_H(string a) const{
			cout << "No real Hamiltonian is constructed for " << Repr() << endl;
			abort();
		}

		const MatrixXd& Get_real_H(int a = 0) const{
			cout << "No real Hamiltonian is constructed for " << Repr() << endl;
			abort();
		}

		virtual ~FloEvolVanilla() {};
};

//=========================================================================================

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

		// Record the two binary states used for even parity states
		vector<vector<int> > even_parity_;

		// Record the two binary states used for odd parity states
		vector<vector<int> > odd_parity_;

		stringstream repr_; // Representation string stream of the model
		string type_; // Type string of the model

		bool eigen_info_; // Whether eigenvectors have been computed

	public:
		// When local dimension is not given
		FloEvolParity(int size): EvolMatrix< ComplexEigenSolver<MatrixXcd> >(size),
			constructed_(false), eigen_info_(false){eigen.resize(2); eigen_name.resize(2,"");}

		// When local dimension is given
		FloEvolParity(int size, int local_dim): 
			EvolMatrix< ComplexEigenSolver<MatrixXcd> >(size, local_dim), 
			constructed_(false), eigen_info_(false){eigen.resize(2);eigen_name.resize(2,"");}

		// Diagnolize time evolution matrix with eigenvectors kept
		void Evol_Diag(){
			if (constructed_){
				eigen[0].compute(evol_op_even_);
				// In case there is no odd sector for small chain
				if (evol_op_odd_.rows()>0) eigen[1].compute(evol_op_odd_);
				eigen_name[0] = "Even";
				eigen_name[1] = "Odd";
				eigen_info_ = true;
			}
			else{
				cout << "The matrix for diagonalization does not exist." <<endl;
				abort();
			}
		}; 

		// Diagnolize time evolution matrix, user can determine whether eigenvectors are kept
		// False is not kept; True is kept
		void Evol_Diag(bool keep){
			if (constructed_){
				eigen[0].compute(evol_op_even_, keep);
				// In case there is no odd sector for small chain
			if (evol_op_odd_.rows()>0) eigen[1].compute(evol_op_odd_, keep);
			eigen_name[0] = "Even";
			eigen_name[1] = "Odd";
			eigen_info_ = keep;
			}
			else{
				cout << "The matrix for diagonalization does not exist." <<endl;
				abort();
			}
		};

		// Return the string format of representation string stream 
		string Repr() const {return repr_.str();} 

		// Return the type of the model
		string Type() const {return type_;} 

		// Return the type of the basis that eigenstates are written in
		string Eigen_Type() const {return "Parity";}

		// Erase the evolutionary operator
		void Evol_Erase() {evol_op_even_.resize(0,0); evol_op_odd_.resize(0,0); 
					       constructed_ = false;}

		// Construct Transition Matrix
		void Transition_Compute(TransitionMatrix&, const string&) const;

		// Return even or odd evol_op according to the string or integer, where 0 is 
		// even and 1 is odd.
		const MatrixXcd& Get_U(string) const;

		const MatrixXcd& Get_U(int) const;

		// No real Hamiltonian to return
		const MatrixXd& Get_real_H(string a) const{
			cout << "No real Hamiltonian is constructed for " << Repr() << endl;
			abort();
		}

		const MatrixXd& Get_real_H(int a = 0) const{
			cout << "No real Hamiltonian is constructed for " << Repr() << endl;
			abort();
		}

		virtual ~FloEvolParity() {};
};

//======================================================================================

/*
 * The general evolution operator for Floquet system with multiple sectors. It can also be 
 * used for Markov Chain dynamics. There are multiple evolutionary operators. 
 */

class FloEvolMultiSec : public EvolMatrix< ComplexEigenSolver<MatrixXcd> >
{		
	protected:
		vector<MatrixXcd> evol_op_; // Set of evolution operators
		map<string,int> op_name_; // Map operator's name to the position in above vector

		bool constructed_; // Whether the matrix has been constructed and not erased

		stringstream repr_; // Representation string stream of the model
		string type_; // Type string of the model

		vector<bool> eigen_info_; // Whether eigenvectors have been computed for each sector

		virtual void Eigen_Name_Construct_() = 0; // Construct eigen_name during diagonalization

	public:
		// When local dimension is not given
		FloEvolMultiSec(int size, int type_num): EvolMatrix< ComplexEigenSolver<MatrixXcd> >(size),
			constructed_(false), eigen_info_(false){evol_op_.resize(type_num); 
				eigen.resize(type_num); eigen_info_.resize(type_num, false); 
				eigen_name.resize(type_num, "");}

		// When local dimension is given
		FloEvolMultiSec(int size, int local_dim, int type_num): 
			EvolMatrix< ComplexEigenSolver<MatrixXcd> >(size, local_dim), 
			constructed_(false), eigen_info_(false){evol_op_.resize(type_num); 
				eigen.resize(type_num);eigen_info_.resize(type_num, false);
				eigen_name.resize(type_num, "");}

		// Diagnolize time evolution matrix. Eigenvectors are kept.
		void Evol_Diag(){
			if (constructed_){
				for (int i=0; i< evol_op_.size(); i++){
					eigen[i].compute(evol_op_[i]);
					eigen_info_[i] = true;
				}
				Eigen_Name_Construct_();
			}
			else{
				cout << "The matrix for diagonalization does not exist." <<endl;
				abort();
			}
		};

		// Diagnolize time evolution matrix, user can determine whether eigenvectors are kept
		// False is not kept; True is kept.  
		void Evol_Diag(bool keep){
			if (constructed_){
				for (int i=0; i< evol_op_.size(); i++){
					eigen[i].compute(evol_op_[i], keep);
					eigen_info_[i] = keep;
				}
				Eigen_Name_Construct_();
			}
			else{
				cout << "The matrix for diagonalization does not exist." <<endl;
				abort();
			}
		};

		// Diagnolize time evolution matrix, user can determine for each sector that whether 
		// eigenvectors are kept. False is not kept; True is kept. If the bool vector has 
		// fewer elements than the number of sectors, the remaining ones eigenvectors are kept.
		void Evol_Diag(vector<bool> keep){
			if (constructed_){
				for (int i=0; i< evol_op_.size(); i++){
					if (i<keep.size()){
						eigen[i].compute(evol_op_[i], keep[i]);
						eigen_info_[i] = keep[i];
					}
					else{
						eigen[i].compute(evol_op_[i]);
						eigen_info_[i] = true;
					}
				}
				Eigen_Name_Construct_();
			}
			else{
				cout << "The matrix for diagonalization does not exist." <<endl;
				abort();
			}
		};

		// Return the string format of representation string stream 
		string Repr() const {return repr_.str();} 

		// Return the type of the model
		string Type() const {return type_;} 

		// Erase the evolutionary operator
		void Evol_Erase() {for (int i=0; i<evol_op_.size();i++) evol_op_[i].resize(0,0); 
					       constructed_ = false;}

		// Return even or odd evol_op according to the string or integer
		const MatrixXcd& Get_U(string) const;

		const MatrixXcd& Get_U(int) const;

		// No real Hamiltonian to return
		const MatrixXd& Get_real_H(string a) const{
			cout << "No real Hamiltonian is constructed for " << Repr() << endl;
			abort();
		}

		const MatrixXd& Get_real_H(int a = 0) const{
			cout << "No real Hamiltonian is constructed for " << Repr() << endl;
			abort();
		}

		virtual ~FloEvolMultiSec() {};
};

#endif
