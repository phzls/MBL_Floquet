#include <complex>
#include <iostream>
#include "constants.h"
#include "flo_evol_model.h"
#include "eigen_output.h"

using namespace std;

/**
 ** Implementation of FloEvolInterRandom class.
 **/

/*
 * Construct the representation string and abstract type of the class.
 */
void FloEvolMarkovInterRandom::Repr_Init_(){
	repr_ << "Markov_Inter_Random_Floquet_L=" << size_ << ",J=" << param_.J 
		  << ",g=" << param_.g << ",h=" << param_.h <<",tau="<< param_.tau;
	type_ = "Markov_Inter_Random_Floquet";
}

/*
 * Initialize the map linking name and position of evol_op
 */
void FloEvolMarkovInterRandom::Op_Name_Init_(){
	op_name_["down"] = 0;
	op_name_["Down"] = 0;
	op_name_["up"] = 1;
	op_name_["Up"] = 1;
}


/*
 * Construct time evolution operator. 
 */
void FloEvolMarkovInterRandom::Evol_Construct(){
	// If the matrix has not been constructed
	if (!constructed_){
		constructed_ = true;
	}
	else{
		cout << "Evolution operator has been constructed." << endl;
		abort();
	}

	/* For random unitary part Ur, we generate angle uniformly random on [-(1-J)*Pi,(1-J)*pi)  
	 * and a random unit vector on 3D sphere, where z uniformly on [0,1] and phi uniformly on 
	 * [0,2*pi). Also J=1 so that no z part. For xxz part Uxxz, the time step tau is taken to be
	 * J*tau. The xxz part natively is written in parity states, so it needs to be converted into
	 * binary basis.
	 */

	FloEvolRandomRotation* floquet_r = new FloEvolRandomRotation(param_.size, param_.tau, 
		0, -(1-param_.J)*Pi, (1-param_.J)*Pi, debug_);

	floquet_r -> Evol_Para_Init();
	floquet_r -> Evol_Construct();

	// Random part of time evolution operator
	MatrixXcd Ur = floquet_r -> Get_U();

	MatrixXcd Uxxz = MatrixXcd::Zero(dim_, dim_);

	Bath_XXZ_Construct_(Uxxz, "down");

	evol_op_[0] = Ur * Uxxz;

	// In the debug mode, print out all matrices
	if (debug_){
		cout << "For bath spin down:" << endl;
		cout<<"Ur:"<<endl;
		complex_matrix_write(Ur);
		cout << endl;

		cout<<"Uxxz:"<<endl;
		complex_matrix_write(Uxxz);
		cout << endl;

		cout<<"Final matrix:"<<endl;
		complex_matrix_write(evol_op_[0]);
		cout<<endl;
	}

	Bath_XXZ_Construct_(Uxxz, "up");

	evol_op_[1] = Ur * Uxxz;

	// In the debug mode, print out all matrices
	if (debug_){
		cout << "For bath spin up:" << endl;
		cout<<"Ur:"<<endl;
		complex_matrix_write(Ur);
		cout << endl;

		cout<<"Uxxz:"<<endl;
		complex_matrix_write(Uxxz);
		cout << endl;

		cout<<"Final matrix:"<<endl;
		complex_matrix_write(evol_op_[1]);
		cout<<endl;
	}

	delete floquet_r;
}

/*
 * Construct the xxz part of the operator
 */
void FloEvolMarkovInterRandom::Bath_XXZ_Construct_(MatrixXcd& U, string type){
	// From the bath spin direction, determines the sign
	int sign;
	if (type == "down" || type == "Down") sign = -1;
	else if (type == "up" || type == "Up") sign = 1;
	else{
		cout << "Type unknown." << endl;
		abort();
	}

	// Effective time step
	const double tau = param_.tau;

	// Construct z part Hamiltonian first, which is diagonal
	U = MatrixXcd::Zero(dim_,dim_);

	SelfAdjointEigenSolver<MatrixXcd> U_eigen;

	for (int i=0; i<dim_; i++){
		// Diagonal terms
		int state = i;
		int prev_spin = 2* (state & 1) - 1 ;
		state = state >> 1;
		
		// Coupling to the bath is not reduced by the coupling J
		U(i,i) += (param_.J*param_.h + sign) * prev_spin; // Contain part from the bath

		for (int j=1;j<size_;j++){
			int current_spin = 2* (state & 1) - 1;
			state = state >> 1;

			U(i,i) += param_.J* (param_.h * current_spin + current_spin * prev_spin);

			prev_spin = current_spin;
		}
	}

	// Write out Hamiltonian matrix
	if (debug_){
		cout << "XXZ Hamiltonian z part with bath " << type << " :" << endl;
		complex_matrix_write(U);
		cout << endl;
	}

	for (int i=0; i<dim_; i++){
		for (int j=0; j< dim_; j++){
			if (i==j){
				U(i,j) = exp(-Complex_I * complex<double>(tau,0) * U(i,j));
			}
			else{
				U(i,j) = complex<double>(0,0);
			}
		}
	}

	if (debug_){
		cout << "XXZ Floquet z part with bath " << type << " :" << endl;
		complex_matrix_write(U);
		cout << endl;
	}

	// x part, which is off-diagonal part
	MatrixXcd Ux = MatrixXcd::Zero(dim_,dim_);

	for (int i=0; i<dim_; i++){
		// Off-diagonal terms
		int site = 1;
		for (int j=0; j<size_; j++){
			int target = i ^ site;
			site = site << 1;
			Ux(target,i) += param_.J * param_.g; 
		}
	}

	// Write out Hamiltonian matrix
	if (debug_){
		cout << "XXZ Hamiltonian x part with bath " << type << " :" << endl;
		complex_matrix_write(Ux);
		cout << endl;
	}

	U_eigen.compute(Ux);

	for (int i=0; i<dim_; i++){
		for (int j=0; j< dim_; j++){
			if (i==j){
				Ux(i,j) = exp(-Complex_I * complex<double>(tau,0) * U_eigen.eigenvalues()(i));
			}
			else{
				Ux(i,j) = complex<double>(0,0);
			}
		}
	}

	Ux = U_eigen.eigenvectors() * Ux * U_eigen.eigenvectors().adjoint();

	if (debug_){
		cout << "XXZ Floquet z part with bath " << type << " :" << endl;
		complex_matrix_write(U);
		cout << endl;
	}

	U = Ux * U;
}

void FloEvolMarkovInterRandom::Eigen_Name_Construct_(){
	eigen_name[op_name_["Down"]] = "Bath Down";
	eigen_name[op_name_["Up"]] = "Bath Up";
}