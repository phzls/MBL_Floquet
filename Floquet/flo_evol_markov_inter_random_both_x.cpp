#include <complex>
#include <iostream>
#include "constants.h"
#include "flo_evol_model.h"
#include "eigen_output.h"

using namespace std;

/**
 ** Implementation of FloEvolInterRandomBoth class.
 **/

/*
 * Construct the representation string and abstract type of the class.
 */
void FloEvolMarkovInterRandomBothX::Repr_Init_(){
	repr_ << "Markov_Inter_Random_Both_X_Floquet_L=" << size_ << ",J=" << param_.J 
		  << ",g=" << param_.g << ",h=" << param_.h <<",tau="<< param_.tau;
	type_ = "Markov_Inter_Random_Both_X_Floquet";
}

/*
 * Initialize the map linking name and position of evol_op
 */
void FloEvolMarkovInterRandomBothX::Op_Name_Init_(){
	op_name_["down"] = 0;
	op_name_["Down"] = 0;
	op_name_["up"] = 1;
	op_name_["Up"] = 1;
	op_name_["isolated"] = 2;
	op_name_["Isolated"] = 2;
}


/*
 * Construct time evolution operator. 
 */
void FloEvolMarkovInterRandomBothX::Evol_Construct(){
	// If the matrix has not been constructed
	if (!constructed_){
		constructed_ = true;
	}
	else{
		cout << "Evolution operator has been constructed." << endl;
		abort();
	}

	/* We generate an inter random floquet operator, which is the isolated one, and then the 
	 * one coupled to the bath is obtained by multiplying an extra exponential at the back
	 */

	FloEvolInterRandom* floquet_iso = new FloEvolInterRandom(param_.size, param_.tau, 
		param_.J, param_.g, param_.h, debug_);

	floquet_iso -> Evol_Para_Init();
	floquet_iso -> Evol_Construct();

	MatrixXcd U = floquet_iso -> Get_U();
	evol_op_[2] = MatrixXcd::Zero(U.rows(), U.cols());

	for (int i=0; i<U.cols();i++){
		for (int j=0; j<U.rows();j++)
			evol_op_[2](j,i) = U(j,i);
	}

	delete floquet_iso;

	// In the debug mode, print out all matrices
	if (debug_){
		cout << "For isolated system:" << endl;
		cout<<"Final matrix:"<<endl;
		complex_matrix_write(evol_op_[2]);
		cout<<endl;
	}

	MatrixXcd Bath = MatrixXcd::Zero(dim_, dim_);

	Bath_Construct_(Bath, "down");

	evol_op_[0] = evol_op_[2] * Bath;

	// In the debug mode, print out all matrices
	if (debug_){
		cout << "For bath spin down:" << endl;

		cout<<"Bath:"<<endl;
		complex_matrix_write(Bath);
		cout << endl;

		cout<<"Final matrix:"<<endl;
		complex_matrix_write(evol_op_[0]);
		cout<<endl;
	}

	Bath_Construct_(Bath, "up");

	evol_op_[1] = evol_op_[2] * Bath;

	// In the debug mode, print out all matrices
	if (debug_){
		cout << "For bath spin up:" << endl;

		cout<<"Bath:"<<endl;
		complex_matrix_write(Bath);
		cout << endl;

		cout<<"Final matrix:"<<endl;
		complex_matrix_write(evol_op_[1]);
		cout<<endl;
	}
}

/*
 * Construct the bath coupling part of the operator, which is exp(-i*tau* \pm 1)
 */
void FloEvolMarkovInterRandomBothX::Bath_Construct_(MatrixXcd& U, string type){
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

	// Coupling strength to the bath
	const double K = param_.K;

	// Construct sigma_x at the right end of the chain
	U = MatrixXcd::Zero(dim_,dim_);

	MatrixXcd sigma_x = MatrixXcd::Zero(dim_,dim_);

	SelfAdjointEigenSolver<MatrixXcd> sigma_eigen;

	for (int i=0; i<dim_; i++){
		int j = i ^ 1; // Flip the rightmost spin
		sigma_x(j,i) = 1;
	}

	// Write out sigma_x matrix
	if (debug_){
		cout << "Rigtmost sigma x matrix:" << endl;
		complex_matrix_write(sigma_x);
		cout << endl;
	}

	sigma_eigen.compute(sigma_x);

	for (int i=0; i<dim_; i++){
		for (int j=0; j< dim_; j++){
			if (i==j){
				sigma_x(i,j) = exp(-Complex_I * complex<double>(tau*K*sign,0) * 
					sigma_eigen.eigenvalues()(i));
			}
			else{
				sigma_x(i,j) = complex<double>(0,0);
			}
		}
	}

	U = sigma_eigen.eigenvectors() * sigma_x * sigma_eigen.eigenvectors().adjoint();

	if (debug_){
		cout << "Coupling to the bath matrix:" << endl;
		complex_matrix_write(U);
		cout << endl;
	}
}

void FloEvolMarkovInterRandomBothX::Eigen_Name_Construct_(){
	eigen_name[op_name_["Up"]] = "Bath Up";
	eigen_name[op_name_["Down"]] = "Bath Down";
	eigen_name[op_name_["Isolated"]] = "Isolated";

	if (debug_){
		cout << "Eigen names:"<<endl;
		for (int i=0; i<3;i++){
			cout << i << "  " << eigen_name[i] << endl;
		}
	}
}