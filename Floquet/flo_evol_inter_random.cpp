#include <complex>
#include <iostream>
#include "constants.h"
#include "flo_evol_model.h"
#include "mtrand.h"
#include "eigen_output.h"

extern MTRand u1rand; // points in [0,1)
extern MTRand_closed u2rand; // points in [0,1]

using namespace std;

/**
 ** Implementation of FloEvolInterRandom class.
 **/

/*
 * Construct the representation string and abstract type of the class.
 */
void FloEvolInterRandom::Repr_Init_(){
	repr_ << "Inter_Random_Floquet_L=" << size_ << ",J=" << param_.J 
		  << ",g=" << param_.g << ",h=" << param_.h <<",tau="<< param_.tau;
	type_ = "Inter_Random_Floquet";
}


/*
 * Construct time evolution operator. 
 */
void FloEvolInterRandom::Evol_Construct(){
	// If the matrix has not been constructed
	if (!constructed_){
		evol_op_ = MatrixXcd::Zero(dim_, dim_);
		constructed_ = true;
	}
	else{
		cout << "Evolution operator has been constructed." << endl;
		abort();
	}

	/* For random unitary part Ur, we generate angle uniformly random on [-J*Pi,J*pi) and a random
	 * unit vector on 3D sphere, where z uniformly on [0,1] and phi uniformly on [0,2*pi). For
	 * xxz part Uxxz, the time step tau is taken to be J*tau.
	 */

	FloEvolXXZ* floquet_xxz = new FloEvolXXZ(param_.size, param_.J * param_.tau, param_.g, 
		param_.h, debug_);

	floquet_xxz -> Evol_Construct();

	// even XXZ part of time evolution operator
	MatrixXcd Uxxz_even = floquet_xxz -> Evol_Op_Even();

	// odd XXZ part of time evolution operator
	MatrixXcd Uxxz_odd = floquet_xxz -> Evol_Op_Odd();
	
	// XXZ part of time evolution operator
	MatrixXcd Uxxz = MatrixXcd::Zero(dim_, dim_);

	for (int i=0; i< Uxxz_even.cols(); i++){
		for (int j=0; j<Uxxz_even.rows(); j++){
			Uxxz(j,i) = Uxxz_even(j,i);
		}
	}

	int even_rank = Uxxz_even.cols();
	if (even_rank != Uxxz_even.rows()){
		cout << "Uxxz_even is not square." << endl;
		abort();
	}

	for (int i=0; i<Uxxz_odd.cols(); i++){
		for (int j=0; j<Uxxz_odd.rows();j++){
			Uxxz(j+even_rank, i+even_rank) = Uxxz_odd(j,i);
		}
	}

	delete floquet_xxz;

	FloEvolRandomRotation* floquet_r = new FloEvolRandomRotation(param_.size, param_.tau, 
		param_.J, -param_.J*Pi, param_.J*Pi, debug_);

	floquet_r -> Evol_Para_Init();
	floquet_r -> Evol_Construct();

	// Random part of time evolution operator
	MatrixXcd Ur = floquet_r -> Evol_Op();

	evol_op_ = Ur * Uxxz;

	// In the debug mode, print out all matrices
	if (debug_){
		cout<<"Ur:"<<endl;
		complex_matrix_write(Ur);
		cout << endl;

		cout<<"Uxxz:"<<endl;
		complex_matrix_write(Uxxz);
		cout << endl;

		cout<<"Final matrix:"<<endl;
		complex_matrix_write(evol_op_);
		cout<<endl;
	}
}

/*
 * Return evol_op; not meant to be called in polymorphism
 */
const MatrixXcd& FloEvolRandom::Evol_Op() const{
	if (constructed_) return evol_op_;
	else{
		cout << Repr() << " has not been constructed yet." << endl;
		abort();
	}
}