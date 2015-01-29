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

	/* For random unitary part Ur, we generate angle uniformly random on [-(1-J)*Pi,(1-J)*pi) and 
	 * a random unit vector on 3D sphere, where z uniformly on [0,1] and phi uniformly on 
	 * [0,2*pi). Also J=1 so that no z part. For xxz part Uxxz, the time step tau is taken to be
	 * J*tau. The xxz part natively is written in parity states, so it needs to be converted into
	 * binary basis.
	 */

	FloEvolXXZ* floquet_xxz = new FloEvolXXZ(param_.size, param_.J * param_.tau, param_.g, 
		param_.h, debug_);

	floquet_xxz -> Evol_Construct();
	floquet_xxz -> Evol_Diag();

	TransitionMatrix transition;
	floquet_xxz -> Transition_Compute(transition, "Basic_Parity");

	// even XXZ part of time evolution operator
	MatrixXcd Uxxz_even = transition.Matrix("Basic_Even") * ( floquet_xxz -> Get_U("Even") )
		* transition.Matrix("Basic_Even").adjoint();

	// odd XXZ part of time evolution operator
	MatrixXcd Uxxz_odd = transition.Matrix("Basic_Odd") * ( floquet_xxz -> Get_U("Odd") )
	 * transition.Matrix("Basic_Odd").adjoint();

	transition.Erase_All();
	
	// XXZ part of time evolution operator
	MatrixXcd Uxxz = Uxxz_even + Uxxz_odd;

	delete floquet_xxz;

	FloEvolRandomRotation* floquet_r = new FloEvolRandomRotation(param_.size, param_.tau, 
		0, -(1-param_.J)*Pi, (1-param_.J)*Pi, debug_);

	floquet_r -> Evol_Para_Init();
	floquet_r -> Evol_Construct();

	// Random part of time evolution operator
	MatrixXcd Ur = floquet_r -> Get_U();

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

	delete floquet_r;
}