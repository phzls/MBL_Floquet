#include <complex>
#include <iostream>
#include "constants.h"
#include "flo_evol_model.h"
#include "randomc.h"

extern CRandomMersenne RanGen_mersenne; // points in [0,1)

using namespace std;

/**
 ** Implementation of FloEvolRandomRotation class.
 **/

/*
 * Construct the representation string and abstract type of the class.
 */
void FloEvolRandomRotation::Repr_Init_(){
	repr_ << "Random_Rotation_Floquet_L=" << size_ << ",J=" << param_.J 
		  << ",angle_min=" << param_.angle_min << ",angle_sup=" << param_.angle_sup
		  <<",tau="<< param_.tau;
	type_ = "Random_Rotation_Floquet";
}

/*
 * Initialize parameters of the model.
 */
void FloEvolRandomRotation::Evol_Para_Init(){
	angle_.resize(size_);
	axis_.resize(size_);

	for (int i=0; i<size_; i++){
		axis_[i].resize(3);
	}

	// We generate angle uniformly random on [angle_min,angle_sup) and a random unit vector on 3D
	// sphere, where z uniformly on [0,1] and phi uniformly on [0,2*pi). It rotates a
	// density matrix on block sphere w.r.t the axis of unit vector by amount of angle
	const double angle_range = param_.angle_sup - param_.angle_min;
	for (int i=0; i<size_; i++){
		angle_[i] = param_.angle_min + angle_range* RanGen_mersenne.Random();

		double z = 2 * RanGen_mersenne.Random() - 1;
		double phi = 2 * Pi * RanGen_mersenne.Random();

		if (debug_) cout << "angle: " << angle_[i] << " phi: " << phi << endl;

		axis_[i][0] = sqrt(1 - z*z) * cos(phi);
		axis_[i][1] = sqrt(1 - z*z) * sin(phi);
		axis_[i][2] = z;
	}

	if (debug_){
		cout<<"All random numbers:"<<endl;
		for (int i=0; i<size_; i++){
			cout<<i<<"  ";
			cout << angle_[i] << " ";
			for (int j=0; j<axis_[i].size(); j++) cout<<axis_[i][j]<<"  ";
			cout<<endl;
		}
	}
}

/*
 * Construct time evolution operator. 
 */
void FloEvolRandomRotation::Evol_Construct(){
	// If the matrix has not been constructed
	if (!constructed_){
		evol_op_ = MatrixXcd::Zero(dim_, dim_);
		constructed_ = true;
	}
	else{
		cout << "Evolution operator has been constructed." << endl;
		abort();
	}

	// Ux part of time evolution operator
	MatrixXcd Ux = MatrixXcd::Zero(dim_, dim_);
	Evol_Site_Construct_(Ux);

	// Uz part of time evolution operator
	MatrixXcd Uz = MatrixXcd::Zero(dim_, dim_);
	Evol_Z_Construct_(Uz);

	evol_op_ = Ux * Uz;

	// In the debug mode, print out all matrices
	if (debug_){
		cout<<"Ux:"<<endl;
		for (int i=0; i< Ux.rows();i++){
			for (int j=0; j< Ux.cols();j++){
				cout<<real(Ux(i,j));
				if (imag(Ux(i,j))<0) cout << imag(Ux(i,j)) << "j  "; 
				else cout << "+" << imag(Ux(i,j)) << "j  ";
			}
			cout<<endl;
		}
		cout<<endl;

		cout<<"Uz:"<<endl;
		for (int i=0; i< Uz.rows();i++){
			for (int j=0; j< Uz.cols();j++){
				cout<<real(Uz(i,j));
				if (imag(Uz(i,j))<0) cout << imag(Uz(i,j)) << "j  "; 
				else cout << "+" << imag(Uz(i,j)) << "j  ";
			}
			cout<<endl;
		}
		cout<<endl;

		cout<<"Final matrix:"<<endl;
		for (int i=0; i< evol_op_.rows();i++){
			for (int j=0; j< evol_op_.cols();j++){
				cout<<real(evol_op_(i,j));
				if (imag(evol_op_(i,j))<0) cout << imag(evol_op_(i,j)) << "j  "; 
				else cout << "+" << imag(evol_op_(i,j)) << "j  ";
			}
			cout<<endl;
		}
		cout<<endl;
	}
}

/*
 * Construct random Ux part.
 */
void FloEvolRandomRotation::Evol_Site_Construct_(MatrixXcd& Ux){
	int current_dim = 1; // The current dimension of the matrix in the Kroneker product

	// Random SU(2) matrix at each site
	MatrixXcd U(2,2);

	Ux(0,0) = Complex_one;

	for (int site = 0; site < size_; site++){

		// Construct U matrix at each site
		U(0,0) = complex<double>(cos(angle_[site]/2), -axis_[site][2]*sin(angle_[site]/2));
		U(0,1) = complex<double>(-sin(angle_[site]/2)*axis_[site][1], 
								 -sin(angle_[site]/2)*axis_[site][0]);
		U(1,0) = -conj(U(0,1));
		U(1,1) = conj(U(0,0));	

		if (debug_){
			cout<<"At site: "<<site<<" U:"<<endl;
			for (int i=0; i< U.rows();i++){
				for (int j=0; j< U.cols();j++){
					cout<<real(U(i,j));
					if (imag(U(i,j))<0) cout << imag(U(i,j)) << "j  "; 
					else cout << "+" << imag(U(i,j)) << "j  ";
				}
				cout<<endl;
			}
			cout<<endl;
		}

		// Doing Kroneker Product operation with U at left
		for (int i=0; i<current_dim;i++){ // The second block in the row
			for (int j=0; j<current_dim; j++){
				Ux(i, j+current_dim) = Ux(i,j) * U(0,1);
			}
		}

		for (int i=0; i<current_dim; i++){ // The second block in the column
			for (int j=0; j<current_dim; j++){
				Ux(i+current_dim, j) = Ux(i,j) * U(1,0); 
			}
		}

		for (int i=0; i<current_dim; i++){ // The second block on the diagonal
			for (int j=0; j<current_dim; j++){
				Ux(i+current_dim, j+current_dim) = Ux(i,j) * U(1,1);
			}
		}

		for (int i=0; i<current_dim; i++){ // The first block
			for (int j=0; j<current_dim; j++){
				Ux(i,j) *= U(0,0);
			}
		}

		current_dim *= 2;
	}

	if (current_dim < dim_){
		cout << "Kroneker Product fails to produce the correct dimension."<<endl;
		abort();
	}
}

/*
 * Construct Uz part.
 */
void FloEvolRandomRotation::Evol_Z_Construct_(MatrixXcd& Uz){
	/*
	 * We are working in S_z basis, and 0 means down while 1 means up when each state is 
	 * encoded in a number. The total state will be encoded by a single number between 0
	 * and dim_, whose binary representation, from left to right, gives each spin
	 * orientation at each site, from left to right. When computing the matrix element, 
	 * -1 is for down spin and 1 is for up spin.
	 */

	for (int i=0; i<dim_; i++){
		int spin = i;
		int current_spin = 2 * (spin % 2) - 1;
		spin /=  2;
		int next_spin;
		double H = 0; // The logarithm of the matrix element

		for (int j=0; j<size_ - 1; j++){
			next_spin = 2 * (spin % 2) - 1;
			spin /= 2;

			H += param_.J * next_spin * current_spin;

			current_spin = next_spin;
		}
		Uz(i,i) = exp(-Complex_I * complex<double>(H,0) );
	}
}

/*
 * Check input parameters
 */
void FloEvolRandomRotation::Para_Check_(){
	if (param_.angle_sup < param_.angle_min){
		cout << "Rotation angle_sup should be no less than angle_min." << endl;
		abort();
	}
	if (param_.angle_sup - param_.angle_min > 2*Pi){
		cout << "Rotation angle should be within 2*pi." << endl;
		abort();
	}
}