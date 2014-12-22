#include <complex>
#include <iostream>
#include "parameter.h"
#include "flo_evol_model.h"
#include "mtrand.h"

extern MTRand u1rand; // points in [0,1)
extern MTRand_closed u2rand; // points in [0,1]

using namespace std;

void FloEvolRandom::Repr_Init_(){
	repr_ << "Random_Floquet_L=" << size_ << ",J=" << param_.J <<",tau="
		  << param_.tau;
}

void FloEvolRandom::Evol_Para_Init(){
	su2_angle_.resize(size_);

	for (int i=0; i<size_; i++){
		su2_angle_[i].resize(3);
	}

	// We generate psi and xi randomly in [0, 2*pi), and (sin(phi))^2 randomly from 
	// [0, 1] where phi belongs to [0, pi/2]
	for (int i=0; i<size_; i++){
		su2_angle_[i][0] = 2*Pi*u1rand(); // Angle psi
		su2_angle_[i][1] = 2*Pi*u1rand(); // Angle xi
		su2_angle_[i][2] = u2rand(); // (sin(phi))^2
	}

	cout<<"All random numbers:"<<endl;
	for (int i=0; i<size_; i++){
		cout<<i<<"  ";
		for (int j=0; j<su2_angle_[i].size(); j++) cout<<su2_angle_[i][j]<<"  ";
		cout<<endl;
	}
}

void FloEvolRandom::Evol_Construct(){
	// Ux part of time evolution operator
	MatrixXcd Ux = MatrixXcd::Zero(dim_, dim_);
	Evol_Site_Construct_(Ux);

	// Uz part of time evolution operator
	MatrixXcd Uz = MatrixXcd::Zero(dim_, dim_);
	Evol_Z_Construct_(Uz);

	evol_op_ = Ux * Uz;

	// In the debug mode, print out all matrices
	if (debug_){
		cout<<"All random numbers:"<<endl;
		for (int i=0; i<size_; i++){
			cout<<i<<"  ";
			for (int j=0; j<su2_angle_[i].size(); j++) cout<<su2_angle_[i][j]<<"  ";
			cout<<endl;
		}
		cout<<endl;

		cout<<"Ux:"<<endl;
		for (int i=0; i< Ux.rows();i++){
			for (int j=0; j< Ux.cols();j++){
				cout<<"("<<real(Ux(i,j))<<","<<imag(Ux(i,j))<<")  ";
			}
			cout<<endl;
		}
		cout<<endl;

		cout<<"Uz:"<<endl;
		for (int i=0; i< Uz.rows();i++){
			for (int j=0; j< Uz.cols();j++){
				cout<<"("<<real(Uz(i,j))<<","<<imag(Uz(i,j))<<")  ";
			}
			cout<<endl;
		}
		cout<<endl;

		cout<<"Final matrix:"<<endl;
		for (int i=0; i< evol_op_.rows();i++){
			for (int j=0; j< evol_op_.cols();j++){
				cout<<"("<<real(evol_op_(i,j))<<","<<imag(evol_op_(i,j))<<")  ";
			}
			cout<<endl;
		}
		cout<<endl;
	}
}

void FloEvolRandom::Evol_Site_Construct_(MatrixXcd& Ux){
	int current_dim = 1; // The current dimension of the matrix in the Kroneker product

	// Random SU(2) matrix at each site
	MatrixXcd U(2,2);

	Ux(0,0) = Complex_one;

	for (int site = 0; site < size_; site++){

		// Construct U matrix at each site
		U(0,0) = exp(Complex_I * complex<double>(su2_angle_[site][0],0)) *
				 complex<double>(sqrt(1-su2_angle_[site][2]),0);
		U(0,1) = exp(Complex_I * complex<double>(su2_angle_[site][1],0)) * 
				 complex<double>(sqrt(su2_angle_[site][2]),0);
		U(1,0) = -conj(U(0,1));
		U(1,1) = exp(-Complex_I * complex<double>(su2_angle_[site][0],0)) * 
				 complex<double>(sqrt(1-su2_angle_[site][2]),0);

		if (debug_){
			cout<<"At site: "<<site<<" U:"<<endl;
			for (int i=0; i< U.rows();i++){
				for (int j=0; j< U.cols();j++){
					cout<<"("<<real(U(i,j))<<","<<imag(U(i,j))<<")  ";
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

void FloEvolRandom::Evol_Z_Construct_(MatrixXcd& Uz){
	/*
	We are working in S_z basis, and 0 means down while 1 means up when each state is 
	encoded in a number. The total state will be encoded by a single number between 0
	and dim_, whose binary representation, from left to right, gives each spin orientation
	at each site, from left to right. When computing the matrix element, -1 is for down
	spin and 1 is for up spin.
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