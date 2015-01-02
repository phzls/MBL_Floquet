#include <complex>
#include <iostream>
#include "constants.h"
#include "flo_evol_model.h"
#include "mtrand.h"

extern MTRand u1rand; // points in [0,1)
extern MTRand_closed u2rand; // points in [0,1]

using namespace std;

/**
 ** Implementation of FloEvolRandom class.
 **/

/*
 * Construct the representation string and abstract type of the class.
 */
void FloEvolRandom::Repr_Init_(){
	repr_ << "Random_Floquet_L=" << size_ << ",J=" << param_.J <<",tau="
		  << param_.tau;
	type_ = "Random_Floquet";
}

/*
 * Initialize parameters of the model.
 */
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

	if (debug_){
		cout<<"All random numbers:"<<endl;
		for (int i=0; i<size_; i++){
			cout<<i<<"  ";
			for (int j=0; j<su2_angle_[i].size(); j++) cout<<su2_angle_[i][j]<<"  ";
			cout<<endl;
		}
	}
}

/*
 * Construct time evolution operator. 
 */
void FloEvolRandom::Evol_Construct(){
	// If the matrix has not been constructed
	if (!constructed_){
		evol_op_ = MatrixXcd::Zero(dim_, dim_);
		constructed_ = true;
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

/*
 * Construct random Ux part.
 */
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

/*
 * Construct Uz part.
 */
void FloEvolRandom::Evol_Z_Construct_(MatrixXcd& Uz){
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

//===============================================================================================

/**
 ** Implementation of FloEvolRandomRotation class.
 **/

/*
 * Construct the representation string and abstract type of the class.
 */
void FloEvolRandomRotation::Repr_Init_(){
	repr_ << "Random_Rotation_Floquet_L=" << size_ << ",J=" << param_.J <<",tau="
		  << param_.tau;
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

	// We generate angle uniformly random on [0,2*pi) and a random unit vector on 3D
	// sphere, where z uniformly on [0,1] and phi uniformly on [0,2*pi). It rotates a
	// density matrix on block sphere w.r.t the axis of unit vector by amount of angle
	const double angle_range = angle_sup_ - angle_min_;
	for (int i=0; i<size_; i++){
		angle_[i] = angle_min_ + angle_range* u1rand();

		double z = 2 * u2rand() - 1;
		double phi = 2 * Pi * u1rand();

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
	if (angle_sup_ < angle_min_){
		cout << "Rotation angle_sup should be no less than angle_min." << endl;
		abort();
	}
	if (angle_min_ < 0 || angle_sup_ > 2*Pi){
		cout << "Rotation angle should be between 0 and 2*pi." << endl;
		abort();
	}
}