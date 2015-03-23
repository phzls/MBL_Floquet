//
// Created by Liangsheng Zhang on 3/23/15.
//

#include <complex>
#include <iostream>
#include "constants.h"
#include "flo_evol_model.h"
#include "eigen_output.h"
#include "randomc.h"

extern CRandomMersenne RanGen_mersenne; // points in [0,1)

using namespace std;

/**
 ** Implementation of FloEvolXXZRandomSimp class
 **/

/*
 * Construct the representation string and abstract type of the class.
 */
void FloEvolXXZRandomSimp::Repr_Init_(){
    repr_ << "XXZ_Random_Floquet_L=" << size_ << ",W=" << W_;
    type_ = "XXZ_Random_Simple_Floquet";
}

/*
 * Initialize random numbers
 */
void FloEvolXXZRandomSimp::Evol_Para_Init() {
    random_h_.resize(size_);

    // Random Fields
    for (int i=0; i<size_;i++){
        double u = RanGen_mersenne.Random();
        random_h_[i] = 1 + W_ * (2*u-1);
    }

    if (debug_){
        cout << "Random longitude field:" << endl;
        for (int i=0; i<size_;i++) cout << random_h_[i] << endl;
    }
}

/*
 * Construct the time evolution matrix
 */
void FloEvolXXZRandomSimp::Evol_Construct() {

    if (!constructed_){
        evol_op_ = MatrixXcd::Zero(dim_, dim_);
        constructed_ = true;
    }
    else{
        cout << "Evolution operator has been constructed." << endl;
        abort();
    }

    MatrixXcd evol_x = MatrixXcd::Zero(dim_,dim_);
    MatrixXcd evol_z = MatrixXcd::Zero(dim_,dim_);

    Evol_X_Construct_(evol_x);
    Evol_Z_Construct_(evol_z);

    evol_op_ = evol_x * evol_z;

    if (debug_){
        cout << "X part of evolution matrix:" << endl;
        complex_matrix_write(evol_x);
        cout << endl;

        cout << "Z part of evolution matrix:" << endl;
        complex_matrix_write(evol_z);
        cout << endl;

        cout << "Evolution matrix:" << endl;
        complex_matrix_write(evol_op_);
        cout << endl;
    }
}

/*
 * Construct X part of the evolution matrix
 */
void FloEvolXXZRandomSimp::Evol_X_Construct_(MatrixXcd & evol_x) {
    for (int i=0; i<dim_; i++){
        int flip_spin = 1;
        for (int j=0;j<size_;j++){
            int pos = i ^ flip_spin;
            evol_x(pos,i) += (1 - W_);
            flip_spin = flip_spin << 1;
        }
    }

    // Check if the matrix is Hermitian
    for (int i=0; i<evol_x.cols(); i++){
        for (int j=i+1; j< evol_x.rows();j++){
            if (abs(evol_x(j,i) - evol_x(i,j)) > 1.0e-10){
                cout << "x part is not Hermitian at row " << j <<" and col " << i << endl;
                cout << "(i,j): " << evol_x(i,j) << endl;
                cout << "(j,i): " << evol_x(j,i) << endl;
                abort();
            }
        }
    }

    SelfAdjointEigenSolver<MatrixXcd> x_eigen;

    x_eigen.compute(evol_x);

    for (int i=0; i<dim_; i++){
        for (int j=0; j< dim_; j++){
            if (i==j){
                evol_x(i,j) = exp( - Complex_I * x_eigen.eigenvalues()[i] );
            }
            else{
                evol_x(i,j) = complex<double>(0,0);
            }
        }
    }

    evol_x = x_eigen.eigenvectors() * evol_x * x_eigen.eigenvectors().adjoint();
}

/*
 *  Construct z part of the evolution matrix
 */
void FloEvolXXZRandomSimp::Evol_Z_Construct_(MatrixXcd & evol_z) {
    for (int i=0; i<dim_;i++){

        int state = i;
        int spin = 0, prev_spin = 0;
        double value = 0; // value which will be exponentiated

        for (int j=0; j<size_;j++){
            if (j>0) prev_spin = spin;
            spin = 2*(state & 1) - 1;
            state = state >> 1;

            value += random_h_[j] * spin;

            if (j>0) value += spin * prev_spin;
        }

        evol_z(i,i) = exp( - Complex_I * complex<double>(value,0) );
    }
}