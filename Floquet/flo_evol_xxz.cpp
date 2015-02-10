#include <complex>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "constants.h"
#include "flo_evol_model.h"
#include "eigen_output.h"

using namespace std;

/**
 ** Implementation of xxz Floquet class
 **/

 /*
 * Construct the representation string and abstract type of the class.
 */
void FloEvolXXZ::Repr_Init_(){
	repr_ << "XXZ_Floquet_L=" << size_ << ",g=" << param_.g << ",h=" << param_.h 
		  <<",tau="<< param_.tau;
	type_ = "XXZ_Floquet";
}

void FloEvolXXZ::Evol_Construct(){
	// If the matrix has not been constructed
	if (!constructed_){
		evol_op_even_ = MatrixXcd::Zero(even_dim_, even_dim_);
		evol_op_odd_ = MatrixXcd::Zero(odd_dim_, odd_dim_);

		even_parity_.resize(even_dim_);
		for (int i=0; i<even_dim_; i++) even_parity_[i].resize(2);

		odd_parity_.resize(odd_dim_);
		for (int i=0; i<odd_dim_; i++) odd_parity_[i].resize(2);
		
		constructed_ = true;
	}
	else{
		cout << "Evolution operator has been constructed." << endl;
		abort();
	}

	Evol_General_Construct_(evol_op_even_, even_dim_);

	if (debug_){
		cout << "Even part of evolution matrix:" << endl;
		complex_matrix_write(evol_op_even_);
		cout << endl;
	}

	if (odd_dim_ > 0){
		Evol_General_Construct_(evol_op_odd_, odd_dim_);

		if (debug_){
			cout << "Odd part of evolution matrix:" << endl;
			complex_matrix_write(evol_op_odd_);
			cout << endl;
		}
	}
	
}

/*
 * A general construction function which can be used for both even and odd part
 */
void FloEvolXXZ::Evol_General_Construct_(MatrixXcd& evol_op, int dim){
	
	if (evol_op.rows() != dim){
		cout << "evol_op size is not consistent with dim in general construction." << endl;
		abort();
	}

	MatrixXcd evol_z = MatrixXcd::Zero(dim, dim);
	MatrixXcd evol_x = MatrixXcd::Zero(dim, dim);
	SelfAdjointEigenSolver<MatrixXcd> x_eigen;

	// Having some trouble passing function pointers, so use dim to check
	if (dim == even_dim_) Evol_Even_Construct_(evol_x, evol_z);
	else Evol_Odd_Construct_(evol_x, evol_z);

	x_eigen.compute(evol_x);

	#pragma omp parallel num_threads(param_.threads_num)
	{	
		#pragma omp for 
		for (int i=0; i<dim; i++){
			for (int j=0; j< dim; j++){
				if (i==j){
					evol_x(i,j) = exp(-Complex_I * complex<double>(param_.tau,0)
						* x_eigen.eigenvalues()(i));
					evol_z(i,j) = exp(-Complex_I * complex<double>(param_.tau,0)
						* evol_z(i,j));
				}
				else{
					evol_x(i,j) = complex<double>(0,0);
					evol_z(i,j) = complex<double>(0,0);
				}
			}
		}
	}

	if (debug_){
        cout << "tau is :" << param_.tau << endl;
        cout << endl;
        
		cout << "Eigenvalues of x part:" << endl;
		cout << x_eigen.eigenvalues();
		cout <<endl;

		cout << "Eigenvectors of x part:" << endl;
		complex_matrix_write(x_eigen.eigenvectors());
		cout << endl;

		complex_matrix_write(x_eigen.eigenvectors()*x_eigen.eigenvectors().adjoint());
		cout << endl;

		cout << "Diagonal x part:" << endl;
		complex_matrix_write(evol_x);
		cout << endl;

		cout << "x part:" << endl;
		complex_matrix_write(x_eigen.eigenvectors() * evol_x * x_eigen.eigenvectors().adjoint());
		cout << endl;

		cout << "z part:" << endl;
		complex_matrix_write(evol_z);
		cout << endl;
	}

	evol_op = x_eigen.eigenvectors() * evol_x * x_eigen.eigenvectors().adjoint() * evol_z;
}

/*
 * Construct some even parts of evolution operator
 */
void FloEvolXXZ::Evol_Even_Construct_(MatrixXcd& evol_x_even, MatrixXcd& evol_z_even){
	int even_counter = 0; //cannot grow more than even_dim_
    int odd_counter =0; 

    // Compute pairs and diagonal elements
    for(int i = 0; i<dim_; i++)
    {
        int pair = 0; // The pair state, which is the reverse of the current state
        int bond = 0; // Contribution from the nn interaction
        int checker = 0; // Double-counting checker

        for(int k = 0; k<even_counter+1; k++)
        {
            if(i>0 && i == even_parity_[k][1]) // For the first state, nothing to check
            {
                checker = 1;
                break;
            }
        }

        if(checker == 0)
        {
       
            if (even_counter >= even_dim_){
            	cout << "Too many even parity states." << endl;
            	abort();
            }

            int site = i; // The last binary digit indicates the rightmost spin
            int current_spin = (site & 1); // The rightmost spin
            int left_spin;

            for(int j=0; j<size_-1; j++)
            {

	            pair += current_spin << (size_-j-1); // Construct the pair state by reversing
                site = site >> 1; // Proceed to the left site
                left_spin = (site & 1);

                if(current_spin == left_spin)
                {
                    bond += 1;
                }
                else
                {
                    bond -= 1;
                }

                current_spin = left_spin;
            }

            pair += current_spin; // The leftmost site

            int spin = __builtin_popcount(i); // Count the number of up spins
            
            // Diagonal element
            evol_z_even(even_counter, even_counter) = complex<double>((bond + 
                                                       param_.h*(2*spin - size_)),0.0);
            
            // Update an even pair. Here the first one is always no larger than the second one
            even_parity_[even_counter][0] = i;
            even_parity_[even_counter][1] = pair;
            even_counter++;

            // Update an odd pair
            if(pair != i)
            {
                odd_parity_[odd_counter][0] = i;
                odd_parity_[odd_counter][1] = pair;
                odd_counter++;
            }
        }
    }

    if (debug_){
    	cout << "Even Hz:" << endl;
    	complex_matrix_write(evol_z_even);
    	cout << endl;
    }

    // Compute off-diagonal element
    for(int i=0; i<even_dim_; i++)
    {
        int counter = 0;
        for(int j=i+1; j<even_dim_; j++)
        {
            // Calculate the number of different bits
            counter = __builtin_popcount(even_parity_[i][0] ^ even_parity_[j][0]);
            if(counter == 1)
            {
                evol_x_even(i,j) += complex<double>(param_.g,0.0);
                evol_x_even(j,i) += complex<double>(param_.g,0.0);
            }

            counter = __builtin_popcount(even_parity_[i][0] ^ even_parity_[j][1]);
            if(counter == 1)
            {
                evol_x_even(i,j) += complex<double>(param_.g,0.0);
                evol_x_even(j,i) += complex<double>(param_.g,0.0);
            }

            if(even_parity_[i][0] == even_parity_[i][1])
            {
                evol_x_even(i,j) = evol_x_even(i,j) / complex<double>(sqrt(2),0.0);
                evol_x_even(j,i) = evol_x_even(i,j);
            }
            if(even_parity_[j][0] == even_parity_[j][1])
            {
                evol_x_even(i,j) = evol_x_even(i,j) / complex<double>(sqrt(2),0.0);
                evol_x_even(j,i) = evol_x_even(i,j);
            }
        }
    }

    if (debug_){
    	cout << "Even Hx:" << endl;
    	complex_matrix_write(evol_x_even);
    	cout << endl;
    }
}

/*
 * Construct some odd parts of evolution operator
 */
void FloEvolXXZ::Evol_Odd_Construct_(MatrixXcd& evol_x_odd, MatrixXcd& evol_z_odd){
	// Compute pairs and diagonal elements of hamiltonian
    for(int i = 0; i<odd_dim_; i++)
    {
        int pair = 0; // The pair state, which is the reverse of the current state
        int bond = 0; // Contribution from the nn interaction
       
        int site = odd_parity_[i][0]; // The last binary digit indicates the rightmost spin
        int current_spin = (site & 1); // The rightmost spin
        int left_spin;

        for(int j=0; j<size_-1; j++)
        {
            site = site >> 1; // Proceed to the left site
            left_spin = (site & 1);

            if(current_spin == left_spin)
            {
                bond += 1;
            }
            else
            {
                bond -= 1;
            }

            current_spin = left_spin;
        }

        // Count the # of up spins (# of 1 bits)
        int spin = __builtin_popcount(odd_parity_[i][0]); 
            
        // Diagonal element
        evol_z_odd(i,i) = complex<double>((bond + param_.h*(2*spin - size_)),0);
    }

    if (debug_){
    	cout << "Odd Hz:" << endl;
    	complex_matrix_write(evol_z_odd);
    	cout << endl;
    }

    // Compute off-diagonal element
    for(int i=0; i<odd_dim_; i++)
    {
        int counter = 0;
        for(int j=i+1; j<odd_dim_; j++)
        {
            counter = __builtin_popcount(odd_parity_[i][0] ^ odd_parity_[j][0]);
            if(counter == 1)
            {
                evol_x_odd(i,j) += complex<double>(param_.g,0);
                evol_x_odd(j,i) += complex<double>(param_.g,0);
            }

            counter = __builtin_popcount(odd_parity_[i][0] ^ odd_parity_[j][1]);
            if(counter == 1)
            {
                evol_x_odd(i,j) -= complex<double>(param_.g,0);
                evol_x_odd(j,i) -= complex<double>(param_.g,0);
            }
        }
    }

    if (debug_){
    	cout << "Odd Hz:" << endl;
    	complex_matrix_write(evol_x_odd);
    	cout << endl;
    }
}