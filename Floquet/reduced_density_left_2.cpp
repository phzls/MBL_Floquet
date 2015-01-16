#include <Eigen/Dense>
#include <complex>
#include <iostream>

using namespace std;
using namespace Eigen;

/**
 ** This file calculates the reduced density matrix for the left part of a spin chain whose
 ** local dimension at each site is 2. It is constructed from state_basic, a vector written
 ** in basic binary basis. size is the total chain length and left_size is the length of the
 ** left part. The binary representation of binary basis states has 0 being down and 1 being
 ** up, and the counting starts from the rightmost spin.
 **/

void reduced_density_left_2(const VectorXcd& state_basic, int size, int left_size, 
MatrixXcd& reduced_density){
	  const int right_size = size - left_size; // Length of right part of the chain
	  const int left_dim = 1 << left_size; // Total dimension of the left part Hiblert space
    const int right_dim = 1 << right_size; // Total dimension of the right part space

  	if (right_size <=0){
  		cout << "Left part for reduced density matrix calculation is too long." << endl;
  		abort();
  	}

  	reduced_density = MatrixXcd::Zero(left_dim, left_dim);

  	/*compute reduced density matrix*/
  	for(int i=0; i<left_dim; i++)
  	{
    	for(int j=0; j<=i; j++)
	  	{
	  		for(int k=0; k<right_dim; k++)
	    	{
	      		reduced_density(j,i) += state_basic(j*right_dim + k) * conj(
	      								state_basic(i*right_dim + k) );
	    	}
	  		if(i != j) reduced_density(i,j) = conj(reduced_density(j,i));
	  	}
  	}
}