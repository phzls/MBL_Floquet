#include <iostream>
#include <cmath>
#include <omp.h>
#include "transition.h"

using namespace std;
using namespace Eigen;

/**
 ** This file implements the construction of transition matrix from full chain parity
 ** eigenstates to parity states. The eigenstates are assumed to be stored column-wise.
 **/

void TransitionMatrix::Parity_Full(const MatrixXd& even_evec, const MatrixXd& odd_evec, 
int threads_N){
	const int even_rank = even_evec.cols();
	const int odd_rank = odd_evec.cols();

	if (Check_Matrix("Even_Full")){
		cout << "Transition matrix even_full has been constructed." << endl;
		abort();
	}

	if (Check_Matrix("Odd_Full")){
		cout << "Transition matrix odd_full has been constructed." << endl;
		abort();
	}

	// Initialize the two matrices
	even_full_ = MatrixXcd::Zero(even_rank, even_rank);
	odd_full_ = MatrixXcd::Zero(odd_rank, odd_rank);

	// Compute the elements of even_full and odd_full
	#pragma omp parallel num_threads(threads_N)
	{
		#pragma omp for

		for (int i=0;i<even_rank;i++)
		{
			for (int j=0;j<even_rank;j++)
			{
				even_full_(i,j) = complex<double>(even_evec(i,j),0);
			}
		}

		for (int i=0;i<odd_rank;i++)
		{
			for (int j=0;j<odd_rank;j++)
			{
				odd_full_(i,j) = complex<double>(odd_evec(i,j),0);
			}
		}
	}

	// Add the two matrices to constructed map
	constructed_type_["Even_Full"] = &even_full_;
	constructed_type_["Odd_Full"] = &odd_full_;
}

void TransitionMatrix::Parity_Full(const MatrixXcd& even_evec, const MatrixXcd& odd_evec, 
int threads_N){
	const int even_rank = even_evec.cols();
	const int odd_rank = odd_evec.cols();

	if (Check_Matrix("Even_Full")){
		cout << "Transition matrix even_full has been constructed." << endl;
		abort();
	}

	if (Check_Matrix("Odd_Full")){
		cout << "Transition matrix odd_full has been constructed." << endl;
		abort();
	}

	// Initialize the two matrices
	even_full_ = MatrixXcd::Zero(even_rank, even_rank);
	odd_full_ = MatrixXcd::Zero(odd_rank, odd_rank);

	// Compute the elements of even_full and odd_full
	#pragma omp parallel num_threads(threads_N)
	{
		#pragma omp for

		for (int i=0;i<even_rank;i++)
		{
			for (int j=0;j<even_rank;j++)
			{
				even_full_(i,j) = even_evec(i,j);
			}
		}

		for (int i=0;i<odd_rank;i++)
		{
			for (int j=0;j<odd_rank;j++)
			{
				odd_full_(i,j) = odd_evec(i,j);
			}
		}
	}

	// Add the two matrices to constructed map
	constructed_type_["Even_Full"] = &even_full_;
	constructed_type_["Odd_Full"] = &odd_full_;
}