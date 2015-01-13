#include <complex>
#include <cmath>
#include <omp.h>
#include <map>
#include "transition.h"

using namespace std;

void TransitionMatrix::Basic_Parity(const vector<vector<int> >& even_parity, 
const vector<vector<int> >& odd_parity, int threads_N ){
	const int even_rank = even_parity.size();
	const int odd_rank = odd_parity.size();
	const int total_rank = even_rank + odd_rank;

	if (Check_Matrix("Basic_Even")){
		cout << "Transition matrix basic_even has been constructed." << endl;
		abort();
	}

	if (Check_Matrix("Basic_Odd")){
		cout << "Transition matrix basic_odd has been constructed." << endl;
		abort();
	}

	// Initialize the two matrices
	basic_even_.resize(total_rank, even_rank);
	basic_odd_.resize(total_rank, odd_rank);
	
	// Compute the elements of basic_even and basic_odd
	#pragma omp parallel num_threads(threads_N)
	{
		#pragma omp for

		for (int i=0;i<even_rank;i++)
		{
			int row1 = even_parity[i][0];
			int row2 = even_parity[i][1];

			if (row1 == row2)
			{
				basic_even_(row1,i) = complex<double>(1.0,0);
			} 
			else
			{
				basic_even_(row1,i) = complex<double>(1.0/sqrt(2),0);
				basic_even_(row2,i) = complex<double>(1.0/sqrt(2),0);
			}
		}

		for (int i=0;i<odd_rank;i++)
		{
			int j = odd_parity[i][0];
			basic_odd_(j,i) = complex<double>(1.0/sqrt(2),0);

			j = odd_parity[i][1];
			basic_odd_(j,i) = complex<double>(-1.0/sqrt(2),0);
		}
	}

	// Add the matrices to constructed map
	constructed_type_["Basic_Even"] = &basic_even_;
	constructed_type_["Basic_Odd"] = &basic_odd_;
}