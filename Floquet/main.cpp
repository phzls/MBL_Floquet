#include <iostream>
#include <ctime>
#include "parameters.h"
#include "mtrand.h"
#include "task.h"

using namespace std;

MTRand u1rand(time(NULL));
MTRand_closed u2rand(time(NULL));

int main(){
	AllPara parameters;
	const bool Flo = true;
	const bool Level = true;

	parameters.generic.size = 6; // System size
	parameters.generic.num_realizations = 4; // Number of realizations
	parameters.generic.threads_N = 4; // Number of threads in openmp
	parameters.generic.evec = false; // Whether compute eigenvectors
	parameters.generic.erase = true; // Whether erase matrix after diagonization

	parameters.output.width = 15; // Width for spacing in output files
	parameters.output.filename_output = true; // Whether print out file names

	parameters.floquet.J_N = 11; // Number of points of coupling strength
	parameters.floquet.J_min = 0; // Minimum J
	parameters.floquet.J_max = 1; // Maximum J
	parameters.floquet.tau = 0.8; // Time step size

	if (Flo){
		if (Level){
			flo_level(parameters);
		}
	}

	return 0;
}