#include <iostream>
#include <ctime>
#include "constants.h"
#include "parameters.h"
#include "mtrand.h"
#include "tasks_models.h"

using namespace std;

MTRand u1rand(time(NULL));
MTRand_closed u2rand(time(NULL));

TasksModels tasks_models; // Record all the tasks and methods

int main(){
	AllPara parameters;
	
	parameters.generic.task = "Flo Rightmost Sigma_z";
	parameters.generic.model = "Random Rotation Flo";

	parameters.generic.size = 6; // System size
	parameters.generic.num_realizations = 100; // Number of realizations
	parameters.generic.threads_N = 1; // Number of threads in openmp
	parameters.generic.evec = false; // Whether compute eigenvectors
	parameters.generic.erase = true; // Whether erase matrix after diagonization
	parameters.generic.debug = false; // Whether output debug information

	parameters.output.width = 30; // Width for spacing in output files
	parameters.output.filename_output = true; // Whether print out file names

	parameters.floquet.J_N = 11; // Number of points of coupling strength
	parameters.floquet.J_min = 0; // Minimum J
	parameters.floquet.J_max = 1; // Maximum J
	parameters.floquet.tau = 0.8; // Time step size
	parameters.floquet.J = 0.9;

	parameters.floquet_random.angle_min = 2*Pi/3; // Minimum angle for the random rotation angle
	parameters.floquet_random.angle_sup = 2*Pi/3; // Supreme angle for the random rotation angle

	tasks_models.Task(parameters.generic.task)(parameters);

	return 0;
}