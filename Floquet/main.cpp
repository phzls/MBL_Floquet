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
	parameters.generic.model = "Inter Random Flo";

	parameters.generic.size = 10; // System size
	parameters.generic.num_realizations = 1; // Number of realizations
	parameters.generic.threads_N = 4; // Number of threads in openmp
	parameters.generic.evec = false; // Whether compute eigenvectors, so far only called in
									 // level statistics calculation
	parameters.generic.erase = true; // Whether erase matrix after diagonization
	parameters.generic.debug = false; // Whether output debug information

	parameters.output.width = 30; // Width for spacing in output files
	parameters.output.filename_output = true; // Whether print out file names

	parameters.floquet.J_N = 11; // Number of points of coupling strength
	parameters.floquet.J_min = 0; // Minimum J
	parameters.floquet.J_max = 1; // Maximum J
	parameters.floquet.tau = 0.8; // Time step size
	parameters.floquet.J = 0.9;

	parameters.floquet_random.angle_min = 5*Pi/6; // Minimum angle for the random rotation angle
	parameters.floquet_random.angle_sup = 5*Pi/6; // Supreme angle for the random rotation angle

	parameters.floquet_xxz.g = 0.9045; // Transverse field strength
	parameters.floquet_xxz.h = 0.8090; // Longitude field strength

	parameters.evolution.time_step = 10; // Number of time steps
	parameters.evolution.step_size = parameters.floquet.tau; // Time step size
	parameters.evolution.init_func_name = "Random Product";
	parameters.evolution.evol_compute["Entropy Per Model"] = true;
	// Initial state name
	parameters.evolution.model_num = 1; // Number of models for evolution
	// If partition the chain to two halves, the size of left part
	parameters.evolution.left_size = parameters.generic.size / 2; 
	parameters.evolution.jump = 1; // jump of time points in evolution

	parameters.evolution.log_time = true; // whehter time changes logarithmically
	parameters.evolution.log_time_jump = 2; // The base for time change logarithmically


	tasks_models.Task(parameters.generic.task)(parameters);

	return 0;
}