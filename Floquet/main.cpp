#include <iostream>
#include <ctime>
#include <string>
#include "constants.h"
#include "parameters.h"
#include "tasks_models.h"
#include "randomc.h"

using namespace std;

CRandomMersenne RanGen_mersenne(time(NULL));

TasksModels tasks_models; // Record all the tasks and methods

int main(){
	AllPara parameters;
	
	parameters.generic.task = "Flo Simple Markov Evolution One Model";
	parameters.generic.model = "Markov Inter Random Both X Flo";

	parameters.generic.size = 2; // System size
	parameters.generic.num_realizations = 1; // Number of realizations
	parameters.generic.threads_N = 1; // Number of threads in openmp
	parameters.generic.evec = false; // Whether compute eigenvectors, so far only called in
									 // level statistics calculation
	parameters.generic.erase = true; // Whether erase matrix after diagonization
	parameters.generic.debug = true; // Whether output debug information

	parameters.output.width = 30; // Width for spacing in output files
	parameters.output.filename_output = true; // Whether print out file names

	parameters.floquet.J_N = 11; // Number of points of coupling strength
	parameters.floquet.J_min = 0; // Minimum J
	parameters.floquet.J_max = 1; // Maximum J
	parameters.floquet.tau = 0.8; // Time step size
	parameters.floquet.J = 0.3;

	parameters.markov.K = 0.8; // Coupling strength to the bath in markov models

	parameters.floquet_random.angle_min = 5*Pi/6; // Minimum angle for the random rotation angle
	parameters.floquet_random.angle_sup = 5*Pi/6; // Supreme angle for the random rotation angle

	parameters.floquet_xxz.g = 0.9045; // Transverse field strength
	parameters.floquet_xxz.h = 0.8090; // Longitude field strength

	parameters.evolution.time_step = 10; // Number of time steps
	parameters.evolution.step_size = parameters.floquet.tau; // Time step size
	parameters.evolution.init_func_name = "Leftmost Spin Z Value"; // Initial state name
	parameters.evolution.evol_compute["Entropy Per Model"] = false;
	parameters.evolution.evol_compute["Leftmost Spin Z Per Model"] = false;
	parameters.evolution.evol_compute["Leftmost Spin Z One Run"] = true;
	parameters.evolution.model_num = 1; // Number of models for evolution
	// If partition the chain to two halves, the size of left part
	parameters.evolution.left_size = parameters.generic.size / 2; 
	parameters.evolution.jump = 1; // jump of time points in evolution

	parameters.evolution.markov_time_jump = 10; // Time jump in the markov process for flipping
										        // True time: time_step * markov_time_jump
	parameters.evolution.markov_jump = true; // Determine whether there will be markov_time_jump
	if (parameters.generic.task.find("Markov") == string::npos){
		// The task has nothing to do with markov process
		parameters.evolution.markov_jump = false;
	}

	parameters.evolution.log_time = false; // whehter time changes logarithmically
	parameters.evolution.log_time_jump = 2; // The base for time change logarithmically

	// The number gives the index of leftmost spin z value
	parameters.evolution.leftmost_spin_z_index = 3;

	Eigen::initParallel();

	tasks_models.Task(parameters.generic.task)(parameters);

	cout << "Calculation finished." << endl;

	return 0;
}
