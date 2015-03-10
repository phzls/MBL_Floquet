#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <vector>
#include <utility>
#include <map>

using namespace std;

/**
 ** This file constains structs that give all necessary parameters to pass to various 
 ** functions. They may contain more parameters than needed for generality.
 **/

/*
 * Generic parameters which may be used by any method/function
 */
struct GenericPara
{
	int num_realizations; // Number of realizations
	int threads_N; // Number of threads 
	int size; // System size
	int local_dimension; // Local dimension at each site
	bool evec; // Whether eigenvectors are computed during diagonization
	bool erase; // Whether erase the Hamiltonian/Evolution Operator after diagonization
	bool debug; // Whether show debug information
	string task; // String for the computation task
	string model; // String for the model
	bool iso_keep; // Whether keep the isolated part; so far only works for 
				   // Flo_Evol_Markov_Inter_Random_Both_X
};

/*
 * Parameters that are related to output and writing to files
 */
struct OutputPara
{
	int width; // Width of output in files
	bool filename_output; // Whether write out filenames
};

/*
 * Parameters that are related to Floquet models.
 */
struct FloPara
{
	double J_min; // Minimum J in a loop
	double J_max; // Maximum J in a loop
	int J_N; // Number points for different coupling strength J, including J_min
			 // and J_max
	double J; // Coupling strength J
	double tau; // Time step
};

/*
 * Parameters that are related to Random Rotation Floquet models.
 */
struct FloRrotationPara
{
	double angle_min; // Minimum angle in rotation. Should be between [0, 2*pi)
	double angle_sup; // Supreme angle in rotation. Should be between [0, 2*pi), and no
					  // less than angle_min
};

/*
 * Parameters that are related to Random Rotation Floquet models.
 */
struct FloXXZPara
{
	double g; // Transverse field strength
	double h; // Longitude field strength
};

struct MatrixPara
{
	string type; // Determine the representation type of the matrix
};

/*
 * Parameters used in time evolutions
 */
struct Evolution
{
	int time_step; // Number of time steps
	double step_size; // Size of time step; would be tau for floquet models
	int model_num; // Number of different models used for time evolution
	int jump; // jump of time points in evolution
	string init_func_name; // Initial state construction function name
	map<string,bool> evol_compute; // Determine which data to compute
	map<string,bool> evol_total_compute; // Determine which data to compute which contains
										 // results for all models
	int left_size; // If partition the chain to two halves, the size of left part

	bool log_time; // Determine whether time increases logarithmically
	int log_time_jump; // The base number which time increases on when logarithmically

	bool markov_jump; // Determine whether there will be markov_time_jump
	int markov_time_jump; // The jump time in markov time evolution

	int leftmost_spin_z_index; // The number gives the index of leftmost spin z value

};

/*
 * Parameters used for Markov models
 */
struct Markov
{
	double K; // Coupling strength to the bath
};

/*
 * Initial information related to multiple sets of parameters under one model of evolution
 * Since so far each initial state construction function only requires one set of parameters,
 * this vector should not be used in any of these functions. Instead, the corresponding one-
 * set-parameter variable should be passed in.
 */
struct MultipleInitPara
{
	// The set of numbers which give indices of leftmost spin z value
	vector<int> leftmost_spin_z_index_set;

	// Whether output full evolution results or just some of them from leftmost_spin_z_index_set.
	// If it is true, then only some are output, and whether they change signs are outputted
	bool easy_full_leftmost_spin_z;

	// A small value beyond which is considered as non-zero
	double non_zero_threshold;
};

/*
 * Parameters used for single model computation
 */
struct SingleModel
{
	map<string,bool> single_model_compute; // Determine which method to be computed
};

/*
 * All parameters.
 */
struct AllPara
{
	// Generic parameters
	GenericPara generic;

	// Output parameters
	OutputPara output;

	// Floquet parameters
	FloPara floquet;

	// Random Rotation Floquet parameters
	FloRrotationPara floquet_random;

	// XXZ Floquet parameters
	FloXXZPara floquet_xxz;

	// Matrix relevant parameters
	MatrixPara matrix_para;

	// Time evolution parameters
	Evolution evolution;

	// Markov model parameters
	Markov markov;

	// Multiple sets of initial conditions
	MultipleInitPara multi_ini_para;

	// Single model computation
	SingleModel single_model;
};

#endif