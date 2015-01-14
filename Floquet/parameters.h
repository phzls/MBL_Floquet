#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>

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
};

#endif