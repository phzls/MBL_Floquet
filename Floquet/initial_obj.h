#ifndef INITIAL_OBJ_H
#define INITIAL_OBJ_H

#include <map>
#include <Eigen/Dense>
#include "transition.h"

using namespace std;
using namespace Eigen;

/**
 ** This file includes classes and functions used for constructing initial states.
 **/

/*
 * Information that is related to constructing initial states. Every initial state construction
 * function just takes this structure, transition matrix class and the initial state vector.
 */

struct InitInfo
{
	int size; // System size
	int dim; // Total dimension of Hilbert space
	double norm_delta; // A small number used to check whether norm is 1
	bool debug; // Whether output debug information
};

// Pointer to all possible initial state construction function which gives a state vector
typedef void (*init_func)(const InitInfo&, const TransitionMatrix&, VectorXcd&);

// Pointer to all possible initial state construction function which gives a density matrix
typedef void (*init_func_C)(const InitInfo&, MatrixXcd&);

class InitObj
{
	private:
		// Map that stores init_func and its name to be called
		map<string, init_func>init_func_map_; 

		// Map that stores init_func_C and its name to be called
		map<string, init_func_C>init_func_C_map_; 

		void map_init_(); // Initialize init_func_map and init_func_C_map

	public:
		InitObj() {map_init_();}

		// Use a string to access different initial construction functions for state vector
		init_func Init_Func(const string&) const;

		// Use a string to access different initial construction functions for density matrix
		init_func_C Init_Func_C(const string&) const;

		// Print out all init_func
		void Print() const;

		// Print out all init_func_C
		void Print_C() const;

		~InitObj(){};
};

/*
 * Initial state functions
 */

void product_random(const InitInfo&, const TransitionMatrix&, VectorXcd&);

void random_product(const InitInfo&, const TransitionMatrix&, VectorXcd&);
void random_product(const InitInfo&, MatrixXcd&);

void random_pure(const InitInfo&, MatrixXcd&);


/*
 * Some useful functions for initial state construction
 */
void random_pure_amplitude(vector<complex<double> >&);
void norm_check(const VectorXcd&, double, const string&);
void state_to_density(const VectorXcd&, MatrixXcd&);

#endif