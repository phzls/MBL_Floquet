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
}

// Pointer to all possible initial state construction function
typedef void (*init_func)(const InitInfo&, const TransitionMatrix&, VectorXcd&);

class InitObj
{
	private:
		// Map that stores init_func and its name to be called
		map<string, init_func>init_func_map_; 

		void map_init_(); // Initialize init_func_map

	public:
		InitObj() {map_init_();}

		// Use a string to access different initial construction functions
		init_func Init_Func(const string&) const;

		// Print out all init_func
		void Print() const;

		~InitObj(){};
};

/*
 * Initial state functions
 */

void product_random(const InitInfo&, const TransitionMatrix&, VectorXcd&);
void random_product(const InitInfo&, const TransitionMatrix&, VectorXcd&);


/*
 * Some useful functions for initial state construction
 */
void random_pure_amplitude(vector<complex<double> >&);
void norm_check(const VectorXcd&, double, const string&);

#endif