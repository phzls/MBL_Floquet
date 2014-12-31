#ifndef FLO_EVOL_MODEL_H
#define FLO_EVOL_MODEL_H

#include <sstream>
#include <vector>
#include "flo_evol.h"

/**
 ** This file contains particular models of Floquet time evolution operators
 **/


/*
 * The random floquet operator. The dimension at each site is 2. Its time 
 * evolutionary operator is Ux * Uz, where Ux has uniformly random 2*2 unitary matrix 
 * at each site, and Uz is constructed from Ising nearest neighbor interactions with
 * interaction strength J.
 */
class FloEvolRandom : public FloEvolVanilla
{
	struct Param // The parameters used in the model
	{
		const double tau; // Time step size
		const double J; // Nearest neighboring

		Param(double tau, double J): tau(tau), J(J) {};
	};

	private:
		const Param param_;

		// The four random angles used in sampling random SU(2) matrices
		vector<vector<double> > su2_angle_; 

		// Construct Ux part in time evolution operator
		void Evol_Site_Construct_(MatrixXcd&);

		// Construct Uz part in time evolution operator
		void Evol_Z_Construct_(MatrixXcd&);

		void Repr_Init_(); // Initialize the representation string stream as well as type

		const bool debug_; // Used for debug output

	public:
		FloEvolRandom(int size, double tau, double J):
			FloEvolVanilla(size), param_(tau, J), debug_(false) { Repr_Init_(); }

		FloEvolRandom(int size, double tau, double J, bool debug):
			FloEvolVanilla(size), param_(tau, J), debug_(debug) { Repr_Init_(); }

		// Initialize four random angles
		void Evol_Para_Init();

		// Construct evolutionary operator
		void Evol_Construct(); 

		virtual ~FloEvolRandom() {};
};


//===============================================================================================


/*
 * The random rotation floquet operator. The dimension at each site is 2. Its time 
 * evolutionary operator is Ux * Uz, where Ux has 2*2 unitary matrix which randomly
 * rotates any state on a bloch sphere at each site, and Uz is constructed from Ising 
 * nearest neighbor interactions with interaction strength J.
 */
class FloEvolRandomRotation : public FloEvolVanilla
{
	struct Param // The parameters used in the model
	{
		const double tau; // Time step size
		const double J; // Nearest neighboring

		Param(double tau, double J): tau(tau), J(J) {};
	};

	private:
		const Param param_;

		// The random unit vector for rotation axis
		vector<vector<double> > axis_;

		// The random angle for rotation about rotation axis
		vector<double> angle_; 

		// Construct Ux part in time evolution operator
		void Evol_Site_Construct_(MatrixXcd&);

		// Construct Uz part in time evolution operator
		void Evol_Z_Construct_(MatrixXcd&);

		void Repr_Init_(); // Initialize the representation string stream as well as type

		const bool debug_; // Used for debug output

	public:
		FloEvolRandomRotation(int size, double tau, double J):
			FloEvolVanilla(size), param_(tau, J), debug_(false) { Repr_Init_(); }

		FloEvolRandomRotation(int size, double tau, double J, bool debug):
			FloEvolVanilla(size), param_(tau, J), debug_(debug) { Repr_Init_(); }

		// Initialize a random angle and a random 3D unit vector
		void Evol_Para_Init();

		// Construct evolutionary operator
		void Evol_Construct(); 

		virtual ~FloEvolRandomRotation() {};
};

#endif