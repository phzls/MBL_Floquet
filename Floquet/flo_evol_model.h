#include <sstream>
#include <vector>
#include "flo_evol.h"

/**
This file contains particular models of Floquet time evolution operators
**/


/*
The random floquet operator. The dimension at each site is 2. Its time evolutionary operator 
is Ux * Uz, where Ux has uniformly random 2*2 unitary matrix at each site, and Uz is 
constructed from Ising nearest neighbor interactions with interaction strength J.
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

		stringstream repr_; // Representation string stream of the model
		void Repr_Init_(); // Initialize the representation string stream

		// Construct Ux part in time evolution operator
		void Evol_Site_Construct_(MatrixXcd&);

		// Construct Uz part in time evolution operator
		void Evol_Z_Construct_(MatrixXcd&);

	public:
		FloEvolRandom(int size, double tau, double J):
			FloEvolVanilla(size), param_(tau, J) { Repr_Init_(); }

		// Initialize four random angles
		void Evol_Para_Init();

		// Construct evolutionary operator
		void Evol_Construct(); 

		// Return the string format of representation string stream 
		string Repr() const {return repr_.str();} 

		virtual ~FloEvolRandom() {};
};