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
		stringstream repr_; // Representation string stream of the model
		void Repr_Init_(); // Initialize the representation string stream

	public:
		FloEvolRandom(int size, double tau, double J):
			FloEvolVanilla(size), param_(tau, J) { Repr_Init_(); }

		// Construct evolutionary operator
		void Evol_Construct(const ModelConstructData<class T>&); 

		// Return the string format of representation string stream 
		string Repr() const {return repr_.str();} 

		virtual ~FloEvolRandom() {};
};

/*
The data structure that is used for this model. 
*/
struct RandomFloData : public ModelConstructData< vector<double> >
{
	void Initialize(); // Generate four random numbers
};