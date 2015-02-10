#ifndef FLO_EVOL_MODEL_H
#define FLO_EVOL_MODEL_H

#include <sstream>
#include <iostream>
#include <vector>
#include "constants.h"
#include "flo_evol.h"

using namespace std;

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
		const double angle_min; // Minimum angle for rotation on bloch sphere
		const double angle_sup; // Supreme angle for rotation on bloch sphere

		Param(double tau, double J, double angle_min, double angle_sup): 
			tau(tau), J(J), angle_min(angle_min), angle_sup(angle_sup) {};
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

		void Para_Check_(); // Check input parameters for the model

		const bool debug_; // Used for debug outputs

	public:
		FloEvolRandomRotation(int size, double tau, double J):
			FloEvolVanilla(size), param_(tau, J, 0, 2*Pi), debug_(false) 
			{ Para_Check_(); Repr_Init_(); }

		FloEvolRandomRotation(int size, double tau, double J, bool debug):
			FloEvolVanilla(size), param_(tau, J, 0, 2*Pi), debug_(debug) 
			{Para_Check_(); Repr_Init_(); }

		// angle_min gives the minimum angle of the rotation angle. angle_sup gives the
		// supreme angle of rotation.
		FloEvolRandomRotation(int size, double tau, double J, double angle_min, double angle_sup):
			FloEvolVanilla(size), param_(tau, J, angle_min, angle_sup), debug_(false) 
			{ Para_Check_(); Repr_Init_(); }

		FloEvolRandomRotation(int size, double tau, double J, double angle_min, double angle_sup,
		bool debug):
			FloEvolVanilla(size), param_(tau, J, angle_min, angle_sup), debug_(debug) 
			{ Para_Check_(); Repr_Init_(); }

		// Initialize a random angle and a random 3D unit vector
		void Evol_Para_Init();

		// Construct evolutionary operator
		void Evol_Construct();

		virtual ~FloEvolRandomRotation() {};
};

//===============================================================================================
/*
 * xxz floquet operator. The dimension at each site is 2. Its time evolutionary operator is 
 * Ux * Uz, where Ux is the exponential of the transverse external field in x direction, and 
 * Uz is constructed from Ising nearest neighbor interactions with interaction strength 1 and
 * longitude field in z direction. This model takes parameters tau, g and h, where g is field
 * strength in x direction and h is the field strength in z direction.
 */

 class FloEvolXXZ : public FloEvolParity
{
	struct Param // The parameters used in the model
	{
		const double tau; // Time step size
		const double g; // Transverse field strength
		const double h; // Longitude field strength
		const int threads_num; // Number of threads used for evolution operator computation

		Param(double tau, double g, double h, int N = 1): 
			tau(tau), g(g), h(h), threads_num(N) {};
	};

	private:
		const Param param_;

		const int even_dim_; // Dimension of even sector
		const int odd_dim_; // Dimension of odd sector

		// Construct even sector of time evolution operator. Even and odd
		// parity vectors are also constructed in it.
		void Evol_Even_Construct_(MatrixXcd&, MatrixXcd&);

		// Construct odd sector of time evolution operator
		void Evol_Odd_Construct_(MatrixXcd&, MatrixXcd&);

		// A general construction function
		void Evol_General_Construct_(MatrixXcd&, int);

		void Repr_Init_(); // Initialize the representation string stream as well as type

		const bool debug_; // Used for debug outputs

	public:
		FloEvolXXZ(int size, double tau, double g, double h, bool debug = false, int N = 1):
			FloEvolParity(size), param_(tau, g, h, N), debug_(debug),
			even_dim_( (dim_>>1) + (1<<( (size+1)/2-1) ) ), 
  			odd_dim_ ( (dim_>>1) - (1<<( (size+1)/2-1) ) )  { Repr_Init_();}

		// Construct evolutionary operator
		void Evol_Construct(); 

		// No parameter to initialize
		void Evol_Para_Init() {cout <<"No parameter to initialize." <<endl;}

		// Return dimension of the even and odd sectors, with even first.
		vector<int> Get_Sector_Dim() const{
			vector<int> dim(2);
			dim[0] = even_dim_; dim[1] = odd_dim_;
			return dim;
		}

		virtual ~FloEvolXXZ() {};
};

//=================================================================================================

/*
 * The random floquet operator which interpolates flo_evol_xxz (when J=1) and 
 * flo_evol_random rotation (when J=0).
 */
class FloEvolInterRandom : public FloEvolVanilla
{
	struct Param // The parameters used in the model
	{
		const double tau; // Time step size
		const double J; // Nearest neighboring
		const double g; // Transverse field strength
		const double h; // Longitude field strength
		const int size; // Size of the chain

		Param(double tau, double J, double g, double h, int L): 
			tau(tau), J(J), g(g), h(h), size(L) {};
	};

	private:
		const Param param_;

		void Repr_Init_(); // Initialize the representation string stream as well as type

		const bool debug_; // Used for debug output

	public:
		FloEvolInterRandom(int size, double tau, double J, double g, double h, bool debug = false):
			FloEvolVanilla(size), param_(tau, J, g, h, size), debug_(debug) { Repr_Init_(); }

		// No parameters to initialize
		void Evol_Para_Init() {};

		// Construct evolutionary operator
		void Evol_Construct(); 

		virtual ~FloEvolInterRandom() {};
};

//=================================================================================================

/*
 * The random floquet operator which interpolates flo_evol_xxz (when J=1) and 
 * flo_evol_random rotation (when J=0) used in Markov time dynamics. The rightmost spin is
 * coupled with the bath, which corresponding to position 0 in binary representation. The
 * default value takes g=0.9045 and h=0.8090. By default tau=0.8. It has 2 evolution operators,
 * 0 or "down" corresponds to the one where bath spin is down, and 1 or "up" corresponds to the
 * one where bath spin is up. They are written in the basic binary basis.
 */
class FloEvolMarkovInterRandom : public FloEvolMultiSec
{
	struct Param // The parameters used in the model
	{
		const double tau; // Time step size
		const double J; // Nearest neighboring
		const double g; // Transverse field strength
		const double h; // Longitude field strength
		const int size; // Size of the chain

		Param(double tau, double J, double g, double h, int L): 
			tau(tau), J(J), g(g), h(h), size(L) {};
	};

	private:
		const Param param_;

		void Repr_Init_(); // Initialize the representation string stream as well as type

		void Op_Name_Init_(); // Initialize the map linking name and position of evol_op

		const bool debug_; // Used for debug output

		void Bath_XXZ_Construct_(MatrixXcd&, string); // Construct the xxz part of the operator

		void Eigen_Name_Construct_(); // Construct eigen_name during diagonalization

	public:
		FloEvolMarkovInterRandom(int size, double J, double tau = 0.8, double g = 0.9045, 
			double h = 0.8090, bool debug = false):
			FloEvolMultiSec(size, 2), param_(tau, J, g, h, size), debug_(debug) { Repr_Init_(); 
			Op_Name_Init_();}

		// No parameters to initialize
		void Evol_Para_Init() {};

		// Construct evolutionary operator
		void Evol_Construct(); 

		// Return the type of the basis that eigenstates are written in
		string Eigen_Type() const {return "Basic";}

		void Transition_Compute(TransitionMatrix& transition, const string& matrix_name) const{
			cout << "Transition matrix computation is not implemented for " << Type() << endl;
			abort();				
		}

		// Return dimension of each sector. Here each sector has dimension dim_, and the number
		// of sectors equal to number of evol_op
		vector<int> Get_Sector_Dim() const{
			vector<int> dim(evol_op_.size());
			for (int i=0; i< dim.size(); i++) dim[i] = dim_;
			return dim;
		}

		virtual ~FloEvolMarkovInterRandom() {};
};

//=================================================================================================

/*
 * This operator constructs both flo_evol_markov_inter_random and flo_evol_inter_random which
 * shares the same flo_evol_random_rotation part. Therefore the latter is just the "isolated"
 * part of the former. In the evol_op, 0 is for coupling to bath of down spin, 1 for coupling 
 * to bath of up spin, 2 for the isolated system
 */
class FloEvolMarkovInterRandomBoth : public FloEvolMultiSec
{
	struct Param // The parameters used in the model
	{
		const double tau; // Time step size
		const double J; // Nearest neighboring
		const double g; // Transverse field strength
		const double h; // Longitude field strength
		const int size; // Size of the chain

		Param(double tau, double J, double g, double h, int L): 
			tau(tau), J(J), g(g), h(h), size(L) {};
	};

	private:
		const Param param_;

		void Repr_Init_(); // Initialize the representation string stream as well as type

		void Op_Name_Init_(); // Initialize the map linking name and position of evol_op

		const bool debug_; // Used for debug output

		void Bath_XXZ_Construct_(MatrixXcd&, string); 
			// Construct the xxz part of the bath operator

		void Eigen_Name_Construct_(); // Construct eigen_name during diagonalization

	public:
		FloEvolMarkovInterRandomBoth(int size, double J, double tau = 0.8, double g = 0.9045, 
			double h = 0.8090, bool debug = false):
			FloEvolMultiSec(size, 3), param_(tau, J, g, h, size), debug_(debug) { Repr_Init_(); 
			Op_Name_Init_();}

		// No parameters to initialize
		void Evol_Para_Init() {};

		// Construct evolutionary operator
		void Evol_Construct(); 

		// Return the type of the basis that eigenstates are written in
		string Eigen_Type() const {return "Basic";}

		void Transition_Compute(TransitionMatrix& transition, const string& matrix_name) const{
			cout << "Transition matrix computation is not implemented for " << Type() << endl;
			abort();				
		}

		// Return dimension of each sector. Here each sector has dimension dim_, and the number
		// of sectors equal to number of evol_op
		vector<int> Get_Sector_Dim() const{
			vector<int> dim(evol_op_.size());
			for (int i=0; i< dim.size(); i++) dim[i] = dim_;
			return dim;
		}

		virtual ~FloEvolMarkovInterRandomBoth() {}
};



#endif