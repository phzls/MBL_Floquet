#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Eigenvalues> 
#include "evol_class.h"

using namespace std;
using namespace Eigen;

/*
The evolution operator for Floquet system with no apparent symmetry to reduce
time evolution operator.
*/

class FloEvolVanilla : public EvolMatrix
{		
	protected:
		MatrixXcd evol_op_; // Time evolution operator

	public:
		ComplexEigenSolver<MatrixXcd> eigen; // Eigenvalues and Eigenvectors

		// When local dimension is not given
		FloEvolVanilla(int size): EvolMatrix(size){
			evol_op_ = MatrixXcd::Zero(dim_, dim_);
		}

		// When local dimension is given
		FloEvolVanilla(int size, int local_dim): EvolMatrix(size, local_dim){
			evol_op_ = MatrixXcd::Zero(dim_, dim_);
		}

		// Diagnolize time evolution matrix with eigenvectors kept
		void Evol_Diag() {eigen.compute(evol_op_);} 

		// Diagnolize time evolution matrix, user can determine whether eigenvectors are kept
		// False is not kept; True is kept
		void Evol_Diag(bool keep) {eigen.compute(evol_op_, keep);} 

		virtual ~FloEvolVanilla() {};
};

/*
The random floquet operator. The dimension at each site is 2. Its time evolutionary operator 
is UxUz, where Ux has uniformly random 2*2 unitary matrix at each site, and Uz is constructed
from Ising nearest neighbor interactions.
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

	public:
		FloEvolRandom(int size, double tau, double J):
			FloEvolVanilla(size), param_(tau, J) {};

		void Evol_Construct(); // Construct evolutionary operator

		virtual ~FloEvolRandom() {};
};