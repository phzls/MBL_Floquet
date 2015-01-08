#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include "evol_class.h"
#include "parameters.h"

using namespace std;
using namespace Eigen;

/**
 ** This file contains the base class and derived classes for generating a pointer to
 ** various evol_class models. To circumvent the template in evol_class, the, up to now,
 ** two possible template specializations are enumerated. A cleverer solution may be used in
 ** the future.
 **/

class ModelFunc
{
 	public:
 		virtual void operator() (const string&, const AllPara&, 
 			EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&) = 0;

 		virtual void operator() (const string&, const AllPara&, 
 			EvolMatrix< EigenSolver<MatrixXd> >*&) = 0; 

 		// Get the type of the model: Floquet, Hamiltonian or All
 		virtual string Get_Type() const = 0; 
};

// For random floquet operator
class FloEvolRandomFunc: public ModelFunc
{
	private:
		const string type_;
	public:
		FloEvolRandomFunc(): type_("Floquet") {};

		void operator() (const string&, const AllPara&, 
 			EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);

		string Get_Type() const;
};


// For random rotation floquet operator
class FloEvolRandomRotationFunc: public ModelFunc
{
	private:
		const string type_;
	public:
		FloEvolRandomRotationFunc(): type_("Floquet") {};

		void operator() (const string&, const AllPara&, 
 			EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);

		string Get_Type() const;
};