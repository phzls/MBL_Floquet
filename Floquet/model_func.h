#ifndef MODEL_FUNC_H
#define MODEL_FUNC_H

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
 	protected:
 		string type_;
 	public:
 		virtual void operator() (const AllPara&, 
 			EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&) = 0;

 		virtual void operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&) = 0; 

 		// Get the type of the model: Floquet, Hamiltonian or All
 		string Type() const {return type_;}

 		virtual ~ModelFunc(){};
};

// For random floquet operator
class FloEvolRandomFunc: public ModelFunc
{
	public:
		void operator() (const AllPara&, EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);

		void operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&);

		virtual ~FloEvolRandomFunc(){};
};


// For random rotation floquet operator
class FloEvolRandomRotationFunc: public ModelFunc
{
	public:
		void operator() (const AllPara&, EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);
		void operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&);

		virtual ~FloEvolRandomRotationFunc(){};
};

// For xxz floquet operator
class FloEvolXXZFunc: public ModelFunc
{
	public:
		void operator() (const AllPara&, EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);
		void operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&);

		virtual ~FloEvolXXZFunc(){};
};

// For inter random floquet operator
class FloEvolInterRandomFunc: public ModelFunc
{
	public:
		void operator() (const AllPara&, EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);
		void operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&);

		virtual ~FloEvolInterRandomFunc(){};
};

// For Markov inter random floquet operator
class FloEvolMarkovInterRandomFunc: public ModelFunc
{
	public:
		void operator() (const AllPara&, EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);
		void operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&);

		virtual ~FloEvolMarkovInterRandomFunc(){};
};

// For Markov inter random both floquet operator
class FloEvolMarkovInterRandomBothFunc: public ModelFunc
{
	public:
		void operator() (const AllPara&, EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);
		void operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&);

		virtual ~FloEvolMarkovInterRandomBothFunc(){};
};

#endif