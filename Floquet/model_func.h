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
 	public:
 		virtual string operator() (const AllPara&, 
 			EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&) = 0;

 		virtual string operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&) = 0;

 		virtual ~ModelFunc(){};
};

// For random floquet operator
class FloEvolRandomFunc: public ModelFunc
{
	public:
		string operator() (const AllPara&, EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);

		string operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&);

		virtual ~FloEvolRandomFunc(){};
};


// For random rotation floquet operator
class FloEvolRandomRotationFunc: public ModelFunc
{
	public:
		string operator() (const AllPara&, EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);
		string operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&);

		virtual ~FloEvolRandomRotationFunc(){};
};

// For xxz floquet operator
class FloEvolXXZFunc: public ModelFunc
{
	public:
		string operator() (const AllPara&, EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);
		string operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&);

		virtual ~FloEvolXXZFunc(){};
};

// For inter random floquet operator
class FloEvolInterRandomFunc: public ModelFunc
{
	public:
		string operator() (const AllPara&, EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);
		string operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&);

		virtual ~FloEvolInterRandomFunc(){};
};

// For Markov inter random floquet operator
class FloEvolMarkovInterRandomFunc: public ModelFunc
{
	public:
		string operator() (const AllPara&, EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);
		string operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&);

		virtual ~FloEvolMarkovInterRandomFunc(){};
};

// For Markov inter random both floquet operator
class FloEvolMarkovInterRandomBothFunc: public ModelFunc
{
	public:
		string operator() (const AllPara&, EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);
		string operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&);

		virtual ~FloEvolMarkovInterRandomBothFunc(){};
};

// For Markov inter random both floquet operator which couples to bath through a sigma_x term
class FloEvolMarkovInterRandomBothXFunc: public ModelFunc
{
	public:
		string operator() (const AllPara&, EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);
		string operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&);

		virtual ~FloEvolMarkovInterRandomBothXFunc(){};
};

// For random xxz floquet operator
class FloEvolXXZRandomFunc: public ModelFunc
{
public:
	string operator() (const AllPara&, EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);
	string operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&);

	virtual ~FloEvolXXZRandomFunc(){};
};

// For Markov XXZ random both floquet operator which couples to bath through a sigma_x term
class FloEvolMarkovXXZRandomBothXFunc: public ModelFunc
{
public:
	string operator() (const AllPara&, EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);
	string operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&);

	virtual ~FloEvolMarkovXXZRandomBothXFunc(){};
};

// For xxz random simple floquet operator
class FloEvolXXZRandomSimpFunc: public ModelFunc
{
public:
	string operator() (const AllPara&, EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);
	string operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&);

	virtual ~FloEvolXXZRandomSimpFunc(){};
};

// For xxz random simple shift floquet operator
class FloEvolXXZRandomSimpShiftFunc: public ModelFunc
{
public:
	string operator() (const AllPara&, EvolMatrix< ComplexEigenSolver<MatrixXcd> >*&);
	string operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&);

	virtual ~FloEvolXXZRandomSimpShiftFunc(){};
};

#endif