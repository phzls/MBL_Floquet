#include <iostream>
#include <algorithm>
#include <string>
#include "model_func.h"
#include "flo_evol_model.h"

using namespace std;

void FloEvolRandomFunc::operator() (const AllPara& parameters, 
EvolMatrix< ComplexEigenSolver<MatrixXcd> >*& model){
	const int size = parameters.generic.size; // System Size
	const double tau = parameters.floquet.tau; // Time step, which seems not used here
	const double J = parameters.floquet.J; // Coupling strength

	const bool debug = parameters.generic.debug;

	model = new FloEvolRandom(size, tau, J, debug);

	type_ = model -> Type();
	replace(type_.begin(), type_.end(), '_', ' ');
}

void FloEvolRandomFunc::operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&){
	cout << "Wrong pointers for model." << endl;
	abort();
}

void FloEvolRandomRotationFunc::operator() (const AllPara& parameters, 
EvolMatrix< ComplexEigenSolver<MatrixXcd> >*& model){
	const int size = parameters.generic.size; // System Size
	const double tau = parameters.floquet.tau; // Time step, which seems not used here
	const double J = parameters.floquet.J; // Coupling strength

	const double angle_min = parameters.floquet_random.angle_min; // Minimum angle
	const double angle_sup = parameters.floquet_random.angle_sup; // Supreme angle

	const bool debug = parameters.generic.debug;

	model = new FloEvolRandomRotation(size, tau, J, angle_min, angle_sup, debug);

	type_ = model -> Type();
	replace(type_.begin(), type_.end(), '_', ' ');
}

void FloEvolRandomRotationFunc::operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&)
{
	cout << "Wrong pointers for model." << endl;
	abort();
}

void FloEvolXXZFunc::operator() (const AllPara& parameters, 
EvolMatrix< ComplexEigenSolver<MatrixXcd> >*& model){
	const int size = parameters.generic.size; // System Size
	const double tau = parameters.floquet.tau; // Time step, which seems not used here
	const double g = parameters.floquet_xxz.g; // Transverse field strength
	const double h = parameters.floquet_xxz.h; // Longitude field strength

	const bool debug = parameters.generic.debug;

	model = new FloEvolXXZ(size, tau, g, h, debug);

	type_ = model -> Type();
	replace(type_.begin(), type_.end(), '_', ' ');
}

void FloEvolXXZFunc::operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&)
{
	cout << "Wrong pointers for model." << endl;
	abort();
}