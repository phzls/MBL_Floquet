#include <iostream>
#include "model_func.h"
#include "flo_evol_model.h"

using namespace std;

void FloEvolRandomFunc::operator() (const AllPara& parameters, 
EvolMatrix< ComplexEigenSolver<MatrixXcd> >*& model){
	const int size = parameters.generic.size; // System Size
	const double tau = parameters.floquet.tau; // Time step, which seems not used here
	const int J = parameters.floquet.J; // Coupling strength

	model = new FloEvolRandom(size, tau, J);
}

void FloEvolRandomFunc::operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&){
	cout << "Wrong pointers for model." << endl;
	abort();
}

void FloEvolRandomRotationFunc::operator() (const AllPara& parameters, 
EvolMatrix< ComplexEigenSolver<MatrixXcd> >*& model){
	const int size = parameters.generic.size; // System Size
	const double tau = parameters.floquet.tau; // Time step, which seems not used here
	const int J = parameters.floquet.J; // Coupling strength

	const double angle_min = parameters.floquet_random.angle_min; // Minimum angle
	const double angle_sup = parameters.floquet_random.angle_sup; // Supreme angle

	model = new FloEvolRandomRotation(size, tau, J, angle_min, angle_sup);
}

void FloEvolRandomRotationFunc::operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&)
{
	cout << "Wrong pointers for model." << endl;
	abort();
}