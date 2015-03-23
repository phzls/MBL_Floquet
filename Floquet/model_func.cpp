#include <iostream>
#include <algorithm>
#include <string>
#include "model_func.h"
#include "flo_evol_model.h"

using namespace std;

string FloEvolRandomFunc::operator() (const AllPara& parameters, 
EvolMatrix< ComplexEigenSolver<MatrixXcd> >*& model){
	const int size = parameters.generic.size; // System Size
	const double tau = parameters.floquet.tau; // Time step, which seems not used here
	const double J = parameters.floquet.J; // Coupling strength

	const bool debug = parameters.generic.debug;

	model = new FloEvolRandom(size, tau, J, debug);

	string type = model -> Type();
	replace(type.begin(), type.end(), '_', ' ');

	return type;
}

string FloEvolRandomFunc::operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&){
	cout << "Wrong pointers for model." << endl;
	abort();
}

string FloEvolRandomRotationFunc::operator() (const AllPara& parameters, 
EvolMatrix< ComplexEigenSolver<MatrixXcd> >*& model){
	const int size = parameters.generic.size; // System Size
	const double tau = parameters.floquet.tau; // Time step, which seems not used here
	const double J = parameters.floquet.J; // Coupling strength

	const double angle_min = parameters.floquet_random.angle_min; // Minimum angle
	const double angle_sup = parameters.floquet_random.angle_sup; // Supreme angle

	const bool debug = parameters.generic.debug;

	model = new FloEvolRandomRotation(size, tau, J, angle_min, angle_sup, debug);

	string type = model -> Type();
	replace(type.begin(), type.end(), '_', ' ');

	return type;
}

string FloEvolRandomRotationFunc::operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&)
{
	cout << "Wrong pointers for model." << endl;
	abort();
}

string FloEvolXXZFunc::operator() (const AllPara& parameters, 
EvolMatrix< ComplexEigenSolver<MatrixXcd> >*& model){
	const int size = parameters.generic.size; // System Size
	const double tau = parameters.floquet.tau; // Time step, which seems not used here
	const double g = parameters.floquet_xxz.g; // Transverse field strength
	const double h = parameters.floquet_xxz.h; // Longitude field strength

	const bool debug = parameters.generic.debug;

	model = new FloEvolXXZ(size, tau, g, h, debug);

	string type = model -> Type();
	replace(type.begin(), type.end(), '_', ' ');

	return type;
}

string FloEvolXXZFunc::operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&)
{
	cout << "Wrong pointers for model." << endl;
	abort();
}

string FloEvolInterRandomFunc::operator() (const AllPara& parameters, 
EvolMatrix< ComplexEigenSolver<MatrixXcd> >*& model){
	const int size = parameters.generic.size; // System Size
	const double J = parameters.floquet.J; // Coupling strength
	const double tau = parameters.floquet.tau; // Time step, which seems not used here
	const double g = parameters.floquet_xxz.g; // Transverse field strength
	const double h = parameters.floquet_xxz.h; // Longitude field strength

	const bool debug = parameters.generic.debug;

	model = new FloEvolInterRandom(size, tau, J, g, h, debug);

	string type = model -> Type();
	replace(type.begin(), type.end(), '_', ' ');

	return type;
}

string FloEvolInterRandomFunc::operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&)
{
	cout << "Wrong pointers for model." << endl;
	abort();
}

string FloEvolMarkovInterRandomFunc::operator() (const AllPara& parameters, 
EvolMatrix< ComplexEigenSolver<MatrixXcd> >*& model){
	const int size = parameters.generic.size; // System Size
	const double J = parameters.floquet.J; // Coupling strength
	const double tau = parameters.floquet.tau; // Time step, which seems not used here
	const double g = parameters.floquet_xxz.g; // Transverse field strength
	const double h = parameters.floquet_xxz.h; // Longitude field strength

	const bool debug = parameters.generic.debug;

	model = new FloEvolMarkovInterRandom(size, J, tau, g, h, debug);

	string type = model -> Type();
	replace(type.begin(), type.end(), '_', ' ');

	return type;
}

string FloEvolMarkovInterRandomFunc::operator() (const AllPara&, 
EvolMatrix< EigenSolver<MatrixXd> >*&)
{
	cout << "Wrong pointers for model." << endl;
	abort();
}

string FloEvolMarkovInterRandomBothFunc::operator() (const AllPara& parameters, 
EvolMatrix< ComplexEigenSolver<MatrixXcd> >*& model){
	const int size = parameters.generic.size; // System Size
	const double J = parameters.floquet.J; // Coupling strength
	const double tau = parameters.floquet.tau; // Time step, which seems not used here
	const double g = parameters.floquet_xxz.g; // Transverse field strength
	const double h = parameters.floquet_xxz.h; // Longitude field strength

	const bool debug = parameters.generic.debug;

	model = new FloEvolMarkovInterRandomBoth(size, J, tau, g, h, debug);

	string type = model -> Type();
	replace(type.begin(), type.end(), '_', ' ');

	return type;

//	type_ = model -> Type();
//	replace(type_.begin(), type_.end(), '_', ' ');
}

string FloEvolMarkovInterRandomBothFunc::operator() (const AllPara&, 
EvolMatrix< EigenSolver<MatrixXd> >*&)
{
	cout << "Wrong pointers for model." << endl;
	abort();
}

string FloEvolMarkovInterRandomBothXFunc::operator() (const AllPara& parameters, 
EvolMatrix< ComplexEigenSolver<MatrixXcd> >*& model){
	const int size = parameters.generic.size; // System Size
	const double J = parameters.floquet.J; // Coupling strength
	const double tau = parameters.floquet.tau; // Time step, which seems not used here
	const double g = parameters.floquet_xxz.g; // Transverse field strength
	const double h = parameters.floquet_xxz.h; // Longitude field strength
	const double K = parameters.markov.K; // Coupling strength

	const bool iso_keep = parameters.generic.iso_keep; // Keep isolated part
	const bool debug = parameters.generic.debug;

	model = new FloEvolMarkovInterRandomBothX(size, J, tau, g, h, K, iso_keep, debug);

	string type = model -> Type();
	replace(type.begin(), type.end(), '_', ' ');

	return type;
}

string FloEvolMarkovInterRandomBothXFunc::operator() (const AllPara&, 
EvolMatrix< EigenSolver<MatrixXd> >*&)
{
	cout << "Wrong pointers for model." << endl;
	abort();
}

string FloEvolXXZRandomFunc::operator() (const AllPara& parameters,
		EvolMatrix< ComplexEigenSolver<MatrixXcd> >*& model){
	const int size = parameters.generic.size; // System Size
	const double tau = parameters.floquet.tau; // Time step, which seems not used here
	const double g = parameters.floquet_xxz.g; // Transverse field strength
	const double h = parameters.floquet_xxz.h; // Longitude field strength
	const double lambda = parameters.floquet.J; // Disorder strength

	const bool debug = parameters.generic.debug;

	model = new FloEvolXXZRandom(size, tau, g, h, lambda, debug);

	string type = model -> Type();
	replace(type.begin(), type.end(), '_', ' ');

	return type;
}

string FloEvolXXZRandomFunc::operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&)
{
	cout << "Wrong pointers for model." << endl;
	abort();
}

string FloEvolMarkovXXZRandomBothXFunc::operator() (const AllPara& parameters,
		EvolMatrix< ComplexEigenSolver<MatrixXcd> >*& model){
	const int size = parameters.generic.size; // System Size
	const double lambda = parameters.floquet.J; // Disorder strength
	const double tau = parameters.floquet.tau; // Time step, which seems not used here
	const double g = parameters.floquet_xxz.g; // Transverse field strength
	const double h = parameters.floquet_xxz.h; // Longitude field strength
	const double K = parameters.markov.K; // Coupling strength

	const bool iso_keep = parameters.generic.iso_keep; // Keep isolated part
	const bool debug = parameters.generic.debug;

	model = new FloEvolMarkovXXZRandomBothX(size, tau, g, h, lambda, K, iso_keep, debug);

	string type = model -> Type();
	replace(type.begin(), type.end(), '_', ' ');

	return type;
}

string FloEvolMarkovXXZRandomBothXFunc::operator() (const AllPara&,
		EvolMatrix< EigenSolver<MatrixXd> >*&)
{
	cout << "Wrong pointers for model." << endl;
	abort();
}

string FloEvolXXZRandomSimpFunc::operator() (const AllPara& parameters,
		EvolMatrix< ComplexEigenSolver<MatrixXcd> >*& model){
	const int size = parameters.generic.size; // System Size
	const double W = parameters.floquet.J; // Disorder strength

	const bool debug = parameters.generic.debug;

	model = new FloEvolXXZRandomSimp(size, W, debug);

	string type = model -> Type();
	replace(type.begin(), type.end(), '_', ' ');

	return type;
}

string FloEvolXXZRandomSimpFunc::operator() (const AllPara&, EvolMatrix< EigenSolver<MatrixXd> >*&)
{
	cout << "Wrong pointers for model." << endl;
	abort();
}