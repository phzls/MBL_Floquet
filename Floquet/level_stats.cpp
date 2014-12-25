#include <fstream>
#include <iostream>
#include <complex>
#include <cmath>
#include "level_stats.h"
#include "sort.h"

using namespace std;

void VanillaFloLevel::Single_Data_Process_(const 
EvolMatrix< ComplexEigenSolver<MatrixXcd> >* const U){

	// The total dimension of Hilbert space
	const int dim = U -> GetDim();

	// The average level spacing of phases of eigenvalues
	const double ave_spacing = 2*pi / dim;

	if (dim != U -> eigen.eigenvalues().rows()){
		cout << "The dimension of eigenvalues of time evolution operator is not correct."<<endl;
		abort();
	}

	if (level_.size() == 0){
		level_.resize(dim);
		for (int i=0; i<dim; i++) level_[i] = 0;
	}
	else if (level_.size() != dim){
		cout << "The dimension of different realizations of eigenvalues do not match."<<endl;
		abort();
	}

	vector<double> phases(dim); // Phases of all eigenvalues in ascending order

	for (int i=0; i<dim; i++){
		phases[i] = arg(U -> eigen.egenvalues()(i));
	}

	// Sort the phases in ascending order
	sort(phases.begin(), phases.end(), Vec_Double_Sort);

	// Compute level statistics
	double local_mean = 0;
	double square_sum = 0;

	double spacing;

	for (int i=0; i<dim-1; i++){
		spacing = ( phases[i+1] - phases[i] ) / ave_spacing;
		level_[i] += spacing;
		sum += spacing;
		square_sum += spacing * spacing;
	}

	
	spacing = ( phases[0] + 2*pi - phases[dim-1] ) / ave_spacing;
	level_[dim-1] += spacing; 
	sum += spacing;
	square_sum += spacing * spacing;

	return pair<double, double>(sum/dim, square_sum/dim);
}

void VanillaFloLevel::Data_Process(const EvolMatrix< ComplexEigenSolver<MatrixXcd> >& U){
	const int num_realization = U.size():

	if (init_){
		cout<< "The object is already holding one set of data." <<endl;
		abort();
	}

	for (int i=0; i<num_realization; i++){

		pair<double, double> mean_square_mean;
		mean_square_mean = Single_Data_Process_(U[i]);

		mean_ += mean_square_mean.first;
		mean_sd_ += mean_square_mean.first * mean_square_mean.first;

		square_mean_ += mean_square_mean.second;
		square_mean_sd_ += mean_square_mean.second * mean_square_mean.second;
	}

	if (num_realization != 0 && level_.size() == 0){
		cout <<"Level spacing calculation fails." << endl;
		abort();
	}

	for (int i=0; i<level_.size();i++) level_[i] /= num_realization;

	mean_ /= num_realization;
	mean_sd_ = sqrt( mean_sd_ / num_realization - mean_ * mean_ );

	square_mean_ /= num_realization;
	square_mean_sd_ = sqrt( square_mean_sd_ / num_realization - square_mean_ * square_mean_ );

	// Construct base filename stringstream
	if (num_realization > 0) base_filename_ << U[0] -> Repr();
	for (int i=1; i< num_realization; ++){
		if (base_filename.str() != U[i] -> Repr()){
			cout <<"String Representations of models are not consistent!"<<endl;
			abort():
		}
	}
	base_filename_ <<",Realizations="<<num_realization;

	init_ = true;
}

void VanillaFloLevel::DataOutput(bool output) const{
	if (!init_){
		cout<< "The object has no data." <<endl;
		abort();
	}

	if (!redirect_ && level_out_){
		stringstream level_stream = base_filename.str() << ",level_spacing.txt";
		ofstream level_output(level_stream.str().c_str());

		if (output) cout << level_stream.str().c_str() << endl;

		for(int i=0; i<level_.size();i++){
			level_output << level_[i]<<endl;
		} 
	}
	if (!redirect_ && mean_out_){
		stringstream mean_stream = base_filename.str() << ",level_spacing_mean.txt";
		ofstream mean_output(mean_stream.str().c_str());

		if (output) cout << mean_stream.str().c_str() << endl;

		mean_output << mean_ <<"  "<< mean_sd_ <<endl;
	}
	if (!redirect_ && square_mean_out_){
		stringstream square_mean_stream = base_filename.str() << ",level_spacing_square_mean.txt";
		ofstream square_mean_output(square_mean_stream.str().c_str());

		if (output) cout << square_mean_stream.str().c_str() << endl;

		square_mean_output << square_mean_ <<"  "<< square_mean_sd_ <<endl;
	}
}

void VanillaFloLevel::DataRedirect(pair< vector<double>, vector<double> >& A) const{
	if (!init_){
		cout << "The object has no data." <<endl;
		abort();
	}
	if (level_.size() == 0){
		cout<< "Level spacing is empty." <<endl;
		abort();
	}

	A.first = level_;
	
	A.second.resize(4);

	A.second[0] = mean_;
	A.second[1] = mean_sd_;
	A.second[2] = square_mean_;
	A.second[2] = square_mean_sd_;
}

void VanillaFloLevel::Reset(){
	level_ = vector<double>();
	mean_ = 0;
	mean_sd_ = 0;
	square_mean_ = 0;
	square_mean_sd_ = 0;
	base_filename_.str("");

	init_ = false;
}

bool VanillaFloLevel::Empty() const{
	return !init_;
}