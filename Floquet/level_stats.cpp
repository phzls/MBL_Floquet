#include <fstream>
#include <iostream>
#include <complex>
#include <cmath>
#include <iomanip>
#include "level_stats.h"
#include "sort.h"
#include "parameter.h"
#include "output_func.h"

using namespace std;

pair<double, double> VanillaFloLevel::Single_Data_Process_(const 
EvolMatrix< ComplexEigenSolver<MatrixXcd> >* const U){

	// The total dimension of Hilbert space
	const int dim = U -> GetDim();

	// The average level spacing of phases of eigenvalues
	const double ave_spacing = 2*Pi / dim;

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
		phases[i] = arg(U -> eigen.eigenvalues()(i));
	}

	// Sort the phases in ascending order
	sort(phases.begin(), phases.end(), Vec_Double_Sort);

	// Compute level statistics
	double local_mean = 0;
	double local_square_mean = 0;

	double spacing;

	for (int i=0; i<dim-1; i++){
		spacing = ( phases[i+1] - phases[i] ) / ave_spacing;
		level_[i] += spacing;
		local_mean += spacing;
		local_square_mean += spacing * spacing;
	}
	
	spacing = ( phases[0] + 2*Pi - phases[dim-1] ) / ave_spacing;
	level_[dim-1] += spacing; 
	local_mean += spacing;
	local_square_mean += spacing * spacing;

	pair<double, double> mean_pair(local_mean/dim, local_square_mean/dim);
	return mean_pair;
}

void VanillaFloLevel::Data_Process(const vector< EvolMatrix< ComplexEigenSolver<MatrixXcd> >* >& U)
{
	const int num_realization = U.size();

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
	for (int i=1; i< num_realization; i++){
		if (base_filename_.str() != U[i] -> Repr()){
			cout <<"String Representations of models are not consistent!"<<endl;
			abort();
		}
	}
	base_filename_ <<",Realizations="<<num_realization;

	init_ = true;
}

void VanillaFloLevel::Data_Output(bool output, int width) const{
	if (!init_){
		cout<< "The object has no data." <<endl;
		abort();
	}

	if (!redirect_ && level_out_){

		string post_string = ",level_spacing.txt";
		ofstream level_output = Of_Construct(base_filename, post_string, output) ;

		for(int i=0; i<level_.size();i++){
			level_output << level_[i]<<endl;
		} 
	}

	if (!redirect_ && mean_out_){

		string post_string = ",level_spacing_mean.txt";
		ofstream mean_output = Of_Construct(base_filename, post_string, output) ;

		vector<double> temp(2);
		temp[0] = mean_;
		temp[1] = mean_sd_;

		Write_File(mean_output, temp, width);
	}

	if (!redirect_ && square_mean_out_){

		string post_string = ",level_spacing_square_mean.txt";
		ofstream square_mean_output = Of_Construct(base_filename, post_string, output) ;

		vector<double> temp(2);
		temp[0] = square_mean_;
		temp[1] = square_mean_sd_;

		Write_File(square_mean_output, temp, width);
	}
}

void VanillaFloLevel::Data_Redirect(pair< vector<double>, vector<double> >& A) const{
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
	A.second[3] = square_mean_sd_;
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