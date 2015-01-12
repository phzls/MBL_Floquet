#include <fstream>
#include <iostream>
#include <complex>
#include <cmath>
#include <iomanip>
#include "level_stats.h"
#include "sort.h"
#include "constants.h"
#include "output_func.h"

using namespace std;

pair<double, double> VanillaFloLevel::Single_Data_Process_(const 
EvolMatrix< ComplexEigenSolver<MatrixXcd> >* const U, const int position, const int dim){

	// The total dimension of Hilbert space
	const int local_dim = U -> Get_Dim();
	if (local_dim != dim){
		cout << "The dimensions are not consistent." << endl;
		abort();
	}

	// The average level spacing of phases of eigenvalues
	const double ave_spacing = 2*Pi / dim;

	// The dimension computed from eigenvalues
	int eigen_dim = 0;

	for (int i=0; i< U -> eigen.size();i++){
		eigen_dim += U -> eigen[i].eigenvalues().rows();
	}

	if (dim != eigen_dim){
		cout << "The dimension of eigenvalues of time evolution operator is not correct."<<endl;
		abort();
	}

	vector<double> phases(dim); // Phases of all eigenvalues in ascending order

	int index = 0;
	for (int i=0; i< U -> eigen.size(); i++){
		for (int j=0; j< U -> eigen[i].eigenvalues.rows(); j++){
			phase[index] = arg( U -> eigen[i].eigenvalues()(j) );
			index ++; 
		}
	}

	// Sort the phases in ascending order
	sort(phases.begin(), phases.end(), Vec_Double_Sort);

	// Compute level statistics
	double local_mean = 0;
	double local_square_mean = 0;

	double spacing; // level spacing

	for (int i=0; i<dim-1; i++){
		spacing = ( phases[i+1] - phases[i] ) / ave_spacing;
		level_[position + i] = spacing;
		local_mean += spacing;
		local_square_mean += spacing * spacing;
	}
	
	spacing = ( phases[0] + 2*Pi - phases[dim-1] ) / ave_spacing; 
	level_[position + dim - 1] = spacing;
	local_mean += spacing;
	local_square_mean += spacing * spacing;

	pair<double, double> mean_pair(local_mean/dim, local_square_mean/dim);
	return mean_pair;
}

void VanillaFloLevel::Data_Process(const vector< EvolMatrix< ComplexEigenSolver<MatrixXcd> >* >& U)
{
	const int num_realization = U.size();
	const int dim = U[0] -> Get_Dim();

	if (init_){
		cout<< "The object is already holding one set of data." <<endl;
		abort();
	}

	level_.resize(dim * num_realization);

	for (int i=0; i<num_realization; i++){

		pair<double, double> mean_square_mean;
		mean_square_mean = Single_Data_Process_(U[i], i*dim, dim);

		mean_ += mean_square_mean.first;
		mean_sd_ += mean_square_mean.first * mean_square_mean.first;

		square_mean_ += mean_square_mean.second;
		square_mean_sd_ += mean_square_mean.second * mean_square_mean.second;
	}

	// if level_ has not been initialized
	if (num_realization != 0 && level_.size() == 0){
		cout <<"Level spacing calculation fails." << endl;
		abort();
	}

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

	// The object now holds the data
	init_ = true;
}

void VanillaFloLevel::Data_Output(bool output, int width) const{
	if (!init_){
		cout<< "The object has no data." <<endl;
		abort();
	}

	if (!redirect_ && level_out_){

		string post_string = ",level_spacing.txt";
		ofstream level_output; 
		Of_Construct(level_output, base_filename_, post_string, output) ;

		for(int i=0; i<level_.size();i++){
			level_output << level_[i]<<endl;
		} 
	}

	if (!redirect_ && mean_out_){

		string post_string = ",level_spacing_mean.txt";
		ofstream mean_output; 
		Of_Construct(mean_output, base_filename_, post_string, output) ;

		vector<double> temp(2);
		temp[0] = mean_;
		temp[1] = mean_sd_;

		Write_File(mean_output, temp, width);
	}

	if (!redirect_ && square_mean_out_){

		string post_string = ",level_spacing_square_mean.txt";
		ofstream square_mean_output;
		Of_Construct(square_mean_output, base_filename_, post_string, output) ;

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