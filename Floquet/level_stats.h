#ifndef LEVEL_STATS_H
#define LEVEL_STATS_H

#include <sstream>
#include <vector>
#include <utility>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "results.h"
#include "evol_class.h"

using namespace std;
using namespace Eigen;

/**
This file includes output classes which output level statistics from Hamiltonian/unitary
time evolution operator
**/

/*
This class outputs level statistics for vanilla Floquet time evolution operators
*/
class VanillaFloLevel: public ResultsOutput< EvolMatrix< ComplexEigenSolver<MatrixXcd> >*, 
pair<vector<double>, vector<double> > >
{
	private:
		vector<double> level_; // The vector that holds all the level spacings
		
		double mean_; // The mean of level spacing
		double mean_sd_; // Standard deviation of the mean

		double square_mean_; // The square mean of level spacing
		double square_mean_sd_; // The standard deviation of square mean

		stringstream base_filename_; // base filename stringstream

		const bool level_out_; // If it is true, level spacing will be outputted
		const bool mean_out_; // If it is true, the mean of level spacings will be outputted
		const bool square_mean_out_; // If it is true, the mean of square will be outputted
		const bool redirect_; // If true, redirect the output of mean and square mean

		// Process data from single realization. Return a pair of doubles. The first is the
		// mean, and the second is the square mean
		pair<double, double> Single_Data_Process_
		(const EvolMatrix< ComplexEigenSolver<MatrixXcd> >* const);

		bool init_; // Check whether the object has processed one set of data

	public:
		VanillaFloLevel(bool level_out, bool mean_out, bool square_mean_out):
			level_out_(level_out), mean_out_(mean_out), square_mean_out_(square_mean_out),
			redirect_(false), mean_(0), mean_sd_(0), square_mean_(0), square_mean_sd_(0),
			init_(false) {};

		VanillaFloLevel():
			level_out_(false), mean_out_(false), square_mean_out_(false), redirect_(true), 
			mean_(0), mean_sd_(0), square_mean_(0), square_mean_sd_(0), init_(false) {};
			
		void Data_Process(const vector< EvolMatrix< ComplexEigenSolver<MatrixXcd> >* >& );
		void Data_Output(bool output) const;

		/*
		Redirect data to a pair. The first vector takes the level spacings. The second takes mean,
		square mean and their standard deviations. The order is: mean, mean sd, square mean, square
		mean sd.
		*/
		void Data_Redirect(pair< vector<double>, vector<double> >&) const;

		bool Empty() const;

		void Reset();

		~VanillaFloLevel() {};
};

#endif