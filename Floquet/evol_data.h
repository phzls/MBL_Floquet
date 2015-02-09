#ifndef EVOL_DATA_H
#define EVOL_DATA_H

#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <Eigen/Dense>
#include "parameters.h"

using namespace std;
using namespace Eigen;

/**
 ** This file includes the class for collecting and outputing data in evolution
 **/

/*
 * This structure gives relevant information at each step that can be used for computation.
 */
struct StepInfo{
	int model; // The index of current model
	int realization; // The index of current realization
	int time; // The index of current time step
	bool debug; // Whether output debug information

	double delta; // A small number

	int left_size; // If partition the chain to two halves, the size of left part
};

class EvolData;

// Functions that initialize data using parameters passed in
typedef void (EvolData::*Data_Init)(const AllPara&);

// Functions that compute data from a state vector and step information
typedef void (EvolData::*Data_Cal)(const VectorXcd&, const StepInfo&);

// Functions that compute data from a complex density matrix and step information
typedef void (EvolData::*Data_Cal_C)(const MatrixXcd&, const StepInfo&);

// Functions that output data using parameters. The string is used as part of the output
// file name
typedef void (EvolData::*Data_Out)(const AllPara&, const string&);

class EvolData
{	
	private:
		map<string,bool> func_status_; // Determine whether a particular data type is calculated
		
		
		map<string, Data_Init> data_init_; // Correlate name with data_init function
		map<string, Data_Cal> data_cal_; // Correlate name with data_cal function
		map<string, Data_Cal_C> data_cal_C_; // Correlate name with data_cal_C function
		map<string, Data_Out> data_out_; // Correlate name with data_out function

		void Data_Func_Map_Init_(); // Initialize data_init_ and data_cal_;
		void Name_Check_() const; // Check names in different maps are consistent

		// Entropy per model. The outer index is for time; the inner index is for realization
		vector<vector<double> > entropy_per_model_;
		void Entropy_Per_Model_Init_(const AllPara&); // Initialize entropy_per_model
		// Compute entropy_per_model given state vector
		void Entropy_Per_Model_Cal_(const VectorXcd&, const StepInfo&); 
		// Compute entropy_per_model given complex density matrix
		void Entropy_Per_Model_Cal_C_(const MatrixXcd&, const StepInfo&); 
		// Output entropy_per_model
		void Entropy_Per_Model_Out_(const AllPara&, const string&);

		// Average leftmost spin z per model. The outer index is for time; the inner index is for
		// realization
		vector<vector<double> > leftmost_spin_z_per_model_;
		// Initialize leftmost_spin_z_per_model
		void Leftmost_Spin_Z_Per_Model_Init_(const AllPara&);
		// Compute leftmost_spin_z_per_model given state vector, not implemented yet
		void Leftmost_Spin_Z_Per_Model_Cal_(const VectorXcd&, const StepInfo&); 
		// Compute leftmost_spin_z_per_model given complex density matrix
		void Leftmost_Spin_Z_Per_Model_Cal_C_(const MatrixXcd&, const StepInfo&);
		// Output leftmost_spin_z_per_model
		void Leftmost_Spin_Z_Per_Model_Out_(const AllPara&, const string&);

		const int size_; // Size of the system
		
	public:
		EvolData(const AllPara&);

		void Print_All_Name() const; // Print all possible names of data type calculation
		void Print_All_Status() const; // Print all possible data type calculation, indicating
									   // whether they will be calculated

		// Compute data at each step given a state vector. The entry which is true in 
		// func_status will be computed.
		void Data_Compute(const VectorXcd&, const StepInfo&); 

		// Compute data at each step given a complex density matrix. The entry which is true in 
		// func_status will be computed.
		void Data_Compute(const MatrixXcd&, const StepInfo&); 

		// Output data to file. All the data that are computed will be outputted
		void Data_Output(const AllPara&, const string&);
};

#endif