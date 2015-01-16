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
	int model;
	int realization;
	int time_step;

	int left_size; // If partition the chain to two halves, the size of left part
};

class EvolData;

typedef void (EvolData::*Data_Init)(const AllPara&);
typedef void (EvolData::*Data_Cal)(const VectorXcd&, const StepInfo&);

class EvolData
{	
	private:
		map<string,bool> func_status_; // Determine whether a particular data type is calculated
		
		
		map<string, Data_Init> data_init_; // Correlate name with data_init function
		map<string, Data_Cal> data_cal_; // Correlate name with data_cal function

		void Data_Func_Map_Init_(); // Initialize data_init_ and data_cal_;
		void Name_Check_() const; // Check names in different maps are consistent

		void Data_Init_(); // Initialize data according to func_stauts

		// Entropy per model. The outer index is for realization; the inner index is for time
		vector<vector<double> > entropy_per_model_;
		void Entropy_Per_Model_Init_(const AllPara&); // Initialize entropy_per_model
		// Compute entropy_per_model 
		void Entropy_Per_Model_Cal_(const VectorXcd&, const StepInfo&); 

		const int size_; // Size of the system
		
	public:
		EvolData(const AllPara&);

		void Print_All_Name() const; // Print all possible names of data type calculation
		void Print_All_Status() const; // Print all possible data type calculation, indicating
									   // whether they will be calculated
};

#endif