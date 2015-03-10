#ifndef EVOL_DATA_TOTAL_H
#define EVOL_DATA_TOTAL_H

#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <Eigen/Dense>
#include "parameters.h"
#include "evol_data.h"

using namespace std;
using namespace Eigen;

/**
 ** This file includes the class for collecting and outputing data in evolution which include results
 ** from multiple models.
 **/

class EvolDataTotal;

// Functions that initialize data using parameters passed in
typedef void (EvolDataTotal::*Data_Total_Init)(const AllPara&);

// Functions that compute data from a state vector and step information
typedef void (EvolDataTotal::*Data_Total_Cal)(const VectorXcd&, const StepInfo&);

// Functions that compute data from a complex density matrix and step information
typedef void (EvolDataTotal::*Data_Total_Cal_C)(const MatrixXcd&, const StepInfo&);

// Functions that output data for all models using parameters. The string is used as part of
// the output file name
typedef void (EvolDataTotal::*Data_Total_Out)(const AllPara&, const string&);

class EvolDataTotal
{
private:
    map<string,bool> func_status_; // Determine whether a particular data type is calculated


    map<string, Data_Total_Init> data_init_; // Correlate name with data_init function
    map<string, Data_Total_Cal> data_cal_; // Correlate name with data_cal function
    map<string, Data_Total_Cal_C> data_cal_C_; // Correlate name with data_cal_C function
    map<string, Data_Total_Out> data_out_total_; // Correlate name with data_out_total function

    void Data_Func_Map_Init_(); // Initialize data_init_ and data_cal_;
    void Name_Check_() const; // Check names in different maps are consistent

    // Average leftmost spin z for multiple models, each with one run. The outer index is
    // for time; the inner index is for model
    vector<vector<double> > leftmost_spin_z_one_run_;
    // Initialize leftmost_spin_z_one_run
    void Leftmost_Spin_Z_One_Run_Init_(const AllPara&);
    // Compute leftmost_spin_z_one_run given state vector, not implemented yet
    void Leftmost_Spin_Z_One_Run_Cal_(const VectorXcd&, const StepInfo&);
    // Compute leftmost_spin_z_one_run given complex density matrix
    void Leftmost_Spin_Z_One_Run_Cal_C_(const MatrixXcd&, const StepInfo&);
    // Output leftmost_spin_z_one_run
    void Leftmost_Spin_Z_One_Run_Out_(const AllPara&, const string&);

    // Average leftmost spin z for multiple initial settings, each with one run. The outer index
    // is for time; the inner index is for initial setting. The values from all initial conditions
    // are put together
    vector<vector<double> > full_leftmost_spin_z_per_model_;
    // Initialize full_leftmost_spin_z_per_model
    void Full_Leftmost_Spin_Z_Per_Model_Init_(const AllPara&);
    // Compute leftmost_spin_z_per_model given state vector, not implemented yet
    void Full_Leftmost_Spin_Z_Per_Model_Cal_(const VectorXcd&, const StepInfo&);
    // Compute leftmost_spin_z_per_model given complex density matrix
    void Full_Leftmost_Spin_Z_Per_Model_Cal_C_(const MatrixXcd&, const StepInfo&);
    // Output leftmost_spin_z_per_model
    void Full_Leftmost_Spin_Z_Per_Model_Out_(const AllPara&, const string&);

    const int size_; // Size of the system

public:
    EvolDataTotal(const AllPara&);

    void Print_All_Name() const; // Print all possible names of data type calculation
    void Print_All_Status() const; // Print all possible data type calculation, indicating
    // whether they will be calculated

    // Compute data at each step given a state vector. The entry which is true in
    // func_status will be computed.
    void Data_Compute(const VectorXcd&, const StepInfo&);

    // Compute data at each step given a complex density matrix. The entry which is true in
    // func_status will be computed.
    void Data_Compute(const MatrixXcd&, const StepInfo&);


    // Output data of all models to file. All the data that are computed will be outputted
    void Data_Total_Output(const AllPara&, const string&);
};

#endif