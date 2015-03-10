#include <cmath>
#include <Eigen/Eigenvalues>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <Python/Python.h>
#include "evol_data_total.h"
#include "methods.h"
#include "generic_func.h"
#include "eigen_output.h"

using namespace std;
using namespace Eigen;

/**
 ** This file contains functions related to full_leftmost_spin_z_per_model. It is the average of
 ** leftmost spin z at a spin chain model for all eigenstates, whose each site has only
 ** spin down and spin up. Spin up is labeled as 1 and spin down is labeled as -1.
 **/

/*
 * Initialize the full_leftmost_spin_z_per_model. The outer index is index for different eigenstate,
 * where the spin z value should be in decreasing order; the inner index is time
 */
void EvolDataTotal::Full_Leftmost_Spin_Z_Per_Model_Init_(const AllPara& parameters){
    const int num_eigenstates = 1 << parameters.generic.size;
    const int time_step = parameters.evolution.time_step;

    full_leftmost_spin_z_per_model_.resize(num_eigenstates);

    for (int i=0; i<num_eigenstates; i++){
        full_leftmost_spin_z_per_model_[i].resize(time_step);

        for (int j=0; j<time_step; j++){
            full_leftmost_spin_z_per_model_[i][j] = 0;
        }
    }
}

/*
 * Compute the leftmost_spin_z at a given time step and initial state, given a complex density matrix
 * in basic binary basis, where 0 is the rightmost position
 */
void EvolDataTotal::Full_Leftmost_Spin_Z_Per_Model_Cal_C_(const MatrixXcd& density_matrix,
        const StepInfo& info){
    const int eigenstate = info.init_condition;
    const int time = info.time;

    cout << "eigenstate num: " << eigenstate << endl;

    // Check whether density_matrix is Hermitian
    if (density_matrix.rows() != density_matrix.cols()){
        cout << "Density matrix passed in entropy_per_model_cal_C is not square." << endl;
        cout << "Rows: " << density_matrix.rows() << endl;
        cout << "Cols: " << density_matrix.cols() << endl;
        abort();
    }

    const int row_num = density_matrix.rows();
    const int size = __builtin_popcount(row_num - 1); // Size of the chain
    const int down = 1 << (size-1); // Number below this will have leftmost spin down

    full_leftmost_spin_z_per_model_[eigenstate][time] = 0;

    for (int i=0; i<row_num; i++){
        if (i<down) full_leftmost_spin_z_per_model_[eigenstate][time] -= real(density_matrix(i,i));
        else full_leftmost_spin_z_per_model_[eigenstate][time] += real(density_matrix(i,i));
    }

    if (info.debug){
        cout << "Average leftmost spin z per model:" << endl;
        cout << full_leftmost_spin_z_per_model_[eigenstate][time] << endl;
        cout << endl;
    }
}

void EvolDataTotal::Full_Leftmost_Spin_Z_Per_Model_Out_(const AllPara& parameters, const string& name){
    const int time_step = parameters.evolution.time_step;
    const int num_eigenstates = 1 << parameters.generic.size;
    const int jump = parameters.evolution.jump;
    const bool output = parameters.output.filename_output;
    const int width = parameters.output.width;
    const double step_size = parameters.evolution.step_size;
    const bool log_time = parameters.evolution.log_time;
    const int log_time_jump = parameters.evolution.log_time_jump;
    const bool markov_jump = parameters.evolution.markov_jump;
    const int markov_time_jump = parameters.evolution.markov_time_jump;
    const bool easy_full_leftmost_spin_z = parameters.multi_ini_para.easy_full_leftmost_spin_z;
    const double non_zero_threshold = parameters.multi_ini_para.non_zero_threshold;
    const bool debug = parameters.generic.debug;

    if (full_leftmost_spin_z_per_model_.size() != num_eigenstates){
        cout << "Not enough eigenstate are computed for leftmost_spin_z." << endl;
        cout << "Expected: " << num_eigenstates << endl;
        cout << "Computed: " << full_leftmost_spin_z_per_model_.size() << endl;
    }

    for(int i=0; i<num_eigenstates; i++){
        if (time_step != full_leftmost_spin_z_per_model_[i].size()){
            cout << "leftmost_spin_z is not computed with enough time steps for eigenstate " << i
                    << endl;
            cout << "Expect: " << time_step << endl;
            cout << "Computed: " << full_leftmost_spin_z_per_model_[i].size() << endl;
            abort();
        }
    }

    stringstream filename, filename_easy;

    filename << name << ",Total=" << time_step;

    if (log_time) filename <<",log_time";

    if (markov_jump) filename <<",markov_jump=" << markov_time_jump;

    filename << ",jump=" << jump;

    if (easy_full_leftmost_spin_z){
        filename_easy << filename.str();
        filename_easy << ",easy_full_leftmost_spin_z,threshold=" << non_zero_threshold
                      << ".txt";
    }

    filename << ",full_left_spin_z_per_model.txt";

    if (output) {
        cout << filename.str() <<endl;

        if (easy_full_leftmost_spin_z) cout << filename_easy.str() << endl;
    }

    ofstream fout( filename.str().c_str() );

    // Compute whether each state has spin z value changed sign during evolution
    // If it is 0, then the initial value is too close to 0. If it is -1, it changes from
    // positive to negative; 1 for changing from negative to positive
    vector<int> spin_z_sign(num_eigenstates);
    if (easy_full_leftmost_spin_z){
        for (int i=0; i < num_eigenstates; i++){
            if (full_leftmost_spin_z_per_model_[i][0] > non_zero_threshold){
                spin_z_sign[i] = 1;
            }
            else if (full_leftmost_spin_z_per_model_[i][0] < -non_zero_threshold){
                spin_z_sign[i] = -1;
            }
            else spin_z_sign[i] = 0;

            if (debug){
                cout << "Initial sign of leftmost spin z:" << endl;
                cout << "For eigenstate " << i <<": ";
                cout << "Value: " << full_leftmost_spin_z_per_model_[i][0] << " sign: "
                     << spin_z_sign[i] << endl;
            }

            if (spin_z_sign[i] != 0){
                for (int t=1; t<time_step; t++){
                    if (spin_z_sign[i] * full_leftmost_spin_z_per_model_[i][t]
                            < - non_zero_threshold){
                        spin_z_sign[i] *= -1;
                        break;
                    }
                }
            }
        }
    }

    if (easy_full_leftmost_spin_z){
        ofstream fout_easy( filename_easy.str().c_str() );

        for (int i=0; i<num_eigenstates;i++) fout_easy << spin_z_sign[i] << endl;
    }

    for (int t=0; t<time_step; t++){
        double time = t * step_size * jump;

        if (markov_jump) time *= markov_time_jump;

        if (log_time){
            long long int power = pow(log_time_jump,t);
            time = step_size * power;
        }

        fout << setw(width) << time;

        if (easy_full_leftmost_spin_z){
            for (int i=0; i<parameters.multi_ini_para.leftmost_spin_z_index_set.size(); i++){
                int state = parameters.multi_ini_para.leftmost_spin_z_index_set[i];

                fout << setw(width) << full_leftmost_spin_z_per_model_[state][t];
            }
        }
        else{
            for (int i=0; i<num_eigenstates; i++)
                fout << setw(width) << full_leftmost_spin_z_per_model_[i][t];
        }

        fout << endl;
    }
}

/*
 * Compute the leftmost_spin_z at a given time step and realization, given a state vector. 
 * Not implemented yet.
 */
void EvolDataTotal::Full_Leftmost_Spin_Z_Per_Model_Cal_(const VectorXcd& init_state,
        const StepInfo& info){
    cout << "Full_Leftmost_Spin_Z_Per_Model_Cal_ Not Implemented Yet." << endl;
    abort();
}