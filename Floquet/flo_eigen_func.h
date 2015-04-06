//
// Created by Liangsheng Zhang on 4/3/15.
//

#ifndef FLOQUET_FLO_EIGEN_FUNC_H
#define FLOQUET_FLO_EIGEN_FUNC_H

#include <iostream>
#include <vector>
#include <map>
#include "parameters.h"
#include "evol_class.h"

using namespace std;

/**
 ** This file contains relevant classes for flo_eigen
 **/

/*
 * The data class used in flo_eigen. It contains all possible outputs.
 */

struct EigenData
{
    // End-to-end sigma_z-sigma_z correlation square for all eigenstates in all realizations
    vector< vector<double> > zz_corr_square;

    // Entropy varaince for all eigenstates
    vector< vector<double> > ent_var;

    // Entropy variance for the eigenstate with smallest phase magnitude among different realizations
    vector< vector<double> > ent_smallest_var;

    // End-to-end sigma_z-sigma_z time four-point correlation
    vector< vector<double> > zz_time_corr;

    // Various components of zz time correlation
    vector< vector<vector<double> > > zz_time_corr_component;

};

/*
 * Local information that helps computation
 */
struct LocalInfo
{
    int realization_index; // The index of realization
    int dim; // The dimension of the system
};

class FloEigenFunc;

// Pointer to functions initializing data
typedef  void (FloEigenFunc::*Flo_init)(const AllPara&);

// Pointer to functions studying properties of a given floquet system
typedef void (FloEigenFunc::*Flo_func)(const AllPara&, const EvolMatrix<ComplexEigenSolver<MatrixXcd> >*,
                                             const LocalInfo&);

// Pointer to functions outputting data. The string gives the first part of filename
typedef void (FloEigenFunc::*Flo_out)(const AllPara&, const string&);

/*
 * The class that contains data and functions which handle it
 */
class FloEigenFunc
{
private:
    // Map that determines which functions are computed
    map<string, bool> flo_func_bool_map_;

    // Map that links strings to initialization functions
    map<string, Flo_init> flo_init_map_;

    // Map that links strings to computation functions
    map<string, Flo_func> flo_func_map_;

    // Map that links strings to output functions
    map<string, Flo_out> flo_out_map_;

    // Initialize all maps
    void map_initialize_(const AllPara&);

    EigenData eigen_data_;

    // For end to end zz correlation square statistics
    void ZZ_corr_square_eigen_init_(const AllPara&);
    void ZZ_corr_square_eigen_compute_(const AllPara&, const EvolMatrix<ComplexEigenSolver<MatrixXcd> >*,
                                       const LocalInfo&);
    void ZZ_corr_square_eigen_out_(const AllPara&, const string&);

public:
    FloEigenFunc(const AllPara& parameters) : flo_func_bool_map_(parameters.transition.flo_transition_compute) {
        map_initialize_(parameters);
    }

    void Compute(const AllPara&, const EvolMatrix<ComplexEigenSolver<MatrixXcd> >*, const LocalInfo&);
    void Output(const AllPara&, const string&);

    ~FloEigenFunc() {};
};


#endif //FLOQUET_FLO_EIGEN_FUNC_H
