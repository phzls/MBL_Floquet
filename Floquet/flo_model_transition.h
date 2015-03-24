//
// Created by Liangsheng Zhang on 3/24/15.
//

#ifndef _FLO_MODEL_TRANSITION_H_
#define _FLO_MODEL_TRANSITION_H_

#include <iostream>
#include <vector>
#include "parameters.h"
#include "evol_class.h"

using namespace std;

/*
 * The data class used in flo_model_transition. It contains all possible outputs.
 */

struct ModelData
{
    // End-to-end sigma_z-sigma_z correlation square averaged over all eigenstates and realizations
    vector< vector<double> > zz_corr_square;
};

/*
 * Local information that helps computation
 */
struct LocalInfo
{
    int J_index; // The index of J
    int realization_index; // The index of realization
};

class FloModelTransition;

// Pointer to functions initializing data
typedef  void (FloModelTransition::*Flo_init)(const AllPara&);

// Pointer to functions studying propertieso of a given floquet system
typedef void (FloModelTransition::*Flo_func)(const AllPara&, const EvolMatrix<ComplexEigenSolver<MatrixXcd> >*,
        LocalInfo&);

// Pointer to functions outputting data. The string gives the first part of filename
typedef void (FloModelTransition::*Flo_out)(const AllPara&, const string&);

/*
 * The class that contains data and functions which handle it
 */
class FloModelTransition
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

    ModelData model_data_;

    // For end to end zz correlation square
    void ZZ_corr_square_init_(const AllPara&);
    void ZZ_corr_square_compute_(const AllPara&, const EvolMatrix<ComplexEigenSolver<MatrixXcd> >*, LocalInfo&);
    void ZZ_corr_square_out_(const AllPara&, const string&);

public:
    FloModelTransition(const AllPara& parameters) : flo_func_bool_map_(parameters.transition.flo_transition_compute) {
        map_initialize_(const AllPara&);
    }

    void Compute(const AllPara&, const EvolMatrix<ComplexEigenSolver<MatrixXcd> >*, const LocalInfo&);
    void Output(const AllPara&, const string&);

    ~FloModelTransition() {};
};

#endif //_FLO_MODEL_TRANSITION_H_
