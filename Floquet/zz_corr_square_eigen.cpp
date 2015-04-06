//
// Created by Liangsheng Zhang on 4/6/15.
//

#include <vector>
#include <complex>
#include <iomanip>
#include "flo_eigen_func.h"
#include "methods.h"

using namespace std;

/**
 ** This file contains functions related to eigenstates and zz_corr, which is the end-to-end zz correlation square
 ** used in studying floquet model transition from thermalization to localization
 **/

/*
 * Initialize zz_corr_square in eigen_data_. The outer index is for realization, and the inner index is for
 * different eigenstates under the same realization. Here assume the model is binary
 */
void FloEigenFunc::ZZ_corr_square_eigen_init_(AllPara const & parameters) {
    const int num_realization = parameters.generic.num_realizations;
    const int dim = 1 << parameters.generic.size;

    eigen_data_.zz_corr_square.resize(num_realization);

    for (int i=0; i<num_realization; i++){
        eigen_data_.zz_corr_square[i].resize(dim);

        for (int j=0; j<dim;j++) eigen_data_.zz_corr_square[i][j] = 0;
    }
}

/*
 * Compute zz_corr_square given realization index
 */
void FloEigenFunc::ZZ_corr_square_eigen_compute_(AllPara const & parameters,
                    const EvolMatrix<ComplexEigenSolver<MatrixXcd> >* floquet, const LocalInfo& local_info) {

    const bool debug = parameters.generic.debug;

    // Vector for eigenvectors in basic binary basis
    vector<vector<complex<double> > > evec_basic;
    const int dim = local_info.dim;
    const int half_dim = dim / 2; // Below which leftmost spin is -1

    evec_basic.resize(dim);
    for (int i=0; i<dim;i++) evec_basic[i].resize(dim);

    // Convert eigenvectors in basic basis
    evec_to_basic(floquet, evec_basic);

    for (int i=0; i<evec_basic.size(); i++){
        // Record values for each eigenstate
        double left_temp = 0;
        double right_temp = 0;
        double left_right_temp = 0;
        double zz_square = 0;

        for (int j=0; j<evec_basic[i].size(); j++){
            int right_spin = -1;
            int left_spin = 1;
            if (j & 1 == 1) right_spin = 1;
            if (j < half_dim ) left_spin = -1;

            left_temp += left_spin * norm(evec_basic[i][j]);
            right_temp += right_spin * norm(evec_basic[i][j]);
            left_right_temp += left_spin * right_spin * norm(evec_basic[i][j]);
        }

        zz_square = (left_right_temp - left_temp*right_temp ) * (left_right_temp - left_temp*right_temp );
        if (debug){
            cout << "Realization " << local_info.realization_index << " eigenstate " << i << ":" << endl;
            cout << "Average left spin: " << left_temp << endl;
            cout << "Average right spin: " << right_temp << endl;
            cout << "Average left X right: " << left_right_temp << endl;
        }

        eigen_data_.zz_corr_square[local_info.realization_index][i] += zz_square;
    }
}

/*
 * Output averages of mean end-to-end zz correlation squares for all eigenstates and realizations
 */
void FloEigenFunc::ZZ_corr_square_eigen_out_(AllPara const & parameters, const string& name) {
    const double J = parameters.floquet.J;
    const int num_realizations = parameters.generic.num_realizations;
    const int version = parameters.generic.version;
    const int width = parameters.output.width;
    const bool output = parameters.output.filename_output;
    const int size = parameters.generic.size;

    if (eigen_data_.zz_corr_square.size() != num_realizations){
        cout << "Not enough number of realizations for zz correlation square." << endl;
        cout << "Expected Number: " << num_realizations << endl;
        cout << "Actual number: " << eigen_data_.zz_corr_square.size() << endl;
        abort();
    }

    stringstream filename;
    filename << name << ",size=" << size << ",Run=" << num_realizations << ",J=" << J << ",all_eigen_zz_corre_square";
    if (version > 0) filename <<",v" << version;
    filename << ".txt";

    if (output) cout << filename.str() <<endl;

    ofstream fout( filename.str().c_str() );

    for (int i=0; i<num_realizations;i++){
        for (int j=0; j<eigen_data_.zz_corr_square[i].size();j++)
            fout << setw(width) << eigen_data_.zz_corr_square[i][j];
        fout << endl;
    }
}


