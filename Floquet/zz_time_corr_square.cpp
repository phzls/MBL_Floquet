//
// Created by Liangsheng Zhang on 3/30/15.
//

#include <vector>
#include <complex>
#include <iomanip>
#include "flo_model_transition.h"
#include "methods.h"
#include "generic_func.h"

using namespace std;

/**
 ** This file contains functions related to zz_time_corr, which is the end-to-end zz time four point correlation square
 ** used in studying floquet model transition from thermalization to localization. The initial system is given a
 ** density matrix of identity matrix
 **/

/*
 * Initialize zz_time_corr_square in model_data_. The outer index is for J, and the inner index is for realization
 */
void FloModelTransition::ZZ_time_corr_square_init_(AllPara const & parameters) {
    const int J_N = parameters.floquet.J_N; // Number of points for J
    const int num_realization = parameters.generic.num_realizations;

    model_data_.zz_time_corr_square.resize(J_N);

    for (int i=0; i<J_N; i++){
        model_data_.zz_time_corr_square[i].resize(num_realization);

        for (int j=0; j<num_realization;j++) model_data_.zz_time_corr_square[i][j] = 0;
    }
}

/*
 * Compute zz_time_corr_square given J index and realization index
 */
void FloModelTransition::ZZ_time_corr_square_compute_(AllPara const & parameters,
                                                      const EvolMatrix<ComplexEigenSolver<MatrixXcd> >* floquet,
                                                      const LocalInfo& local_info) {

    const bool debug = parameters.generic.debug;

    // Vector for eigenvectors in basic binary basis
    vector<vector<complex<double> > > evec_basic;
    const int dim = floquet -> Get_Dim();
    const int half_dim = dim / 2; // Below which leftmost spin is -1

    evec_basic.resize(dim);
    for (int i=0; i<dim;i++) evec_basic[i].resize(dim);

    // Convert eigenvectors in basic basis
    evec_to_basic(floquet, evec_basic);

    // Record values for all eigenstate
    double zz_square = 0;
    double z_left_square = 0;
    double z_right_square = 0;
    double z_left_ave_right_ave = 0;

    for (int i=0; i<evec_basic.size(); i++){
        // Record values for each eigenstate
        double left_temp = 0;
        double right_temp = 0;
        double left_right_temp = 0;

        for (int j=0; j<evec_basic[i].size(); j++){
            int right_spin = -1;
            int left_spin = 1;
            if (j & 1 == 1) right_spin = 1;
            if (j < half_dim ) left_spin = -1;

            left_temp += left_spin * norm(evec_basic[i][j]);
            right_temp += right_spin * norm(evec_basic[i][j]);
            left_right_temp += left_spin * right_spin * norm(evec_basic[i][j]);
        }

        zz_square += left_right_temp * left_right_temp;
        z_left_square += left_temp * left_temp;
        z_right_square += right_temp * right_temp;
        z_left_ave_right_ave += left_temp * right_temp;

        if (debug){
            cout << "Realization " << local_info.realization_index << " eigenstate " << i << ":" << endl;
            cout << "Average left spin: " << left_temp << endl;
            cout << "Average right spin: " << right_temp << endl;
            cout << "Average left X right: " << left_right_temp << endl;
        }
    }

    zz_square /= double(dim);
    z_left_square /= double(dim);
    z_right_square /= double(dim);
    z_left_ave_right_ave /= double(dim);

    double zz_time_corr = zz_square - z_left_square * z_right_square - z_left_ave_right_ave * z_left_ave_right_ave;

    model_data_.zz_time_corr_square[local_info.J_index][local_info.realization_index] += zz_time_corr * zz_time_corr;
    if (debug){
        cout << "Realization " << local_info.realization_index << " average time correlation: "
        << zz_time_corr << endl;
        cout << "Realization " << local_info.realization_index << " average time correlation square: "
        << zz_time_corr * zz_time_corr << endl;
        cout << endl;
    }
}

/*
 * Output averages of mean end-to-end zz time four-point correlation squares and their standard deviations among all
 * realizations for each J
 */
void FloModelTransition::ZZ_time_corr_square_out_(AllPara const & parameters, const string& name) {
    const double J_min = parameters.floquet.J_min;
    const double J_max = parameters.floquet.J_max;
    const int J_N = parameters.floquet.J_N;
    const int num_realizations = parameters.generic.num_realizations;
    const int version = parameters.generic.version;
    const int width = parameters.output.width;
    const bool output = parameters.output.filename_output;
    const int size = parameters.generic.size;

    if (model_data_.zz_time_corr_square.size() != J_N){
        cout << "Not enough number of J for zz correlation square." << endl;
        cout << "Expected Number: " << J_N << endl;
        cout << "Actual number: " << model_data_.zz_time_corr_square.size() << endl;
        abort();
    }

    for (int i=0; i< J_N; i++){
        if (model_data_.zz_time_corr_square[i].size() != num_realizations){
            cout << "Not enough number of realizations at " << i <<"th J for zz correlation square." << endl;
            cout << "Expected Number: " << num_realizations << endl;
            cout << "Actual Number: " << model_data_.zz_time_corr_square[i].size() << endl;
            abort();
        }
    }

    stringstream filename;
    filename << name << ",size=" << size << ",Run=" << num_realizations << ",J_N=" << J_N << ",J_min=" << J_min
    << ",J_max=" << J_max << ",zz_time_corr_square";
    if (version > 0) filename <<",v" << version;
    filename << ".txt";

    if (output) cout << filename.str() <<endl;

    ofstream fout( filename.str().c_str() );

    for (int i=0; i<J_N;i++){
        double J;
        double mean,sd;
        if (J_N > 1) J = J_min + i * (J_max - J_min)/(J_N-1);
        else J = J_min;

        generic_mean_sd(model_data_.zz_time_corr_square[i], mean, sd);
        fout << setw(10) << J << setw(width) << mean << setw(width) << sd << endl;
    }
}

