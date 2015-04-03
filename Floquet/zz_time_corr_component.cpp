//
// Created by Liangsheng Zhang on 4/3/15.
//

#include <vector>
#include <complex>
#include <iomanip>
#include "flo_model_transition.h"
#include "methods.h"
#include "generic_func.h"

using namespace std;

/**
 ** This file contains functions related to zz_time_corr_component, which reports individual components of the
 ** end-to-end zz time four point correlation used in studying floquet model transition from thermalization to
 ** localization. The initial system is given a density matrix of identity matrix
 **/

/*
 * Initialize zz_time_corr_component in model_data_. The outer index is for J, and the middle index is for different
 * components, and the inner index is for realization
 */
void FloModelTransition::ZZ_time_corr_component_init_(AllPara const & parameters) {
    const int J_N = parameters.floquet.J_N; // Number of points for J
    const int num_realization = parameters.generic.num_realizations;

    model_data_.zz_time_corr_component.resize(J_N);

    for (int i=0; i<J_N; i++){
        model_data_.zz_time_corr_component[i].resize(3); // Three components

        for (int j=0; j<3;j++) model_data_.zz_time_corr_component[i][j].resize(num_realization);
    }
}

/*
 * Compute zz_time_corr_component given J index and realization index
 */
void FloModelTransition::ZZ_time_corr_component_compute_(AllPara const & parameters,
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
    double zz = 0;
    double z_left = 0;
    double z_right = 0;
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

        zz += left_right_temp * left_right_temp;
        z_left += left_temp * left_temp;
        z_right += right_temp * right_temp;
        z_left_ave_right_ave += left_temp * right_temp;

        if (debug){
            cout << "Realization " << local_info.realization_index << " eigenstate " << i << ":" << endl;
            cout << "Average left spin: " << left_temp << endl;
            cout << "Average right spin: " << right_temp << endl;
            cout << "Average left X right: " << left_right_temp << endl;
        }
    }

    zz /= double(dim);
    z_left /= double(dim);
    z_right /= double(dim);
    z_left_ave_right_ave /= double(dim);


    model_data_.zz_time_corr_component[local_info.J_index][0][local_info.realization_index] += zz;
    model_data_.zz_time_corr_component[local_info.J_index][1][local_info.realization_index] += z_left * z_right;
    model_data_.zz_time_corr_component[local_info.J_index][2][local_info.realization_index] +=
            z_left_ave_right_ave * z_left_ave_right_ave;
}

/*
 * Output averages of mean end-to-end zz time four-point correlation components and their standard deviations among all
 * realizations for each J
 */
void FloModelTransition::ZZ_time_corr_component_out_(AllPara const & parameters, const string& name) {
    const double J_min = parameters.floquet.J_min;
    const double J_max = parameters.floquet.J_max;
    const int J_N = parameters.floquet.J_N;
    const int num_realizations = parameters.generic.num_realizations;
    const int version = parameters.generic.version;
    const int width = parameters.output.width;
    const bool output = parameters.output.filename_output;
    const int size = parameters.generic.size;

    if (model_data_.zz_time_corr_component.size() != J_N){
        cout << "Not enough number of J for zz time correlation components." << endl;
        cout << "Expected Number: " << J_N << endl;
        cout << "Actual number: " << model_data_.zz_time_corr_component.size() << endl;
        abort();
    }

    for (int i=0; i< J_N; i++){
        for (int j=0; j<3; j++){
            if (model_data_.zz_time_corr_component[i][j].size() != num_realizations){
                cout << "Not enough number of realizations at " << i <<"th J for zz time correlation" <<
                        "component " << j << endl;
                cout << "Expected Number: " << num_realizations << endl;
                cout << "Actual Number: " << model_data_.zz_time_corr_component[i][j].size() << endl;
                abort();
            }
        }
    }

    stringstream filename;
    filename << name << ",size=" << size << ",Run=" << num_realizations << ",J_N=" << J_N << ",J_min=" << J_min
    << ",J_max=" << J_max << ",zz_time_corr_component";
    if (version > 0) filename <<",v" << version;
    filename << ".txt";

    if (output) cout << filename.str() <<endl;

    ofstream fout( filename.str().c_str() );

    vector<string> component_name(3);
    component_name[0] = "<zz>^2";
    component_name[1] = "<z>^2<z>^2";
    component_name[2] = "(<z><z>)^2";

    fout << setw(10) << "J";
    for (int i=0; i<3; i++){
        fout << setw(width) << component_name[i] << " mean" << setw(width) << component_name[i] << " sd";
    }
    fout << endl;

    for (int i=0; i<J_N;i++){
        double J;
        double mean,sd;
        if (J_N > 1) J = J_min + i * (J_max - J_min)/(J_N-1);
        else J = J_min;

        fout << setw(10) << J;

        for (int j=0; j<3; j++){
            generic_mean_sd(model_data_.zz_time_corr_component[i][j], mean, sd);
            fout << setw(width) << mean << setw(width) << sd;
        }
        fout << endl;
    }
}



