//
// Created by Liangsheng Zhang on 3/26/15.
//

#include <vector>
#include <utility>
#include <complex>
#include <iomanip>
#include <cmath>
#include <Eigen/Eigenvalues>
#include "flo_model_transition.h"
#include "methods.h"
#include "generic_func.h"
#include "eigen_output.h"
#include "sort.h"

using namespace std;
using namespace Eigen;

/**
 ** This file contains functions related to ent_smallest_var, which is the entropy variance for
 ** the eigenstate with smallest phase magnitude among all realizations. The entropy is half chain entanglement
 ** entropy, so the spin chain must have even length.
 **/

/*
 * Initialize ent_smallest_var in model_data_. The outer index is for J, and the inner index is for realization
 */
void FloModelTransition::Ent_smallest_var_init_(AllPara const & parameters) {
    const int J_N = parameters.floquet.J_N; // Number of points for J
    const int num_realization = parameters.generic.num_realizations;

    model_data_.ent_smallest_var.resize(J_N);

    for (int i=0; i<J_N; i++){
        model_data_.ent_smallest_var[i].resize(num_realization);

        for (int j=0; j<num_realization;j++) model_data_.ent_smallest_var[i][j] = 0;
    }
}

/*
 * Compute ent_smallest_var given J index and realization index
 */
void FloModelTransition::Ent_smallest_var_compute_(AllPara const & parameters,
        const EvolMatrix<ComplexEigenSolver<MatrixXcd> >* floquet, const LocalInfo& local_info) {

    const bool debug = parameters.generic.debug;

    // Vector for eigenvectors in basic binary basis
    vector<vector<complex<double> > > evec_basic;
    const int dim = floquet -> Get_Dim();
    const int size = floquet -> Get_Size(); // Size of the chain

    if (size %2 != 0){
        cout << "Length of the chain must be even for entropy variance calculation in flo_transition." << endl;
        cout << "Spin chain length: " << size << endl;
        abort();
    }

    evec_basic.resize(dim);
    for (int i=0; i<dim;i++) evec_basic[i].resize(dim);

    // Convert eigenvectors in basic basis
    evec_to_basic(floquet, evec_basic);

    VectorXcd evec(dim); // Vector for one eigenstate in binary basis
    MatrixXcd reduced_d; // Reduced density matrix
    double ent = 0; // Entropy
    vector<pair<double,int> > eigen_mag(dim); // First position is magnitude of energy; second position

    int index = 0;
    for (int i=0; i<floquet -> eigen.size(); i++){
        for (int j=0; j<floquet -> eigen[i].eigenvalues().rows(); j++){
            eigen_mag[index].first = abs( arg(floquet -> eigen[i].eigenvalues()[j]) );
            eigen_mag[index].second = index;
            index ++;
        }
    }

    if (index != dim){
        cout << "Total dimension of eigenenergy is not correct." << endl;
        cout << "Expected number: " << dim << endl;
        cout << "Actual number: " << index << endl;
    }

    // Sort according to the energy magnitude in ascending order
    sort(eigen_mag.begin(), eigen_mag.end(), Vec_Pair_Double_First_Sort<int>);

    int pos = eigen_mag[0].second;

    for (int j=0; j<evec_basic[pos].size();j++) evec[j] = evec_basic[pos][j];

    // Compute left reduce density matrix
    reduced_density_left_2(evec, size, size/2, reduced_d);

    SelfAdjointEigenSolver<MatrixXcd> density_eigen; // Eigen for reduced density matrix
    density_eigen.compute(reduced_d, false); // Eigenvectors not computed

    // Compute entropy for this state
    for (int j=0; j<density_eigen.eigenvalues().rows();j++){
        double eval = density_eigen.eigenvalues()(j);
        if (abs(eval)>1.0e-15)
        {
            if (eval<0){
                cout << "Density matrix has significant negative eigenvalues." << endl;
                cout << eval << endl;
                abort();
            }
            ent += -eval*log2(eval);
        }
    }

    if (debug){
        cout << "Realization " << local_info.realization_index << " eigenstate " << pos << ":" << endl;
        cout << "Reduced Density Matrix: " << endl;
        complex_matrix_write(reduced_d);
        cout << "Entropy: " << ent << endl;
    }

    model_data_.ent_smallest_var[local_info.J_index][local_info.realization_index] = ent;
}

/*
 * Output averages of mean entropy variance for eigenstate with smallest phase magnitude among all
 * realizations for each J
 */
void FloModelTransition::Ent_smallest_var_out_(AllPara const & parameters, const string& name) {
    const double J_min = parameters.floquet.J_min;
    const double J_max = parameters.floquet.J_max;
    const int J_N = parameters.floquet.J_N;
    const int num_realizations = parameters.generic.num_realizations;
    const int version = parameters.generic.version;
    const int width = parameters.output.width;
    const bool output = parameters.output.filename_output;
    const int size = parameters.generic.size;

    if (model_data_.ent_smallest_var.size() != J_N){
        cout << "Not enough number of J for entropy variance of smallest" <<
                " energy magnitude eigenstate." << endl;
        cout << "Expected Number: " << J_N << endl;
        cout << "Actual number: " << model_data_.ent_smallest_var.size() << endl;
        abort();
    }

    for (int i=0; i< J_N; i++){
        if (model_data_.ent_smallest_var[i].size() != num_realizations){
            cout << "Not enough number of realizations at " << i <<"th J for entropy varaince "
                    << "for smallest eigenstate" << endl;
            cout << "Expected Number: " << num_realizations << endl;
            cout << "Actual Number: " << model_data_.ent_smallest_var[i].size() << endl;
            abort();
        }
    }

    stringstream filename;
    filename << name << ",size=" << size << ",Run=" << num_realizations << ",J_N=" << J_N << ",J_min=" << J_min
            << ",J_max=" << J_max << ",entropy_variance_smallest";
    if (version > 0) filename <<",v" << version;
    filename << ".txt";

    if (output) cout << filename.str() <<endl;

    ofstream fout( filename.str().c_str() );

    for (int i=0; i<J_N;i++){
        double J;
        double mean,sd;
        if (J_N > 1) J = J_min + i * (J_max - J_min)/(J_N-1);
        else J = J_min;

        generic_mean_sd(model_data_.ent_smallest_var[i], mean, sd);
        fout << setw(10) << J << setw(width) << sd << endl;
    }
}


