//
// Created by Liangsheng Zhang on 4/6/15.
//

#include <vector>
#include <complex>
#include <iomanip>
#include "flo_eigen_func.h"
#include "methods.h"
#include "eigen_output.h"

using namespace std;

/**
 ** This file contains functions output eigenvectors
 **/

/*
 * evec in eigen_data_ are initialized during running time. Only 1 realization is possible
 */
void FloEigenFunc::Evec_eigen_init_(AllPara const & parameters) {

    const int num_realization = parameters.generic.num_realizations;
    if (num_realization > 1){
        cout << "Only the first realization is used for output eigenvectors." << endl;
    }
}

/*
 * Compute evec. It is done for realization_index 0
 */
void FloEigenFunc::Evec_eigen_compute_(AllPara const & parameters,
                                                 const EvolMatrix<ComplexEigenSolver<MatrixXcd> >* floquet,
                                                 const LocalInfo& local_info) {

    if (local_info.realization_index == 0){
        const bool debug = parameters.generic.debug;

        const int sec_num = floquet -> eigen.size();

        eigen_data_.evec.resize(sec_num);

        for (int i=0; i < sec_num; i++){
            cout << "sector: " << i << endl;
            // Each eigenvector
            const int eigen_num = floquet -> eigen[i].eigenvectors().cols(); // Num of eigenvectors
            const int eigen_dim = floquet -> eigen[i].eigenvectors().rows(); // Dim of each eigenvector
            eigen_data_.evec[i].resize(eigen_num);
            for (int j = 0; j < eigen_num; j++){
                eigen_data_.evec[i][j].resize(eigen_dim);
                for (int k=0; k < eigen_dim; k++){
                    eigen_data_.evec[i][j][k] = floquet -> eigen[i].eigenvectors()(k,j);
                }
            }
        }
    }
}

/*
 * Output all eigenvectors for the first realization. Each sector outputs to a seperate file. Each row
 * is an eigenvector
 */
void FloEigenFunc::Evec_eigen_out_(AllPara const & parameters, const string& name) {
    const double J = parameters.floquet.J;
    const int num_realizations = parameters.generic.num_realizations;
    const int version = parameters.generic.version;
    const int width = parameters.output.width;
    const bool output = parameters.output.filename_output;
    const int size = parameters.generic.size;

    stringstream filename_base;
    filename_base << name << ",size=" << size << ",Run_index=0,J=" << J << ",evec_sector_";

    for (int i=0; i < eigen_data_.evec.size(); i++){
        stringstream filename;
        filename << filename_base.str() << i;
        if (version > 0) filename <<",v" << version;
        filename << ".txt";

        if (output) cout << filename.str() <<endl;

        ofstream fout( filename.str().c_str() );

        for (int j=0; j < eigen_data_.evec[i].size(); j++){
            for (int k=0; k < eigen_data_.evec[i][j].size(); k++)
                fout << setw(width) << complex_write_return(eigen_data_.evec[i][j][k]);
            fout << endl;
        }
    }
}


