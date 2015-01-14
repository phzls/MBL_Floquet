#include "eigen_output.h"

using namespace std;
using namespace Eigen;

/**
 ** Implement functions in eigen_output.h
 **/

void complex_matrix_write(const MatrixXcd& M){
	for (int i=0; i<M.rows(); i++){
		for (int j=0; j<M.cols(); j++){
			complex_write(M(i,j));
		}
		cout << endl;
	}
}