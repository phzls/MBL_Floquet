#include <vector>
#include <Eigen/Dense>
#include "evol_class.h"
#include "transition.h"

using namespace std;
using namespace Eigen;

template<class T1, class T2>
void evec_to_basic(EvolMatrix<T1>* const & evol, vector<vector<T2> >& evec){
	if (evol -> Eigen_Type() == "Basic"){
		int index = 0;
		for (int i=0; i<evol -> eigen.size(); i++){
			if (evol -> eigen[i].eigenvectors().rows() != evec[0].size()){
				cout << "The length of eigenvector in sector " << i << " is incompatible." << endl;
				abort();
			}

			for (int j=0; j < evol -> eigen[i].eigenvectors().cols(); j++){
				for (int k=0; k < evol -> eigen[i].eigenvectors().rows(); k++){
					evec[index][k] = evol -> eigen[i].eigenvectors()(k,j);
				}
				index ++; 
			}
		}
	}

	else if (evol -> Eigen_Type() == "Parity"){
		// Assume only 2 sectors, 0 is even and 1 is odd
		// In the evec, even eigenvectors come first
		if (evol -> eigen.size() != 2){
			cout << "The position of even and odd sector in parity basis is not clear." << endl;
			abort();
		}

		TransitionMatrix transition;
		evol -> Transition_Compute(transition, "Basic_Parity");

		if (evec[0].size() != transition.Matrix("Basic_Even").rows()){
			cout << "The dimension of eigenvector acceptor is not correct." << endl;
			cout << "Assumed vector length: " << evec[0].size() << endl;
			cout << "Would obtained vector length: " << transition.Matrix("Basic_Even").rows() 
				 << endl;
			abort();
		}

		Matrix<T2, Dynamic, 1> single_evec;

		int index = 0;
		for (int i=0; i< evol -> eigen[0].eigenvectors().cols(); i++){
			single_evec = transition.Matrix("Basic_Even") * evol -> eigen[0].eigenvectors().col(i);

			for (int j=0; j< single_evec.size();j++)
				evec[index][j] = single_evec(j);

			index ++;
		}

		for (int i=0; i< evol -> eigen[1].eigenvectors().cols(); i++){
			single_evec = transition.Matrix("Basic_Odd") * evol -> eigen[1].eigenvectors().col(i);
			for (int j=0; j< single_evec.size();j++)
				evec[index][j] = single_evec(j);

			index ++;
		}
	}

	else{
		cout << "Eigen type " << evol -> Eigen_Type() << " unknown." << endl;
		abort();
	}
}