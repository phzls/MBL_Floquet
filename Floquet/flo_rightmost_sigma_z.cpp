#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include "flo_evol_model.h"
#include "methods.h"
#include "output_func.h"
#include "tasks_models.h"
#include "parameters.h"
#include "sort.h"

/*
 * Compute entry and norm representations of the right most sigma matrix at the end of chain
 * for floquet systems.
 */

using namespace std;

extern TasksModels tasks_models; // Record all the tasks and methods. Defined in main.

void flo_rightmost_sigma_z(const AllPara& parameters){

	// System Size
	const int size = parameters.generic.size; 
	// Number of realizations.
	const int num_realization = parameters.generic.num_realizations;
	const string model = parameters.generic.model;
	const bool erase = parameters.generic.erase; // Whether erase the evolution matrix
	const bool debug = parameters.generic.debug; // Whether print debug information

	const double tau = parameters.floquet.tau; // Time step, which seems not used here
	const double J = parameters.floquet.J; // Coupling strength

	const double angle_min = parameters.floquet_random.angle_min; // Minimum angle
	const double angle_sup = parameters.floquet_random.angle_sup; // Supreme angle	

	const int width = parameters.output.width; // Width in output file

	bool output_init = false; // Whether output filenames have been initialized

	EvolMatrix<ComplexEigenSolver<MatrixXcd> >* floquet;

	MatrixXcd sigma_z_entry(0,0); // Record entry
	MatrixXd sigma_z_norm(0,0); // Record norm for each entry

	ofstream entry_out; // Output entry representation
	ofstream norm_out; // Output norm of each entry representation
	stringstream base_name;

	// The first entry is phases of eigenvalues, and for the second pair,
	// the first is its position in sector; the second is its position in 
	// that sector
	vector<pair<double, pair<int, int> > > eval_pos(1 << size);

	// Eigenvectors ordered in terms of phases of eigenvalues
	vector<vector<complex<double> > > evec_basis(1 << size);

	for (int i=0; i< evec_basis.size(); i++) evec_basis[i].resize(1 << size);

	for (int i=0; i< num_realization; i++){

		cout << i << "th run:" << endl;

		cout << "Initialize Model." << endl;
		tasks_models.Model(model, parameters, floquet);

		if (!output_init){
			base_name << floquet -> Type()<<"_L=" << size << ",tau=" << tau
						  << ",Realizations=" << num_realization << ",J="<< J;

			int found = floquet -> Type().find("Rotation");
			if (found != string::npos){
				// It is the random rotation floquet
				base_name << ",angle_min=" << angle_min <<",angle_sup=" << angle_sup;
			}

			base_name<<",right_most_sigma_z";

			string post_string =  "_entry.txt";
			Of_Construct(entry_out, base_name, post_string, true) ;

			post_string = "_norm.txt";
			Of_Construct(norm_out, base_name, post_string, true) ;

			output_init = true;
		}

		cout << "Diagonalize Evolution Operator." << endl;
		floquet -> Evol_Para_Init();
		floquet -> Evol_Construct();
		floquet -> Evol_Diag();
		if (erase) floquet -> Evol_Erase();

		cout << "Prepare eigenvector basis." << endl;
		int eigen_dim = 0; // Number of eigenvalues
		for (int j=0; j< floquet -> eigen.size(); j++){
			eigen_dim += floquet -> eigen[j].eigenvalues().rows();
		}

		if (eval_pos.size() != eigen_dim){
			cout << "Vector size does not match number of eigenvalues." << endl;
			cout << "Vector size: " << eval_pos.size() << endl;
			cout << "Number of eigenvalues: " << eigen_dim << endl;
			abort();
		}

		int index = 0;
		for (int j=0; j < floquet -> eigen.size(); j++){
			for (int k=0; k<floquet -> eigen[j].eigenvalues.rows(); k++){
				eval_pos[index].first = arg(floquet -> eigen[j].eigenvalues()(k));
				eval_pos[index].second.first = j;
				eval_pos[index].second.second = k;
			}
		}

		// Sort according to the phase magnitude
		sort(eval_pos.begin(), eval_pos.end(), Vec_Pair_Double_Int_First_Sort);

		eigen_dim = 0; // Compute number of eigenvectors
		for (int j=0; j< floquet -> eigen.size(); j++){
			eigen_dim += floquet -> eigen[j].eigenvectors().cols();
		}

		if (evec_basis.size() != eigen_dim){
			cout << "Vector size does not match number of eigenvectors." << endl;
			cout << "Vector size: " << evec_basis.size() << endl;
			cout << "Number of eigenvectors: " << eigen_dim << endl;
			abort();
		}

		for (int j=0; j<evec_basis.size();j++){
			int sector = eval_pos[j].second.first;
			int sec_pos = eval_pos[j].second.second

			for (int k=0; k < evec_basis[j].size(); k++){
				evec_basis[j][k] = floquet -> eigen[sector].eigenvectors()(k, sec_pos);
			}
		}

		cout << "Construct rightmost sigma z matrix in entry." << endl;
		rightmost_sigma_z_sum(sigma_z_entry, evec_basis, "entry");

		if (debug){
			cout << "Entry:" << endl;
			for (int j=0; j< sigma_z_entry.rows(); j++){
				for (int k=0; k<sigma_z_entry.cols(); k++){
					cout << real(sigma_z_entry(j,k));
					if (imag(sigma_z_entry(j,k))<0){
						cout << imag(sigma_z_entry(j,k)) <<"j  ";
					}
					else cout << "+" << imag(sigma_z_entry(j,k)) << "j  ";
				}
				cout << endl;
			}
		}

		cout << "Construct rightmost sigma z matrix in norm." << endl;
		rightmost_sigma_z_sum(sigma_z_norm, evec_basis, "norm");

		if (debug){
			cout << "Norm:" << endl;
			for (int j=0; j< sigma_z_norm.rows(); j++){
				for (int k=0; k<sigma_z_norm.cols(); k++){
					cout << sigma_z_norm(j,k) << "  ";
				}
				cout << endl;
			}
		}

		cout << endl;

		delete floquet;
		floquet = NULL;
	}

	if (sigma_z_entry.rows() != sigma_z_norm.rows()){
		cout << "Two sigma z matrices have different rows." << endl;
		cout << "sigma z entry rows: " << sigma_z_entry.rows() << endl;
		cout << "sigma z norm rows: " << sigma_z_norm.rows() << endl;
		abort();
	}

	if (sigma_z_entry.cols() != sigma_z_norm.cols()){
		cout << "Two sigma z matrices have different cols." << endl;
		cout << "sigma z entry cols: " << sigma_z_entry.cols() << endl;
		cout << "sigma z norm cols: " << sigma_z_norm.cols() << endl;
		abort();
	}

	// Vectors used to record each row of the tow matrices
	vector<complex<double> > entry_row(sigma_z_entry.cols());
	vector<double> norm_row(sigma_z_norm.cols());

	cout << "Output files." << endl;
	for (int i=0; i< sigma_z_entry.rows(); i++){
		for (int j=0; j< sigma_z_entry.cols(); j++){
			entry_row[j] = sigma_z_entry(i,j) / double(num_realization);
			norm_row[j] = sigma_z_norm(i,j) / num_realization;
		}

		Write_File(entry_out, entry_row, width);
		Write_File(norm_out, norm_row, width);
	}

}