#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <iomanip>
#include "flo_evol_model.h"
#include "methods.h"
#include "output_func.h"
#include "tasks_models.h"
#include "parameters.h"
#include "sort.h"
#include "eigen_output.h"

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

	const int width = parameters.output.width; // Width in output file

	bool output_init = false; // Whether output filenames have been initialized

	EvolMatrix<ComplexEigenSolver<MatrixXcd> >* floquet;

	MatrixXcd sigma_z_entry(0,0); // Record entry
	MatrixXd sigma_z_norm(0,0); // Record norm for each entry

	vector<double> ave_phases(1 << size, 0); // Use to store average phases

	ofstream entry_out; // Output entry representation
	ofstream norm_out; // Output norm of each entry representation
	ofstream phase_out; // Output phases of the corresponding eigenvectors used in 
						// representing the matrix. The two orders match.
	stringstream base_name_m; // Used for matrix related output
	stringstream base_name_p; // Used for phases related output

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
			base_name_m << floquet -> Repr();
			base_name_p << floquet -> Repr();

			base_name_m<<",Realizations="<<num_realization<<",right_most_sigma_z";
			base_name_p<<",Realizations="<<num_realization<<",eval_phases";

			string post_string =  "_entry.txt";
			Of_Construct(entry_out, base_name_m, post_string, true) ;

			post_string = "_norm.txt";
			Of_Construct(norm_out, base_name_m, post_string, true) ;

			post_string = ".txt";
			Of_Construct(phase_out, base_name_p, post_string, true);

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
			for (int k=0; k<floquet -> eigen[j].eigenvalues().rows(); k++){
				eval_pos[index].first = arg(floquet -> eigen[j].eigenvalues()(k));
				eval_pos[index].second.first = j;
				eval_pos[index].second.second = k;
				index ++;
			}
		}

		cout << "Eval finished." << endl;

		// Sort according to the phase magnitude
		sort(eval_pos.begin(), eval_pos.end(), Vec_Pair_Double_First_Sort<pair<int,int> >);

		for (int i=0; i<eval_pos.size(); i++) ave_phases[i] += eval_pos[i].first;

		eigen_dim = 0; // Compute number of eigenvectors
		for (int j=0; j< floquet -> eigen.size(); j++){
			eigen_dim += floquet -> eigen[j].eigenvectors().cols();
		}

		cout << "Total dimension: " << eigen_dim << endl;

		if (evec_basis.size() != eigen_dim){
			cout << "Vector size does not match number of eigenvectors." << endl;
			cout << "Vector size: " << evec_basis.size() << endl;
			cout << "Number of eigenvectors: " << eigen_dim << endl;
			abort();
		}

		evec_to_basic(floquet, evec_basis);
		vector<int> sec_dim = floquet -> Get_Sector_Dim(); // Dimension of each sector

		// Now sec_dim[i] computes the total number of vectors up to sector i
		for (int j=1; j<sec_dim.size();j++) sec_dim[j] += sec_dim[j-1];

		// This vector records the original position (in the whole space) of the eigenvector
		// that is currently held at every site
		vector<int> index_set(evec_basis.size());
		for (int j=0; j<index_set.size();j++) index_set[j] = j;

		// We assume originally evec_basis[i] holds eigenvector whose position in the whole
		// space is also i
		for (int j=0; j<eval_pos.size();j++){
			int sector = eval_pos[j].second.first;
			int sec_pos = eval_pos[j].second.second;

			int pos; // The position of the vector recorded in eval_pos in the whole space
			if (sector == 0) pos = 0;
			else pos = sec_dim[sector-1];
			pos += sec_pos; 

			if (pos >= index_set.size()){
				cout << "The global position of eigenvector from eval_pos[" << j << "]"
					 << " is too large." << endl;
				abort();
			}

			int swap_pos = pos; // The position to swap vector
			while (index_set[swap_pos] != pos ) swap_pos = index_set[swap_pos];

			evec_basis[j].swap(evec_basis[swap_pos]);

			index_set[swap_pos] = index_set[j];
			index_set[j] = pos;
		}

		cout << "Construct rightmost sigma z matrix in entry." << endl;
		rightmost_sigma_z_sum(sigma_z_entry, evec_basis, "entry");

		if (debug){
			cout << "Entry:" << endl;
			complex_matrix_write(sigma_z_entry);
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

	for (int i=0; i<ave_phases.size();i++){
		ave_phases[i] /= double(num_realization);
		phase_out << setw(width) << ave_phases[i] << endl;
	}

}