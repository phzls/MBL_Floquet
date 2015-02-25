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
 * Compute entry and norm representations of the left and right sigma matrix at the end of chain
 * for floquet systems. This function is to be called under the flo_single_model function. Since
 * there is only one model, the number of realization is always 1.
 */

using namespace std;

void flo_chain_end_sigma_z_under_one_model(const AllPara& parameters,
EvolMatrix<ComplexEigenSolver<MatrixXcd> >* floquet){

	// System Size
	const int size = parameters.generic.size; 

	// Number of realizations.
	const int num_realization = parameters.generic.num_realizations;

	const string model = parameters.generic.model;
	const bool erase = parameters.generic.erase; // Whether erase the evolution matrix
	const bool debug = parameters.generic.debug; // Whether print debug information	

	const int width = parameters.output.width; // Width in output file

	EvolMatrix<ComplexEigenSolver<MatrixXcd> >* floquet_iso; // Possible isolated part
	floquet -> Get_Iso(floquet_iso);
	
	// Isolated part does not exist
	if (floquet_iso == NULL) floquet_iso = floquet;
	
	if (!floquet_iso -> Get_Eigen_Computed()) floquet_iso -> Evol_Diag();

	MatrixXcd sigma_z_left_entry(0,0); // Record entry for leftmost sigma z
	MatrixXd sigma_z_left_norm(0,0); // Record norm for each entry for leftmost sigma z
	MatrixXcd sigma_z_right_entry(0,0); // Record entry for rightmost sigma z
	MatrixXd sigma_z_right_norm(0,0); // Record norm for each entry for rightmost sigma z

	vector<double> ave_phases(1 << size, 0); // Use to store average phases

	ofstream left_entry_out; // Output entry representation for leftmost matrix
	ofstream left_norm_out; // Output norm of each entry representation for leftmost matrix
	ofstream right_entry_out; // Output entry representation for rightmost matrix
	ofstream right_norm_out; // Output norm of each entry representation for rightmost matrix
	ofstream phase_out; // Output phases of the corresponding eigenvectors used in 
						// representing the matrix. The two orders match.
	stringstream base_name_m_l; // Used for leftmost matrix related output
	stringstream base_name_m_r; // Used for rightmost matrix related output
	stringstream base_name_p; // Used for phases related output

	// The first entry is phases of eigenvalues, and for the second pair,
	// the first is its position in sector; the second is its position in 
	// that sector
	vector<pair<double, pair<int, int> > > eval_pos(1 << size);

	// Eigenvectors ordered in terms of phases of eigenvalues
	vector<vector<complex<double> > > evec_basis(1 << size);

	for (int i=0; i< evec_basis.size(); i++) evec_basis[i].resize(1 << size);

	base_name_m_l << floquet_iso -> Repr();
	base_name_m_r << floquet_iso -> Repr();
	base_name_p << floquet_iso -> Repr();

	base_name_m_l<<",Realizations=1,leftmost_sigma_z_under_one_model";
	base_name_m_r<<",Realizations=1,rightmost_sigma_z_under_one_model";
	base_name_p<<",Realizations=1,eval_phases_under_one_model";

	string post_string =  "_entry.txt";
	Of_Construct(left_entry_out, base_name_m_l, post_string, true) ;
	Of_Construct(right_entry_out, base_name_m_r, post_string, true) ;

	post_string = "_norm.txt";
	Of_Construct(left_norm_out, base_name_m_l, post_string, true) ;
	Of_Construct(right_norm_out, base_name_m_r, post_string, true) ;

	post_string = ".txt";
	Of_Construct(phase_out, base_name_p, post_string, true);

	cout << "Prepare eigenvector basis." << endl;
	int eigen_dim = 0; // Number of eigenvalues
	for (int j=0; j< floquet_iso -> eigen.size(); j++){
		eigen_dim += floquet_iso -> eigen[j].eigenvalues().rows();
	}

	if (eval_pos.size() != eigen_dim){
		cout << "Vector size does not match number of eigenvalues." << endl;
		cout << "Vector size: " << eval_pos.size() << endl;
		cout << "Number of eigenvalues: " << eigen_dim << endl;
		abort();
	}

	int index = 0;
	for (int j=0; j < floquet_iso -> eigen.size(); j++){
		for (int k=0; k<floquet_iso -> eigen[j].eigenvalues().rows(); k++){
			eval_pos[index].first = arg(floquet_iso -> eigen[j].eigenvalues()(k));
			eval_pos[index].second.first = j;
			eval_pos[index].second.second = k;
			index ++;
		}
	}

	cout << "Eval finished." << endl;

	if (debug){
		cout << "eigenvalues: " << endl;
		for (int i=0; i<floquet_iso -> eigen.size();i++){
			complex_matrix_write(floquet_iso -> eigen[i].eigenvalues());
		}
		cout << endl;
	}

	// Sort according to the phase magnitude
	sort(eval_pos.begin(), eval_pos.end(), Vec_Pair_Double_First_Sort<pair<int,int> >);

	for (int i=0; i<eval_pos.size(); i++) ave_phases[i] += eval_pos[i].first;

	eigen_dim = 0; // Compute number of eigenvectors
	for (int j=0; j< floquet_iso -> eigen.size(); j++){
		eigen_dim += floquet_iso -> eigen[j].eigenvectors().cols();
	}

	cout << "Total dimension: " << eigen_dim << endl;

	if (evec_basis.size() != eigen_dim){
		cout << "Vector size does not match number of eigenvectors." << endl;
		cout << "Vector size: " << evec_basis.size() << endl;
		cout << "Number of eigenvectors: " << eigen_dim << endl;
		abort();
	}

	evec_to_basic(floquet_iso, evec_basis);
	vector<int> sec_dim = floquet_iso -> Get_Sector_Dim(); // Dimension of each sector

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

	cout << "Construct leftmost sigma z matrix in entry." << endl;
	leftmost_sigma_z_sum(sigma_z_left_entry, evec_basis, "entry");

	if (debug){
		cout << "Entry:" << endl;
		complex_matrix_write(sigma_z_left_entry);
	}

	cout << "Construct leftmost sigma z matrix in norm." << endl;
	leftmost_sigma_z_sum(sigma_z_left_norm, evec_basis, "norm");

	if (debug){
		cout << "Norm:" << endl;
		for (int j=0; j< sigma_z_left_norm.rows(); j++){
			for (int k=0; k<sigma_z_left_norm.cols(); k++){
				cout << sigma_z_left_norm(j,k) << "  ";
			}
			cout << endl;
		}
	}

	cout << endl;

	cout << "Construct rightmost sigma z matrix in entry." << endl;
	rightmost_sigma_z_sum(sigma_z_right_entry, evec_basis, "entry");

	if (debug){
		cout << "Entry:" << endl;
		complex_matrix_write(sigma_z_right_entry);
	}

	cout << "Construct rightmost sigma z matrix in norm." << endl;
	rightmost_sigma_z_sum(sigma_z_right_norm, evec_basis, "norm");

	if (debug){
		cout << "Norm:" << endl;
		for (int j=0; j< sigma_z_right_norm.rows(); j++){
			for (int k=0; k<sigma_z_right_norm.cols(); k++){
				cout << sigma_z_right_norm(j,k) << "  ";
			}
			cout << endl;
		}
	}

	cout << endl;

	if (sigma_z_left_entry.rows() != sigma_z_left_norm.rows()){
		cout << "Two leftmost sigma z matrices have different rows." << endl;
		cout << "Leftmost sigma z entry rows: " << sigma_z_left_entry.rows() << endl;
		cout << "Leftmost sigma z norm rows: " << sigma_z_left_norm.rows() << endl;
		abort();
	}

	if (sigma_z_right_entry.rows() != sigma_z_right_norm.rows()){
		cout << "Two rightmost sigma z matrices have different rows." << endl;
		cout << "Rightmost sigma z entry rows: " << sigma_z_right_entry.rows() << endl;
		cout << "Rightmost sigma z norm rows: " << sigma_z_right_norm.rows() << endl;
		abort();
	}

	if (sigma_z_left_entry.cols() != sigma_z_left_norm.cols()){
		cout << "Two leftmost sigma z matrices have different cols." << endl;
		cout << "Leftmost sigma z entry cols: " << sigma_z_left_entry.cols() << endl;
		cout << "Leftmost sigma z norm cols: " << sigma_z_left_norm.cols() << endl;
		abort();
	}

	if (sigma_z_right_entry.cols() != sigma_z_right_norm.cols()){
		cout << "Two rightmost sigma z matrices have different cols." << endl;
		cout << "Rightmost sigma z entry cols: " << sigma_z_right_entry.cols() << endl;
		cout << "Rightmost sigma z norm cols: " << sigma_z_right_norm.cols() << endl;
		abort();
	}

	// Vectors used to record each row of the matrices
	vector<complex<double> > left_entry_row(sigma_z_left_entry.cols());
	vector<double> left_norm_row(sigma_z_left_norm.cols());

	vector<complex<double> > right_entry_row(sigma_z_right_entry.cols());
	vector<double> right_norm_row(sigma_z_right_norm.cols());

	cout << "Output files." << endl;
	for (int i=0; i< sigma_z_left_entry.rows(); i++){
		for (int j=0; j< sigma_z_left_entry.cols(); j++){
			left_entry_row[j] = sigma_z_left_entry(i,j) / double(num_realization);
			left_norm_row[j] = sigma_z_left_norm(i,j) / num_realization;
		}

		Write_File(left_entry_out, left_entry_row, width);
		Write_File(left_norm_out, left_norm_row, width);
	}

	for (int i=0; i< sigma_z_right_entry.rows(); i++){
		for (int j=0; j< sigma_z_right_entry.cols(); j++){
			right_entry_row[j] = sigma_z_right_entry(i,j) / double(num_realization);
			right_norm_row[j] = sigma_z_right_norm(i,j) / num_realization;
		}

		Write_File(right_entry_out, right_entry_row, width);
		Write_File(right_norm_out, right_norm_row, width);
	}

	for (int i=0; i<ave_phases.size();i++){
		ave_phases[i] /= double(num_realization);
		phase_out << setw(width) << ave_phases[i] << endl;
	}

}