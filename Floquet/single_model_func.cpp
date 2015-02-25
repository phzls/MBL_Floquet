#include <map>
#include <string>
#include "single_model_func.h"
#include "tasks_models.h"
#include "parameters.h"

using namespace std;

/**
 ** This file implements the function which calls other functions using a single evolution
 ** model matrix.
 **/

extern TasksModels tasks_models; // Record all the tasks and methods. Defined in main.

// Pointer to single model functions for floquet systems
typedef void (*flo_single_model)(const AllPara&, EvolMatrix<ComplexEigenSolver<MatrixXcd> >*);

// Pointer to single model functions for Hamiltonian systems
typedef void (*ham_single_model)(const AllPara&, EvolMatrix<EigenSolver<MatrixXd> >*);

// Initialize flo_single_model_func_map
void flo_single_model_func_map_initialize(map<string, flo_single_model>&);

void single_model(const AllPara& parameters){
	// Map that links strings to flo_single_model functions
	map<string, flo_single_model> flo_single_model_func_map;

	// Map that links strings to ham_single_model functions
	map<string, ham_single_model> ham_single_model_func_map;

	// Map that determines which functions are computed
	map<string, bool> single_model_bool_map;

	single_model_bool_map = parameters.single_model.single_model_compute;

	flo_single_model_func_map_initialize(flo_single_model_func_map);

	// Check the number of function
	if ( single_model_bool_map.size() != flo_single_model_func_map.size() +
		 ham_single_model_func_map.size()){

		cout << "Number of registered functions in single_model is not the same as"
			 << " number of registered functions that can be called." << endl;
		cout << "Functions in parameters:" << endl;
		for (map<string, bool>::iterator it = single_model_bool_map.begin(); 
			it != single_model_bool_map.end(); it++){
			cout << it -> first << endl;
		}
		cout << "Total Number: " << single_model_bool_map.size() << endl;

		cout << "Functions that can be called for floquet:" << endl;
		for (map<string, flo_single_model>::iterator it = flo_single_model_func_map.begin(); 
			it != flo_single_model_func_map.end(); it ++){
			cout << it -> first << endl;
		}
		cout << "Total Number:" << flo_single_model_func_map.size() << endl;

		cout << "Functions that can be called for Hamiltonian:" << endl;
		for (map<string, ham_single_model>::iterator it = ham_single_model_func_map.begin(); 
			it != ham_single_model_func_map.end(); it ++){
			cout << it -> first << endl;
		}
		cout << "Total Number:" << ham_single_model_func_map.size() << endl;

		abort();
	}

	EvolMatrix<ComplexEigenSolver<MatrixXcd> >* floquet;

	const string model = parameters.generic.model;

	tasks_models.Model(model, parameters, floquet);

	cout << "Diagonalize Evolution Operator." << endl;
	// For safety now, the evolution matrix and its eigenvectors are both kept
	floquet -> Evol_Para_Init();
	floquet -> Evol_Construct();
	floquet -> Evol_Diag();

	// Call functions
	for (map<string, bool>::iterator it = single_model_bool_map.begin(); 
	it != single_model_bool_map.end(); it++){
		if (it -> second){
			if (it -> first.find("Flo") != string::npos){
				flo_single_model_func_map[it -> first](parameters, floquet);
			}
			else if (it -> first.find("Ham") != string::npos){
				cout << "Not implemented yet" << endl;
				abort();
			}
			else{
				cout << "According to the name " << it -> first << " it is not clear "
					 << "its evolution system." << endl;
				abort();
			}
		}
	}

	delete floquet;
	floquet = NULL;
}

void flo_single_model_func_map_initialize(map<string, flo_single_model>& func_map){
	map<string, flo_single_model>::iterator func_it;

	// sigma matrices at both ends of the spin chain
	string name1 = "Flo Chain End Sigma Z";
	flo_single_model func1 = flo_chain_end_sigma_z_under_one_model;

	// Make sure the name has not been used before
	func_it = func_map.find(name1);
	if (func_it != func_map.end()){
		cout << name1 << " for single_model has appeared before." << endl;
		abort();
	}

	func_map[name1] = func1;

	// simple markov evolution with possible multiple initial parameters
	string name2 = "Flo Evolution Simple Markov";
	flo_single_model func2 = flo_evolution_simple_markov_under_one_model;

	// Make sure the name has not been used before
	func_it = func_map.find(name2);
	if (func_it != func_map.end()){
		cout << name2 << " for single_model has appeared before." << endl;
		abort();
	}

	func_map[name2] = func2;

}