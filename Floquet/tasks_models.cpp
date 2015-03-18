#include <iostream>
#include <utility>
#include <iostream>
#include "tasks_models.h"
#include "task_func.h"

using namespace std;

void TasksModels::Map_Construct_(){
	// Compute level statistics for Floquet system
	string task_name1;
	string task_type1;
	task_func task_function1;

	task_name1 = "Flo Level";
	task_type1 = "Floquet";
	task_function1 = &flo_level;
	Task_Map_Insert(task_name1, task_type1, task_function1);

	// Compute the representation of rightmost sigma_z operator for a floquet
	// system in the basis of eigenstates
	string task_name2;
	string task_type2;
	task_func task_function2;

	task_name2 = "Flo Rightmost Sigma_z";
	task_type2 = "Floquet";
	task_function2 = &flo_rightmost_sigma_z;
	Task_Map_Insert(task_name2, task_type2, task_function2);

	// Compute time evolution for a Floquet system of pure states
	string task_name3;
	string task_type3;
	task_func task_function3;

	task_name3 = "Flo Evolution";
	task_type3 = "Floquet";
	task_function3 = &flo_evolution;
	Task_Map_Insert(task_name3, task_type3, task_function3);

	// Compute Markov time evolution for a Floquet system
	string task_name4;
	string task_type4;
	task_func task_function4;

	task_name4 = "Flo Simple Markov Evolution";
	task_type4 = "Floquet";
	task_function4 = &flo_evolution_simple_markov;
	Task_Map_Insert(task_name4, task_type4, task_function4);

	// Compute the representation of leftmost sigma_z operator for a floquet
	// system in the basis of eigenstates
	string task_name5;
	string task_type5;
	task_func task_function5;

	task_name5 = "Flo Leftmost Sigma_z";
	task_type5 = "Floquet";
	task_function5 = &flo_leftmost_sigma_z;
	Task_Map_Insert(task_name5, task_type5, task_function5);

	// Compute the representation of leftmost and rightmost sigma_z operator for 
	// a floquet system in the basis of eigenstates
	string task_name6;
	string task_type6;
	task_func task_function6;

	task_name6 = "Flo Chain End Sigma_z";
	task_type6 = "Floquet";
	task_function6 = &flo_chain_end_sigma_z;
	Task_Map_Insert(task_name6, task_type6, task_function6);

	// Compute Markov time evolution for a Floquet system with one model and possible
	// multiple sets of initial parameters
	string task_name7;
	string task_type7;
	task_func task_function7;

	task_name7 = "Flo Simple Markov Evolution One Model";
	task_type7 = "Floquet";
	task_function7 = &flo_evolution_simple_markov_one_model;
	Task_Map_Insert(task_name7, task_type7, task_function7);

	// Implement single model calculations, in which all calculations use one
	// realization of the model
	string task_name8;
	string task_type8;
	task_func task_function8;

	task_name8 = "Single Model";
	task_type8 = "All";
	task_function8 = &single_model;
	Task_Map_Insert(task_name8, task_type8, task_function8);

	// Compute time evolution for a Floquet system of density matrices
	string task_name9;
	string task_type9;
	task_func task_function9;

	task_name9 = "Flo Evolution Density";
	task_type9 = "Floquet";
	task_function9 = &flo_evolution_density;
	Task_Map_Insert(task_name9, task_type9, task_function9);

	// Random Floquet Operator
	string model_name1;
	string model_type1;
	ModelFunc* model_function1;

	model_name1 = "Random Flo";
	model_type1 = "Random Floquet";
	model_function1 = new FloEvolRandomFunc();
	Model_Map_Insert(model_name1, model_type1, model_function1);

	// Random Rotation Floquet Operator
	string model_name2;
	string model_type2;
	ModelFunc* model_function2;

	model_name2 = "Random Rotation Flo";
	model_type2 = "Random Rotation Floquet";
	model_function2 = new FloEvolRandomRotationFunc();
	Model_Map_Insert(model_name2, model_type2, model_function2);

	// XXZ Floquet Operator
	string model_name3;
	string model_type3;
	ModelFunc* model_function3;

	model_name3 = "XXZ Flo";
	model_type3 = "XXZ Floquet";
	model_function3 = new FloEvolXXZFunc();
	Model_Map_Insert(model_name3, model_type3, model_function3);

	// Inter Random Floquet Operator
	string model_name4;
	string model_type4;
	ModelFunc* model_function4;

	model_name4 = "Inter Random Flo";
	model_type4 = "Inter Random Floquet";
	model_function4 = new FloEvolInterRandomFunc();
	Model_Map_Insert(model_name4, model_type4, model_function4);

	// Markov Inter Random Floquet Operator
	string model_name5;
	string model_type5;
	ModelFunc* model_function5;

	model_name5 = "Markov Inter Random Flo";
	model_type5 = "Markov Inter Random Floquet";
	model_function5 = new FloEvolMarkovInterRandomFunc();
	Model_Map_Insert(model_name5, model_type5, model_function5);

	// Markov Inter Random Both Floquet Operator
	string model_name6;
	string model_type6;
	ModelFunc* model_function6;

	model_name6 = "Markov Inter Random Both Flo";
	model_type6 = "Markov Inter Random Both Floquet";
	model_function6 = new FloEvolMarkovInterRandomBothFunc();
	Model_Map_Insert(model_name6, model_type6, model_function6);

	// Markov Inter Random Both Floquet Operator which couples to the bath through a
	// sigma_x term
	string model_name7;
	string model_type7;
	ModelFunc* model_function7;

	model_name7 = "Markov Inter Random Both X Flo";
	model_type7 = "Markov Inter Random Both X Floquet";
	model_function7 = new FloEvolMarkovInterRandomBothXFunc();
	Model_Map_Insert(model_name7, model_type7, model_function7);

	// Random XXZ Floquet Operator
	string model_name8;
	string model_type8;
	ModelFunc* model_function8;

	model_name8 = "XXZ Random Flo";
	model_type8 = "XXZ Random Floquet";
	model_function8 = new FloEvolXXZRandomFunc();
	Model_Map_Insert(model_name8, model_type8, model_function8);

	// Markov XXZ Random Both Floquet Operator which couples to the bath through a
	// sigma_x term
	string model_name9;
	string model_type9;
	ModelFunc* model_function9;

	model_name9 = "Markov XXZ Random Both X Flo";
	model_type9 = "Markov XXZ Random Both X Floquet";
	model_function9 = new FloEvolMarkovXXZRandomBothXFunc();
	Model_Map_Insert(model_name9, model_type9, model_function9);
}

void TasksModels::Task_Map_Insert(const string& task_name, const string& task_type, 
const task_func& task_function){
	pair<string, task_func> task_type_op;

	if (tasks_.find(task_name) == tasks_.end()){ // No duplicate names
		task_type_op = make_pair(task_type, task_function);
		tasks_[task_name] = task_type_op;
	}
	else{
		cout << "Duplicate names for tasks." << endl;
		abort();
	}
}

void TasksModels::Model_Map_Insert(const string& model_name, const string& model_type, 
ModelFunc* const & model_function){
	pair<string, ModelFunc*> model_type_op;

	if (models_.find(model_name) == models_.end()){ // No duplicate names
		model_type_op = make_pair(model_type, model_function);
		models_[model_name] = model_type_op;
	}
	else{
		cout << "Duplicate names for models." << endl;
		abort();
	}
}

bool TasksModels::Task_Look_Up(const string& task_name) const {
	if (tasks_.find(task_name) == tasks_.end()){
		return false;
	}
	else return true;
}

bool TasksModels::Model_Look_Up(const string& model_name) const {
	if (models_.find(model_name) == models_.end()){
		return false;
	}
	else return true;
}

task_func TasksModels::Task(const string& task_name) {
	map<string, pair<string, task_func> >::const_iterator it;
	it = tasks_.find(task_name);
	if (it == tasks_.end()){
		cout << "The task desired is not found." << endl;
		Print_Task();
		abort();
	}
	else{
		task_type_ = it -> second.first;
		return it -> second.second;
	}
}

void TasksModels::Print_Task() const {
	cout << "The tasks are: " << endl;
	map<string, pair<string, task_func> >::const_iterator it;

	for (it = tasks_.begin(); it != tasks_.end(); it++){
		cout << "Name: " << it -> first <<"  Type: " << it -> second.first << endl;
	}
	cout << endl;
}

void TasksModels::Print_Model() const {
	cout << "The models are: " << endl;
	map<string, pair<string, ModelFunc*> >::const_iterator it;

	for (it = models_.begin(); it != models_.end(); it++){
		cout << "Name: " << it -> first <<"  Type: " << it -> second.first << endl;
	}
	cout << endl;
}

TasksModels::~TasksModels(){
	map<string, pair<string, ModelFunc*> >::iterator it;

	for (it = models_.begin(); it != models_.end(); it++){
		if (it -> second.second != NULL){
			delete it -> second.second;
			it -> second.second = NULL;
		}
	}
}
