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

	// Compute level statistics for Floquet system
	string task_name2;
	string task_type2;
	task_func task_function2;

	task_name2 = "Flo Rightmost Sigma_z";
	task_type2 = "Floquet";
	task_function2 = &flo_rightmost_sigma_z;
	Task_Map_Insert(task_name2, task_type2, task_function2);

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
