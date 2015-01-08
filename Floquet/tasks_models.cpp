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
	task_function1 = flo_level;
	Task_Map_Insert(task_name1, task_type1, task_function1);

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
		abort();
	}
	else{
		task_type_ = it -> second.first;
		return it -> second.second;
	}
}
