#include <iostream>
#include <utility>
#include <iostream>
#include "tasks_models.h"
#include "task_func.h"

using namespace std;

void TasksModels::Map_Construct_(){
	// Compute level statistics for Floquet system
	string task_name1;
	task_func task_function1;

	task_name1 = "Flo Level";
	task_function1 = flo_level;
	Task_Map_Insert(task_name1, task_function1);

	// Random Floquet Operator
	pair<string, ModelFunc*> model_name_op1;
	string model_name1;
	ModelFunc* model_function1;

	model_name1 = "Random Flo";
	model_function1 = new FloEvolRandomFunc();
	Model_Map_Insert(model_name1, model_function1);

	// Random Rotation Floquet Operator
	pair<string, ModelFunc*> model_name_op2;
	string model_name2;
	ModelFunc* model_function2;

	model_name2 = "Random Flo Rotation";
	model_function2 = new FloEvolRandomRotationFunc();
	model_name_op2 = make_pair(model_name2, model_function2);
	models_.insert(model_name_op2);
}

void TasksModels::Task_Map_Insert(const string& task_name, const task_func& task_function){
	pair<string, task_func> task_name_op;
	if (tasks_.find(task_name) == tasks_.end()){ // No duplicate names
		task_name_op = make_pair(task_name, task_function);
		tasks_.insert(task_name_op);
	}
	else{
		cout << "Duplicate names for tasks." << endl;
		abort();
	}
}

void TasksModels::Model_Map_Insert(const string& model_name, ModelFunc* const & model_function){
	pair<string, ModelFunc*> model_name_op;
	if (models_.find(model_name) == models_.end()){ // No duplicate names
		model_name_op = make_pair(model_name, model_function);
		models_.insert(model_name_op);
	}
	else{
		cout << "Duplicate names for models." << endl;
		abort();
	}
}
