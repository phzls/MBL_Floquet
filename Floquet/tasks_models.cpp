#include <iostream>
#include <utility>
#include "tasks_models.h"
#include "task_func.h"

using namespace std;

void TasksModels::Map_Construct_(){
	// Compute level statistics for Floquet system
	pair<string, task_func> task_name_op1;
	string task_name1;
	task_func task_function1;

	task_name1 = "Flo Level";
	task_function1 = flo_level;
	task_name_op1 = make_pair(task_name1, task_function1);
	tasks_.insert(task_name_op1);

	// Random Floquet Operator
	pair<string, ModelFunc*> model_name_op1;
	string model_name1;
	ModelFunc* model_function1;

	model_name1 = "Random Flo";
	model_function1 = new FloEvolRandomFunc();
	model_name_op1 = make_pair(model_name1, model_function1);
	models_.insert(model_name_op1);

	// Random Rotation Floquet Operator
	pair<string, ModelFunc*> model_name_op2;
	string model_name2;
	ModelFunc* model_function2;

	model_name2 = "Random Flo Rotation";
	model_function2 = new FloEvolRandomRotationFunc();
	model_name_op2 = make_pair(model_name2, model_function2);
	models_.insert(model_name_op2);
}
