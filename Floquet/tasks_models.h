#include <map>
#include <utility>
#include <iostream>
#include <string>
#include "parameters.h"
#include "evol_class.h"
#include "model_func.h"

using namespace std;

/**
 ** This file contains a class which associates names with functions in task_func.h and 
 ** derived classes from evol_class.h. Note it is assumed that all functions in task_func.h 
 ** just take parameter class as argument.
 **/

 typedef void (*task_func)(const AllPara&);

 class TasksModels
 {	
 	private:
 		// Map for task names. The key is its name, and the value is the task type and 
 		// corresponding function call
 		map<string, pair<string, task_func> > tasks_;

 		// Map for model names. The key is its name, and the value is the model type, which
 		// should be the same as the type obtained from the model class, and the corresponding 
 		// function call
 		map<string, pair<string, ModelFunc*> > models_;

 		// The recorded task type
 		string task_type_;

 		// Construct the two maps
 		void Map_Construct_();

 		// Insert task pair in map
 		void Task_Map_Insert(const string&, const string&, const task_func&);

 		// Insert model pair in map
 		void Model_Map_Insert(const string&, const string&, ModelFunc* const &);

 	public:
 		TasksModels() {Map_Construct_();}
 		~TasksModels(); 

 		// Look up a name in tasks. If it exists, return true, otherwise false
 		bool Task_Look_Up(const string&) const;

 		// Look up a name in models. If it exists, return true, otherwise false
 		bool Model_Look_Up(const string&) const;


 		// According to the name of the task, return the task function pointer
 		task_func Task(const string&); 

 		// According to the name of the model, return one model pointer
 		template<class T>
 		void Model(const string&, const AllPara&, EvolMatrix<T>*&) const;

 		// Print all tasks names
 		void Print_Task() const;

 		// Print all models names
 		void Print_Model() const;
 };

template <class T>
void TasksModels::Model(const string& model_name, const AllPara& parameters, 
EvolMatrix<T>*& model) const {
	map<string, pair<string, ModelFunc*> >::const_iterator it;
	it = models_.find(model_name);

	if (it == models_.end()){
		cout << "The model desired is not found." << endl;
		Print_Model();
		abort();
	}
	else{
		if (it -> second.first != it -> second.second -> Type()){
			cout << "Model type is not consistent." << endl;
			abort();
		}

		if (it -> second.first.find(task_type_) == string::npos){
			cout << "Model type and task type are not consistent." << endl;
			abort();
		}

		*(it -> second.second)(parameters, model);
	}
}

