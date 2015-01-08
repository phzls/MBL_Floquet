#include <map>
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
 		// Map for task names. The key is its name, and the value is the corresponding 
 		// function call
 		map<string, task_func> tasks_;

 		// Map for model names. The key is its name, and the value is the corresponding
 		// function call
 		map<string, ModelFunc*> models_;

 		// Construct the two maps
 		void Map_Construct_();

 		// The model type that may be induced by the method. It can be "Floquet", 
 		// "Hamiltonian" or "All"
 		string method_type_();

 		// Insert task pair in map
 		void Task_Map_Insert(const string&, const task_func&);

 		// Insert model pair in map
 		void Model_Map_Insert(const string&, ModelFunc* const &);

 	public:
 		TasksModels() {Map_Construct_();}
 		~TasksModels(); 

 		// Look up a name. If it exists, return true, otherwise false
 		bool Look_Up(const string&) const;

 		// According to the name of the task, return the task function pointer
 		task_func Task(const string&) const; 

 		// According to the name of the model, return one model pointer
 		template<class T>
 		void Model(const string&, const AllPara&, EvolMatrix<T>*&);

 		// Print all tasks names
 		void Print_Task() const;

 		// Print all models names
 		void Print_Model() const;
 };

