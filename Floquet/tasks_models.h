#include <map>
#include <vector>
#include "parameters.h"
#include "evol_class.h"

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
 		// Map for task names. The key is its name, and the value is its description
 		map<string, string> tasks_;

 		// Map for model names. The key is its name, and the value is its description
 		map<string, string> models_;

 		// Construct the two maps
 		void Map_Construct_();

 	public:
 		TasksModels() {Map_Construct_();}
 		~TasksModels(){}; 

 		// Look up a name. If it exists, the description is returned; otherwise return
 		// "Not Found"
 		string Look_Up(const string&) const;

 		// According to the name of the task, return the task function pointer
 		task_func Task(const string&) const; 

 		// According to the name of the model, return one model pointer
 		template<class T>
 		void Model(const string&, const AllPara&, EvolMatrix<T>*&);

 		// Print all tasks and their descriptions
 		void Print_Task() const;

 		// Print all models and their descriptions
 		void Print_Model() const;
 };

