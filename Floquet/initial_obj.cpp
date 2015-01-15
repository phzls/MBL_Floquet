#include <iostream>
#include <map>
#include "initial_obj.h"

using namespace std;

init_func InitObj::Init_Func(const string& init_func_name) const {
	map<string, init_func>::const_iterator it = map_init_.find(init_func_name);

	if (it == map_init_.end()){
		cout << "The requested init_func does not exist." << endl;
		Print();
		cout << "Requested function: " << init_func_name << endl;
		abort();
	}
	else return it -> second;
}

void InitObj::Print() const {
	map<string, init_func>::const_iterator it;

	cout << "Initial State Construction Functions: "<< endl;
	for (it = map_init_.begin(); it != map_init_.end(); it++){
		cout << it -> first << endl;
	}
}

void InitObj::map_init_(){
	map<string, init_func>::const_iterator it;

	string name1 = "Random Product";
	init_func func1 = random_product;

	it = init_map_.find(name1);
	if (it != init_map_.end()){
		cout << "init_func " << name1 << " already exists." << endl;
		abort();
	}
	init_map_[name1] = func1;

	string name2 = "Product Random";
	init_func func2 = product_random;

	it = init_map_.find(name2);
	if (it != init_map_.end()){
		cout << "init_func " << name2 << " already exists." << endl;
		abort();
	}
	init_map_[name2] = func2;

}

