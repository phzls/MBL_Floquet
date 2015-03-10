#include <iostream>
#include <map>
#include <sstream>
#include "initial_obj.h"

using namespace std;

InitInfo::InitInfo(const InitInfo& init_info){
	size = init_info.size;
	dim = init_info.dim;
	norm_delta = init_info.norm_delta;
	debug = init_info.debug;
	complex_eigen = init_info.complex_eigen;
	real_eigen = init_info.real_eigen;
	leftmost_spin_z_index = init_info.leftmost_spin_z_index;
	multi_ini_para = init_info.multi_ini_para;
	multi_ini_para_num = init_info.multi_ini_para_num;
}

init_func InitObj::Init_Func(const string& init_func_name) const {
	map<string, init_func>::const_iterator it = init_func_map_.find(init_func_name);

	if (it == init_func_map_.end()){
		cout << "The requested init_func does not exist." << endl;
		Print();
		cout << "Requested function: " << init_func_name << endl;
		abort();
	}
	else return it -> second;
}

init_func_C InitObj::Init_Func_C(const string& init_func_name) const {
	map<string, init_func_C>::const_iterator it = init_func_C_map_.find(init_func_name);

	if (it == init_func_C_map_.end()){
		cout << "The requested init_func does not exist." << endl;
		Print_C();
		cout << "Requested function: " << init_func_name << endl;
		abort();
	}
	else return it -> second;
}

void InitObj::Print() const {
	map<string, init_func>::const_iterator it;

	cout << "Initial State Construction Functions: "<< endl;
	for (it = init_func_map_.begin(); it != init_func_map_.end(); it++){
		cout << it -> first << endl;
	}
}

void InitObj::Print_C() const {
	map<string, init_func_C>::const_iterator it;

	cout << "Initial State Construction Functions: "<< endl;
	for (it = init_func_C_map_.begin(); it != init_func_C_map_.end(); it++){
		cout << it -> first << endl;
	}
}

void InitObj::map_init_(){
	map<string, init_func>::const_iterator it;
	map<string, init_func_C>::const_iterator it_C;

	string name1 = "Random Product";
	init_func func1 = random_product;
	init_func_C func_C1 = random_product;

	it = init_func_map_.find(name1);
	if (it != init_func_map_.end()){
		cout << "init_func " << name1 << " already exists." << endl;
		abort();
	}
	init_func_map_[name1] = func1;

	it_C = init_func_C_map_.find(name1);
	if (it_C != init_func_C_map_.end()){
		cout << "init_func " << name1 << " already exists." << endl;
		abort();
	}
	init_func_C_map_[name1] = func_C1;

	string name2 = "Product Random";
	init_func func2 = product_random;

	it = init_func_map_.find(name2);
	if (it != init_func_map_.end()){
		cout << "init_func " << name2 << " already exists." << endl;
		abort();
	}
	init_func_map_[name2] = func2;

	string name3 = "Random Pure";
	init_func_C func_C3 = random_pure;

	it_C = init_func_C_map_.find(name3);
	if (it_C != init_func_C_map_.end()){
		cout << "init_func " << name3 << " already exists." << endl;
		abort();
	}
	init_func_C_map_[name3] = func_C3;

	string name4 = "Largest Leftmost Spin Z Value";
	init_func_C func_C4 = largest_leftmost_spin_z_complex_eigenstate;

	it_C = init_func_C_map_.find(name4);
	if (it_C != init_func_C_map_.end()){
		cout << "init_func " << name4 << " already exists." << endl;
		abort();
	}
	init_func_C_map_[name4] = func_C4;

	string name5 = "Leftmost Spin Z Value";
	init_func_C func_C5 = leftmost_spin_z_complex_eigenstate;

	it_C = init_func_C_map_.find(name5);
	if (it_C != init_func_C_map_.end()){
		cout << "init_func " << name5 << " already exists." << endl;
		abort();
	}
	init_func_C_map_[name5] = func_C5;

	string name6 = "Full Leftmost Spin Z Value";
	init_func_C func_C6 = leftmost_spin_z_complex_eigenstate;

	it_C = init_func_C_map_.find(name6);
	if (it_C != init_func_C_map_.end()){
		cout << "init_func " << name6 << " already exists." << endl;
		abort();
	}
	init_func_C_map_[name6] = func_C6;

}

void InitObj::Multi_Num_Init(const string& init_name, InitInfo& init_info) const {
	if (init_name == "Leftmost Spin Z Value"){
		init_info.multi_ini_para_num = init_info.multi_ini_para.leftmost_spin_z_index_set.size();
	}
	else if (init_name == "Full Leftmost Spin Z Value"){
		// Here the leftmost_spin_z_index_set is also initialized
		init_info.multi_ini_para_num = init_info.dim;

		init_info.multi_ini_para.leftmost_spin_z_index_set.resize(init_info.dim);
		for (int i=0; i<init_info.dim; i++) init_info.multi_ini_para.leftmost_spin_z_index_set[i] = i+1;
	}
	else{
		cout << init_name << " is not implemented with multiple sets of initial parameters." 
			 << endl;
		init_info.multi_ini_para_num = 0;
	}
}

string InitObj::Init_Para_String(const string& init_name, const InitInfo& init_info) const {
	stringstream output;
	if (init_name == "Leftmost Spin Z Value"){
		output << ",Leftmost_spin_z_index=" << init_info.leftmost_spin_z_index;
		return output.str();
	}
	else return "";
}


string InitObj::Init_Para_String(const string& init_name, const AllPara& para) const {
	stringstream output;
	if (init_name == "Leftmost Spin Z Value"){
		output << "Spin_z_index=" << para.evolution.leftmost_spin_z_index;
		return output.str();
	}
	else return "";
}
