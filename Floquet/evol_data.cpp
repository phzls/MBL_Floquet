#include "evol_data.h"

using namespace std;

/**
 ** This file contains some common functions in evol_data.h
 **/

EvolData::EvolData(const AllPara& parameters): size_(parameters.generic.size){
	func_status_ = parameters.evolution.evol_compute;
	Data_Func_Map_Init_();

	if ( func_status_.size() != data_init_.size()){
		cout << "Number of registered functions in parameters for evolution is not the same as"
			 << " number of registered functions that can be called." << endl;
		cout << "Functions in parameters:" << endl;
		for (map<string, bool>::iterator it = func_status_.begin(); it != func_status_.end(); it++){
			cout << it -> first << endl;
		}
		cout << "Total Number: " << func_status_.size() << endl;

		cout << "Functions that can be called:" << endl;
		for (map<string, Data_Init>::iterator it = data_init_.begin(); it != data_init_.end();
			it ++){
			cout << it -> first << endl;
		}
		cout << "Total Number:" << data_init_.size() << endl;

		abort();
	}

	Name_Check_();

	for (map<string, bool>::iterator it = func_status_.begin(); it != func_status_.end(); it++){
		if (it -> second) ( this ->* (data_init_[it -> first]) ) (parameters);
	}
}

void EvolData::Print_All_Name() const{
	map<string, bool>::const_iterator it;

	for (it = func_status_.begin(); it != func_status_.end(); it ++){
		cout << it -> first << endl;
	}
}

void EvolData::Print_All_Status() const{
	map<string, bool>::const_iterator it;

	for (it = func_status_.begin(); it != func_status_.end(); it ++){
		cout << it -> first << "   ";

		if (it -> second) cout << "Would be computed." << endl;
		else cout << "Would not be computed." << endl;
	}
}

void EvolData::Name_Check_() const{
	map<string, bool>::const_iterator para_it;
	map<string, Data_Init>::const_iterator data_it;

	for (para_it = func_status_.begin(); para_it != func_status_.end(); para_it ++){
		data_it = data_init_.find(para_it -> first);
		if (data_it == data_init_.end()){
			cout << "Names in evolution are not consistent." << endl;

			cout << "Names in parameters:" << endl;
			map<string, bool>::const_iterator it;
			for (it = func_status_.begin(); it!= func_status_.end(); it++){
				cout << it -> first << endl;
			}

			cout << "Names in EvolData:" << endl;
			for (data_it = data_init_.begin(); data_it != data_init_.end(); data_it ++){
				cout << data_it -> first << endl;
			}
			abort();
		}
	}
}

void EvolData::Data_Func_Map_Init_(){
	map<string, Data_Init>::iterator init_it;
	map<string, Data_Cal>::iterator cal_it;

	string name1 = "Entropy Per Model";
	Data_Init init_func1 = &EvolData::Entropy_Per_Model_Init_;
	Data_Cal cal_func1 = &EvolData::Entropy_Per_Model_Cal_;

	init_it = data_init_.find(name1);
	cal_it = data_cal_.find(name1);
	if (init_it != data_init_.end() || cal_it != data_cal_.end()){
		cout << name1 << " for evolution has appeared before." << endl;
		abort();
	}

	data_init_[name1] = init_func1;
	data_cal_[name1] = cal_func1;

	// Check data_init_ and data_cal_ have the same size
	if (data_init_.size() != data_cal_.size()){
		cout << "Number of initializations in evolution is not the same as number of calculations."
			 << endl;
		cout << "Registered initializations:" << endl;
		for (init_it = data_init_.begin(); init_it != data_init_.end(); init_it ++){
			cout << init_it -> first << endl;
		}
		cout << "Total Number: " << data_init_.size();

		for (cal_it = data_cal_.begin(); cal_it != data_cal_.end(); cal_it ++){
			cout << cal_it -> first << endl;
		}
		cout << "Total Number: " << data_cal_.size(); 
	}

}