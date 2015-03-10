#include "evol_data_total.h"

using namespace std;

/**
 ** This file contains some common functions in evol_data.h
 **/

EvolDataTotal::EvolDataTotal(const AllPara& parameters): size_(parameters.generic.size){
    func_status_ = parameters.evolution.evol_total_compute;
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
        for (map<string, Data_Total_Init>::iterator it = data_init_.begin(); it != data_init_.end();
             it ++){
            cout << it -> first << endl;
        }
        cout << "Total Number:" << data_init_.size() << endl;

        abort();
    }

    Name_Check_();

    // Initialize all quantities which will be computed
    for (map<string, bool>::iterator it = func_status_.begin(); it != func_status_.end(); it++){
        if (it -> second) ( this ->* (data_init_[it -> first]) ) (parameters);
    }
}

void EvolDataTotal::Print_All_Name() const{
    map<string, bool>::const_iterator it;

    for (it = func_status_.begin(); it != func_status_.end(); it ++){
        cout << it -> first << endl;
    }
}

void EvolDataTotal::Print_All_Status() const{
    map<string, bool>::const_iterator it;

    for (it = func_status_.begin(); it != func_status_.end(); it ++){
        cout << it -> first << "   ";

        if (it -> second) cout << "Would be computed." << endl;
        else cout << "Would not be computed." << endl;
    }
}

void EvolDataTotal::Name_Check_() const{
    map<string, bool>::const_iterator para_it;
    map<string, Data_Total_Init>::const_iterator data_it;

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

void EvolDataTotal::Data_Compute(const VectorXcd& state, const StepInfo& info){
    for (map<string, bool>::iterator it = func_status_.begin(); it != func_status_.end(); it++){
        if (it -> second) ( this ->* (data_cal_[it -> first]) ) (state, info);
    }
}

void EvolDataTotal::Data_Compute(const MatrixXcd& state_density, const StepInfo& info){
    for (map<string, bool>::iterator it = func_status_.begin(); it != func_status_.end(); it++){
        if (it -> second) ( this ->* (data_cal_C_[it -> first]) ) (state_density, info);
    }
}

void EvolDataTotal::Data_Total_Output(const AllPara& parameters, const string& type_name){
    for (map<string, bool>::iterator it = func_status_.begin(); it != func_status_.end(); it++){
        if (it -> second){
            map<string, Data_Total_Out>::iterator out_it = data_out_total_.find(it -> first);
            if (out_it != data_out_total_.end())
                ( this ->* (data_out_total_[it -> first]) ) (parameters, type_name);
        }
    }
}

void EvolDataTotal::Data_Func_Map_Init_(){
    map<string, Data_Total_Init>::iterator init_it;
    map<string, Data_Total_Cal>::iterator cal_it;
    map<string, Data_Total_Cal_C>::iterator cal_C_it;
    map<string, Data_Total_Out>::iterator out_total_it;

    // Leftmost Spin Z One Run data
    string name1 = "Leftmost Spin Z One Run";
    Data_Total_Init init_func1 = &EvolDataTotal::Leftmost_Spin_Z_One_Run_Init_;
    Data_Total_Cal cal_func1 = &EvolDataTotal::Leftmost_Spin_Z_One_Run_Cal_;
    Data_Total_Cal_C cal_C_func1 = &EvolDataTotal::Leftmost_Spin_Z_One_Run_Cal_C_;
    Data_Total_Out out_total_func1 = &EvolDataTotal::Leftmost_Spin_Z_One_Run_Out_;

    // Make sure the name has not been used before
    init_it = data_init_.find(name1);
    cal_it = data_cal_.find(name1);
    out_total_it = data_out_total_.find(name1);
    if (init_it != data_init_.end() || cal_it != data_cal_.end() ||
            out_total_it != data_out_total_.end()){
        cout << name1 << " for evolution has appeared before." << endl;
        abort();
    }

    data_init_[name1] = init_func1;
    data_cal_[name1] = cal_func1;
    data_cal_C_[name1] = cal_C_func1;
    data_out_total_[name1] = out_total_func1;

    // Full Leftmost Spin Z Per Model data
    string name2 = "Full Leftmost Spin Z Per Model";
    Data_Total_Init init_func2 = &EvolDataTotal::Full_Leftmost_Spin_Z_Per_Model_Init_;
    Data_Total_Cal cal_func2 = &EvolDataTotal::Full_Leftmost_Spin_Z_Per_Model_Cal_;
    Data_Total_Cal_C cal_C_func2 = &EvolDataTotal::Full_Leftmost_Spin_Z_Per_Model_Cal_C_;
    Data_Total_Out out_total_func2 = &EvolDataTotal::Full_Leftmost_Spin_Z_Per_Model_Out_;

    // Make sure the name has not been used before
    init_it = data_init_.find(name2);
    cal_it = data_cal_.find(name2);
    out_total_it = data_out_total_.find(name2);
    if (init_it != data_init_.end() || cal_it != data_cal_.end() ||
            out_total_it != data_out_total_.end()){
        cout << name2 << " for evolution has appeared before." << endl;
        abort();
    }

    data_init_[name2] = init_func2;
    data_cal_[name2] = cal_func2;
    data_cal_C_[name2] = cal_C_func2;
    data_out_total_[name2] = out_total_func2;

    // Check data_init_ and data_cal_ have the same size
    if (data_init_.size() != data_cal_.size()){
        cout << "Number of initializations in evolution is not the same as number of calculations."
                << endl;
        cout << "Registered initializations:" << endl;
        for (init_it = data_init_.begin(); init_it != data_init_.end(); init_it ++){
            cout << init_it -> first << endl;
        }
        cout << "Total Number: " << data_init_.size();

        cout << "Registered calculations:" << endl;
        for (cal_it = data_cal_.begin(); cal_it != data_cal_.end(); cal_it ++){
            cout << cal_it -> first << endl;
        }
        cout << "Total Number: " << data_cal_.size();
    }

    // Check data_init_ and data_out_ have the same size
    if (data_init_.size() !=  data_out_total_.size()){
        cout << "Number of initializations in evolution is not the same as number of output."
                << endl;
        cout << "Registered initializations:" << endl;
        for (init_it = data_init_.begin(); init_it != data_init_.end(); init_it ++){
            cout << init_it -> first << endl;
        }
        cout << "Total Number: " << data_init_.size();

        cout << "Registered output for total:" << endl;
        for (out_total_it = data_out_total_.begin(); out_total_it != data_out_total_.end();
             out_total_it ++){
            cout << out_total_it -> first << endl;
        }
        cout << "Total Number: " << data_out_total_.size();
    }

}