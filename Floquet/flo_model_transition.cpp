//
// Created by Liangsheng Zhang on 3/24/15.
//

#include <iostream>
#include <cstdlib>
#include "flo_model_transition.h"

using namespace std;

void FloModelTransition::map_initialize_(const AllPara& parameters) {

    map<string, Flo_init >::iterator init_it;
    map<string, Flo_func >::iterator cal_it;
    map<string, Flo_out >::iterator out_it;

    // end-to-end z-z correlation square
    string name1 = "ZZ Correlation Square";
    Flo_init init_func1 = &FloModelTransition::ZZ_corr_square_init_;
    Flo_func cal_func1 = &FloModelTransition::ZZ_corr_square_compute_;
    Flo_out out_func1 = &FloModelTransition::ZZ_corr_square_out_;

    // Make sure the name has not been used before
    init_it = flo_init_map_.find(name1);
    cal_it = flo_func_map_.find(name1);
    out_it = flo_out_map_.find(name1);
    if (init_it != flo_init_map_.end() || cal_it != flo_func_map_.end() || out_it != flo_out_map_.end()){
        cout << name1 << " for floquet transition has appeared before." << endl;
        abort();
    }

    flo_init_map_[name1] = init_func1;
    flo_func_map_[name1] = cal_func1;
    flo_out_map_[name1] = out_func1;

    // entropy variance
    string name2 = "Entropy Variance";
    Flo_init init_func2 = &FloModelTransition::Ent_var_init_;
    Flo_func cal_func2 = &FloModelTransition::Ent_var_compute_;
    Flo_out out_func2 = &FloModelTransition::Ent_var_out_;

    // Make sure the name has not been used before
    init_it = flo_init_map_.find(name2);
    cal_it = flo_func_map_.find(name2);
    out_it = flo_out_map_.find(name2);
    if (init_it != flo_init_map_.end() || cal_it != flo_func_map_.end() || out_it != flo_out_map_.end()){
        cout << name2 << " for floquet transition has appeared before." << endl;
        abort();
    }

    flo_init_map_[name2] = init_func2;
    flo_func_map_[name2] = cal_func2;
    flo_out_map_[name2] = out_func2;

    // entropy variance
    string name3 = "Entropy Variance Smallest";
    Flo_init init_func3 = &FloModelTransition::Ent_smallest_var_init_;
    Flo_func cal_func3 = &FloModelTransition::Ent_smallest_var_compute_;
    Flo_out out_func3 = &FloModelTransition::Ent_smallest_var_out_;

    // Make sure the name has not been used before
    init_it = flo_init_map_.find(name3);
    cal_it = flo_func_map_.find(name3);
    out_it = flo_out_map_.find(name3);
    if (init_it != flo_init_map_.end() || cal_it != flo_func_map_.end() || out_it != flo_out_map_.end()){
        cout << name3 << " for floquet transition has appeared before." << endl;
        abort();
    }

    flo_init_map_[name3] = init_func3;
    flo_func_map_[name3] = cal_func3;
    flo_out_map_[name3] = out_func3;

    // Check the number of function
    if ( flo_init_map_.size() != flo_func_map_.size() ){

        cout << "Number of initializing functions in flo_transition is not the same as number of registered functions"
                << "that can be called." << endl;
        cout << "Functions in parameters:" << endl;
        for (map<string, Flo_init >::iterator it = flo_init_map_.begin();
             it != flo_init_map_.end(); it++){
            cout << it -> first << endl;
        }
        cout << "Total Number: " << flo_init_map_.size() << endl;

        cout << "Functions that can be calledt:" << endl;
        for (map<string, Flo_func>::iterator it = flo_func_map_.begin();
             it != flo_func_map_.end(); it ++){
            cout << it -> first << endl;
        }
        cout << "Total Number:" << flo_func_map_.size() << endl;

        abort();
    }

    if ( flo_out_map_.size() != flo_func_map_.size() ){

        cout << "Number of output functions in flo_transition is not the same as number of registered functions"
                << "that can be called." << endl;
        cout << "Functions in parameters:" << endl;
        for (map<string, Flo_out>::iterator it = flo_out_map_.begin();
             it != flo_out_map_.end(); it++){
            cout << it -> first << endl;
        }
        cout << "Total Number: " << flo_out_map_.size() << endl;

        cout << "Functions that can be calledt:" << endl;
        for (map<string, Flo_func>::iterator it = flo_func_map_.begin();
             it != flo_func_map_.end(); it ++){
            cout << it -> first << endl;
        }
        cout << "Total Number:" << flo_func_map_.size() << endl;

        abort();
    }

    if ( flo_func_bool_map_.size() > flo_func_map_.size() ){

        cout << "Number of registered functions in flo_transition is larger than the number of registered functions"
                << "that can be called." << endl;
        cout << "Functions in parameters:" << endl;
        for (map<string, bool>::iterator it = flo_func_bool_map_.begin();
             it != flo_func_bool_map_.end(); it++){
            cout << it -> first << endl;
        }
        cout << "Total Number: " << flo_func_bool_map_.size() << endl;

        cout << "Functions that can be calledt:" << endl;
        for (map<string, Flo_func>::iterator it = flo_func_map_.begin();
             it != flo_func_map_.end(); it ++){
            cout << it -> first << endl;
        }
        cout << "Total Number:" << flo_func_map_.size() << endl;

        abort();
    }

    // Initialize data structures
    for (map<string, bool>::iterator it = flo_func_bool_map_.begin(); it != flo_func_bool_map_.end(); it++){
        if (it -> second) ( this ->*(flo_init_map_[it -> first]) ) (parameters);
    }
}

void FloModelTransition::Compute(AllPara const & parameters,
        EvolMatrix <ComplexEigenSolver<MatrixXcd> > const * floquet,
        LocalInfo const & local_info) {

    for (map<string, bool>::iterator it = flo_func_bool_map_.begin(); it != flo_func_bool_map_.end(); it++){
        if (it -> second) ( this ->*(flo_func_map_[it -> first]) )(parameters, floquet, local_info);
    }
}

void FloModelTransition::Output(AllPara const & parameters, const string& name) {
    for (map<string, bool>::iterator it = flo_func_bool_map_.begin(); it != flo_func_bool_map_.end(); it++){
        if (it -> second) ( this ->*(flo_out_map_[it -> first]) )(parameters, name);
    }
}

