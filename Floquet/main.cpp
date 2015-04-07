#include <iostream>
#include <ctime>
#include <string>
#include "constants.h"
#include "parameters.h"
#include "tasks_models.h"
#include "randomc.h"

using namespace std;

CRandomMersenne RanGen_mersenne(time(NULL));

TasksModels tasks_models; // Record all the tasks and methods

int main() {
    AllPara parameters;




//===================================================================================





    parameters.generic.task = "Flo Transition";
    parameters.generic.model = "XXZ Random Simple Flo";

    parameters.generic.size = 6; // System size
    parameters.generic.num_realizations = 100; // Number of realizations
    parameters.generic.threads_N = 4; // Number of threads in openmp
    parameters.generic.evec = false; // Whether compute eigenvectors, so far only called in
    // level statistics calculation
    parameters.generic.erase = true; // Whether erase matrix after diagonization
    parameters.generic.debug = false; // Whether output debug information
    parameters.generic.iso_keep = true; // Whether isolated part is kept
    parameters.generic.version = 1; // Version of the output
    parameters.generic.time = false; // Whether the program is timed

    parameters.output.width = 30; // Width for spacing in output files
    parameters.output.filename_output = true; // Whether print out file names

    parameters.floquet.J_N = 1; // Number of points of coupling strength
    parameters.floquet.J_min = 0.6; // Minimum J
    parameters.floquet.J_max = 0.9; // Maximum J
    parameters.floquet.tau = 0.8; // Time step size
    parameters.floquet.J = 0.6;

    parameters.markov.K = 0.8; // Coupling strength to the bath in markov models

    parameters.floquet_random.angle_min = 5 * Pi / 6; // Minimum angle for the random rotation angle
    parameters.floquet_random.angle_sup = 5 * Pi / 6; // Supreme angle for the random rotation angle

    parameters.floquet_xxz.g = 0.9045; // Transverse field strength
    parameters.floquet_xxz.h = 0.8090; // Longitude field strength




//===============================================================================




    parameters.evolution.time_step = 20; // Number of time steps
    parameters.evolution.step_size = parameters.floquet.tau; // Time step size
    parameters.evolution.init_func_name = "Leftmost Spin Random State"; // Initial state name

    parameters.evolution.evol_compute["Entropy Per Model"] = false;
    parameters.evolution.evol_compute["Leftmost Spin Z Per Model"] = true;
    parameters.evolution.evol_compute["Leftmost Spin X Per Model"] = true;
    parameters.evolution.evol_compute["Leftmost Spin Y Per Model"] = true;

    parameters.evolution.evol_total_compute["Leftmost Spin Z One Run"] = false;
    parameters.evolution.evol_total_compute["Full Leftmost Spin Z Per Model"] = false;

    parameters.evolution.model_num = 1; // Number of models for evolution
    // If partition the chain to two halves, the size of left part
    parameters.evolution.left_size = parameters.generic.size / 2;
    parameters.evolution.jump = 1; // jump of time points in evolution

    parameters.evolution.markov_time_jump = 10; // Time jump in the markov process for flipping
    // True time: time_step * markov_time_jump
    parameters.evolution.markov_jump = true; // Determine whether there will be markov_time_jump
    if (parameters.generic.task.find("Markov") == string::npos &&
        parameters.generic.task.find("Single Model") == string::npos) {
        // The task has nothing to do with markov process
        parameters.evolution.markov_jump = false;
    }

    parameters.evolution.log_time = false; // whehter time changes logarithmically
    parameters.evolution.log_time_jump = 2; // The base for time change logarithmically

    // The number gives the index of leftmost spin z value
    parameters.evolution.leftmost_spin_z_index = 3;

    // Multiple sets of initial conditions
    int leftmost_spin_z_index_set[] = {1, 2, 3, 4, 5, 63, 64};
    parameters.multi_ini_para.leftmost_spin_z_index_set.assign(leftmost_spin_z_index_set,
                                                               leftmost_spin_z_index_set +
                                                               sizeof(leftmost_spin_z_index_set) / sizeof(int));

    // Whether output the evolution of leftmot spin z for all eigenstates in full_leftmost_spin_z
    parameters.multi_ini_para.easy_full_leftmost_spin_z = false;

    // Threshold in time evolution of leftmost spin z for non-zero values
    parameters.multi_ini_para.non_zero_threshold = 0.001;




//=================================================================================





    // Methods to be called under single model
    parameters.single_model.single_model_compute["Flo Chain End Sigma Z"] = false;
    parameters.single_model.single_model_compute["Flo Evolution Simple Markov"] = true;




//===================================================================================




    // Methods to be called for studying transition of floquet systems from thermal to localization
    parameters.transition.flo_transition_compute["ZZ Correlation Square"] = false; // End-to-end sigma_z X sigma_z
    // correlation square

    parameters.transition.flo_transition_compute["Entropy Variance"] = false; // Entropy variance for all eigenstates

    parameters.transition.flo_transition_compute["Entropy Variance Smallest"] = false; // Entropy variance eigenstates with
    // smallest phase magnitude among all
    // realizations

    parameters.transition.flo_transition_compute["ZZ Time Correlation"] = false; // End-to-end sigma_z X sigma_z
    // time correlation

    parameters.transition.flo_transition_compute["ZZ Time Correlation Components"] = false; // End-to-end
    // sigma_z X sigma_z time correlation components. The first row in output is text header

    parameters.transition.flo_transition_compute["ZZ All Correlation Square"] = true; // zz correlation square
                                                                                      // at all distances



//============================================================================================



    // Methods to be called for studying eigenstate properties of floquet systems
    parameters.eigenvec.flo_eigen_compute["ZZ Correlation Square"] = false; // End-to-end sigma_z sigma_z correlation
                                                                            // square
    parameters.eigenvec.flo_eigen_compute["Eigenvectors"] = true; // End-to-end sigma_z sigma_z correlation
                                                                            // square



//============================================================================================




    Eigen::initParallel();

    tasks_models.Task(parameters.generic.task)(parameters);

    cout << "Calculation finished." << endl;

    return 0;
}