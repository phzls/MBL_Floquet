#include <ctime>
#include <string>
#include <cstdlib>
#include <iostream>
#include <Eigen/Core>
#include <utility>
#include <sstream>
#include "flo_evol_model.h"
#include "level_stats.h"
#include "output_func.h"
#include "methods.h"
#include "parameters.h"
#include "tasks_models.h"

using namespace std;

/*
 * Compute level spacing for a Floquet operator.
 */

extern TasksModels tasks_models; // Record all the tasks and methods. Defined in main.

void flo_level(const AllPara& parameters){

	// System Size
	const int size = parameters.generic.size; 
	// Number of realizations.
	const int num_realization = parameters.generic.num_realizations;
	const string model = parameters.generic.model;

	const double tau = parameters.floquet.tau; // Time step, which seems not used here
	const int J_N = parameters.floquet.J_N; // Number of points for J
	const double J_max = parameters.floquet.J_max; // Maximum J
	const double J_min = parameters.floquet.J_min; // Minimum J

	const double angle_min = parameters.floquet_random.angle_min; // Minimum angle
	const double angle_sup = parameters.floquet_random.angle_sup; // Supreme angle	

	const int width = parameters.output.width; // Output spacing

	AllPara local_parameters(parameters); // Local parameters which can be changed

	stringstream base_filename; // Filename

	ofstream level_out; // Output for level spacings
	ofstream mean_out; // Output for mean of level spacings
	ofstream square_mean_out; // Output for square mean of level spacings

	bool output_init = false; // Whether output objects have been constructed

	ResultsOutput< EvolMatrix< ComplexEigenSolver<MatrixXcd> >*, 
				  pair<vector<double>, vector<double> > >* 
		result = new VanillaFloLevel(); // Output Object

	// The vector that holds processed data. In the pair, the first vector holds
	// level spacing, and the second holds: mean, mean sd, square mean, square mean sd.
	vector<pair<vector<double>, vector<double> > > data(J_N); 

	vector<double> temp(2); // Used to hold data for output

	Eigen::initParallel();

	for (int i=0; i<J_N; i++){
		double J;
		if (J_N > 1) J = J_min + i * (J_max - J_min)/(J_N-1);
		else J = J_min;

		local_parameters.floquet.J = J;

		cout << local_parameters.floquet.J << endl;

		vector<EvolMatrix<ComplexEigenSolver<MatrixXcd> >* > floquet(num_realization);

		cout << "Initialize Evolution Operators." <<endl;
		for (int k=0; k<num_realization; k++){
			//floquet[k] = new FloEvolRandom(size, tau, J);
			//floquet[k] = new FloEvolRandomRotation(size, tau, J, angle_min, angle_sup);
			tasks_models.Model(model, local_parameters, floquet[k]);
		}

		if (!output_init){
			base_filename << floquet[0] -> Type()<<"_L=" << size << ",tau=" << tau
						  << ",Realizations=" << num_realization << ",J_min="
						  << J_min << ",J_max=" << J_max << ",J_N=" << J_N;

			int found = floquet[0] -> Type().find("Rotation");
			if (found != string::npos){
				// It is the random rotation floquet
				base_filename << ",angle_min=" << angle_min <<",angle_sup=" << angle_sup;
			}

			string post_string =  ",level_spacing.txt";
			Of_Construct(level_out, base_filename, post_string, true) ;

			post_string = ",level_spacing_mean.txt";
			Of_Construct(mean_out, base_filename, post_string, true) ;

			post_string = ",level_spacing_square_mean.txt";
			Of_Construct(square_mean_out, base_filename, post_string, true) ;

			output_init = true;
		}

		level_cal(parameters, floquet, result, data[i]);

		cout << "Output Eigenvalues." <<endl;

		Write_File(level_out, J, data[i].first, width);

		temp[0] = data[i].second[0];
		temp[1] = data[i].second[1];
		Write_File(mean_out, J, temp, width);

		temp[0] = data[i].second[2];
		temp[1] = data[i].second[3];
		Write_File(square_mean_out, J, temp, width);

		for (int k=0; k<num_realization;k++){
			delete floquet[k];
			floquet[k] = NULL;
		}
	}

	delete result;
	result = NULL;
}