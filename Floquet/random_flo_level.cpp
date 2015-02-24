#include <ctime>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <Eigen/Core>
#include <utility>
#include <sstream>
#include "mtrand.h"
#include "flo_evol_model.h"
#include "level_stats.h"
#include "output_func.h"

using namespace std;



int main(){
	const int size = 2; // System Size
	const double tau = 0.8; // Time step, which seems not used here
	const int threads_N = 4; // Number of threads
	const int num_realization = 400; // Number of realizations.
	const int J_N = 10; // Number of points for J, with min=0, max=1.

	stringstream base_filename; // Filename

	bool output_init = false; // Whether output objects have been constructed

	ResultsOutput< EvolMatrix< ComplexEigenSolver<MatrixXcd> >*, 
				  pair<vector<double>, vector<double> > >* 
		result = new VanillaFloLevel(); // Output Object

	// The vector that holds processed data. In the pair, the first vector holds
	// level spacing, and the second holds: mean, mean sd, square mean, square mean sd.
	vector<pair<vector<double>, vector<double> > > data(J_N); 

	Eigen::initParallel();

//	vector<double> temp(2); // Used to hold data for output

	for (int i=0; i<J_N; i++){
		double J = (i+1) * (1.0/J_N);

		cout << J << endl;

		vector<EvolMatrix<ComplexEigenSolver<MatrixXcd> >* > floquet(num_realization);

		cout << "Initialize Evolution Operators." <<endl;
		for (int k=0; k<num_realization; k++){
			floquet[k] = new FloEvolRandom(size, tau, J);
		}

		if (!output_init){
			base_filename << floquet[0] -> Type()<<"_L=" << size << ",tau=" << tau;

			string post_string =  ",level_spacing.txt";
			ofstream level_out;
			Of_Construct(level_out, base_filename, post_string, true) ;

			post_string = ",level_spacing_mean.txt";
			ofstream mean_out;
			Of_Construct(mean_out, base_filename, post_string, true) ;

			post_string = ",level_spacing_square_mean.txt";
			ofstream square_mean_out;
			Of_Construct(square_mean_out, base_filename, post_string, true) ;

			output_init = true;
		}

/*		for (int k=0; k<num_realization; k++){
			floquet[k] -> Evol_Para_Init();
		}

		cout << "Diagonolize Evolution Operators." <<endl;

		#pragma omp parallel num_threads(threads_N)
		{
			#pragma omp for
			for (int k=0; k<num_realization;k++){
				floquet[k] -> Evol_Construct();
				floquet[k] -> Evol_Diag(false);
			}
		}

		cout << "Process Eigenvalues." <<endl;

		if (!result -> Empty()) result -> Reset();

		result -> Data_Process(floquet); 

		result -> Data_Redirect(data[i]);*/

		cout << "Output Eigenvalues." <<endl;

		Write_File(level_out, J, data[i].first, 15);

		temp[0] = data[i].second[0];
		temp[1] = data[i].second[1];
		Write_File(mean_out, J, temp, 15);

		temp[0] = data[i].second[2];
		temp[1] = data[i].second[3];
		Write_File(square_mean_out, J, temp, 15);
	}

	delete result;

	return 0;
}