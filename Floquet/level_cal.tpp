#include <iostream>
#include <omp.h>
#include <utility>
#include "evol_class.h"
#include "results.h"
#include "output_func.h"
#include "parameters.h"

using namespace std;

template <class T1, class T2>
void level_cal(const AllPara& parameters, vector<EvolMatrix<T1>*>& floquet, 
	ResultsOutput<EvolMatrix<T1>*, T2>* result, T2& data){

	const int threads_N = parameters.generic.threads_N;
	const bool evec = parameters.generic.evec;
	const bool erase = parameters.generic.erase;

	for (int k=0; k<floquet.size(); k++){
		floquet[k] -> Evol_Para_Init();
	}

	cout << "Diagonalize Evolution Operators." <<endl;

	#pragma omp parallel num_threads(threads_N)
	{
		#pragma omp for
		for (int k=0; k<floquet.size();k++){
			floquet[k] -> Evol_Construct();
			floquet[k] -> Evol_Diag(evec);
			if (erase) floquet[k] -> Evol_Erase();
		}
	}

	cout << "Process Eigenvalues." <<endl;

	if (!result -> Empty()) result -> Reset();

	result -> Data_Process(floquet); 

	result -> Data_Redirect(data);
}


template <class T1, class T2>
void level_cal(const AllPara& parameters, vector<EvolMatrix<T1>*>& floquet, 
	ResultsOutput<EvolMatrix<T1>*, T2>* result){

	vector<double> temp(2); // Used to hold data for output

	const int threads_N = parameters.generic.threads_N;
	const bool evec = parameters.generic.evec;
	const bool erase = parameters.generic.erase;

	const int width = parameters.output.width;
	const bool filename_output = parameters.output.filename_output;


	for (int k=0; k<floquet.size(); k++){
		floquet[k] -> Evol_Para_Init();
	}

	cout << "Diagonalize Evolution Operators." <<endl;

	#pragma omp parallel num_threads(threads_N)
	{
		#pragma omp for
		for (int k=0; k<floquet.size();k++){
			floquet[k] -> Evol_Construct();
			floquet[k] -> Evol_Diag(evec);
			if (erase) floquet[k] -> Evol_Erase();
		}
	}

	cout << "Process Eigenvalues." <<endl;

	if (!result -> Empty()) result -> Reset();

	result -> Data_Process(floquet); 

	result -> Data_Output(filename_output, width);
}