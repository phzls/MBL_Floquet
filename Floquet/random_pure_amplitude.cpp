#include <complex>
#include <vector>
#include "constants.h"
#include "mtrand.h"

using namespace std;

/**
 ** This function creates a vector of random complex amplitudes which can be used for random
 ** pure state
 **/

extern MTRand u1rand;

void random_pure_amplitude(vector<complex<double> >& amplitude){
	double sum = 0;
	for (int i=0;i<amplitude.size();i++)
	{
		double U1 = u1rand();
		double U2 = u1rand();

		double real = sqrt(-2*log(U1))*cos(2*Pi*U2);
		double imag = sqrt(-2*log(U1))*sin(2*Pi*U2);

		amplitude[i] = complex<double>(real,imag);

		sum += norm(amplitude[i]);
	}

	for (int i=0;i<amplitude.size();i++)
	{
		amplitude[i] /= complex<double>(sqrt(sum),0);
	}
}