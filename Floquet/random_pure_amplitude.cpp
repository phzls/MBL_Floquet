#include <complex>
#include <vector>
#include "constants.h"
#include "randomc.h"

/**
 ** This function creates a vector of random complex amplitudes which can be used for random
 ** pure state
 **/

extern CRandomMersenne RanGen_mersenne; // points in [0,1)

void random_pure_amplitude(vector<complex<double> >& amplitude){
	double sum = 0;
	for (int i=0;i<amplitude.size();i++)
	{
		double U1 = RanGen_mersenne.Random();
		double U2 = RanGen_mersenne.Random();

		double real = sqrt(-2*log(1-U1))*cos(2*Pi*U2);
		double imag = sqrt(-2*log(1-U1))*sin(2*Pi*U2);

		amplitude[i] = complex<double>(real,imag);

		sum += norm(amplitude[i]);
	}

	for (int i=0;i<amplitude.size();i++)
	{
		amplitude[i] /= complex<double>(sqrt(sum),0);
	}
}