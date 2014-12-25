#ifndef RESULTS_H
#define RESULTS_H

#include <vector>

using namespace std;

/**
This file defines a pure base class for processing raw data and output it to a file. The
raw data, which may come from multiple instances of calculations, should be passed in as
a vector. The object can also redirect the processed data to outside. The object can only
process data once and only store one set of data. If it is used to process multiple data sets,
then it must be reset in between, which loses all information about the previous data set.
**/

template <class T1, class T2>
class ResultsOutput
{
	public:
		/* The function which processes all the data. It may call different functions specific
		to particular derived classes to hand different data, and different bools to control
		excution of calculations. */
		virtual void Data_Process(const vector<T1>& ) = 0;

		// Output data to a file. The bool decides whether filenames are outputed. The
		// integer determines the width
		virtual void Data_Output(bool, int) const = 0; 
		virtual void Data_Redirect(T2&) const = 0; // Redirect data

		// Reset the object so that it can be used to process another set of data
		virtual void Reset() = 0;

		// Check whether it has processed one set of data, so that it cannot process another set
		virtual bool Empty() const = 0; 

		virtual ~ResultsOutput() {};
};

#endif
