#ifndef RESULTS_H
#define RESULTS_H

/**
This file defines a pure base class for processing raw data and output it to a file.
**/

template <class T>
class ResultsOutput
{
	public:
		/* The function which processes all the data. It may call different functions specific
		to particular derived classes to hand different data, and different bools to control
		excution of calculations. */
		virtual void Data_Process(const T&) = 0;

		virtual void Data_Output() = 0; // Output data to a file

		virtual ~ResultsOutput() {};
};

#endif
