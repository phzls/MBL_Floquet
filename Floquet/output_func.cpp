#include <iomanip>
#include <iostream>
#include "output_func.h"

using namespace std;

void Of_Construct(ofstream& fout, const stringstream& base, const string& post_string, bool output){
	stringstream filename;
	filename << base.str() <<post_string;

	if (output) cout << filename.str() <<endl;

	fout.open( filename.str().c_str() );
}

void Of_Construct(ofstream& fout, const string& pre_string, const stringstream& base, bool output){
	stringstream filename;
	filename << pre_string << base.str();

	if (output) cout << filename.str() <<endl;

	fout.open( filename.str().c_str() );
}

void Of_Construct(ofstream& fout, const string& pre_string, const stringstream& base, 
const string& post_string, bool output){

	stringstream filename;
	filename << pre_string << base.str() << post_string;

	if (output) cout << filename.str() <<endl;

	fout.open( filename.str().c_str() );
}

void Write_File(ofstream& fout, const vector<double>& data, int width){
	for (int i=0; i<data.size(); i++){
		fout << setw(width) << data[i];
	}
	fout<<endl;
}

void Write_File(ofstream& fout, double para, const vector<double>& data, int width){
	fout << setw(width) << para;
	for (int i=0; i<data.size(); i++){
		fout << setw(width) << data[i];
	}
	fout<<endl;
}

void Write_File(ofstream& fout, const vector<complex<double> >& data, int width){
	for (int i=0; i<data.size(); i++){
		double re = real(data[i]);
		double im = imag(data[i]);

		fout << setw(width) << re;

		if (im < 0) fout << "-";
		else fout << "+";

		fout << im << "j";
	}
	fout<<endl;
}

void Write_File(ofstream& fout, double para, const vector<complex<double> >& data, int width){
	fout << setw(width) << para;
	for (int i=0; i<data.size(); i++){
		double re = real(data[i]);
		double im = imag(data[i]);

		fout << setw(width) << re;

		if (im < 0) fout << "-";
		else fout << "+";

		fout << im << "j";
	}
	fout<<endl;
}