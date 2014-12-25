#include <iomanip>
#include <iostream>
#include "output_func.h"

using namespace std;

ofstream Of_Construct(const stringstream& base, const string& post_string, bool output){
	stringstream filename;
	filename << base.str() <<"post_string";

	if (output) cout << filename.str() <<endl;

	ofstream fout(filename.str().c_str());

	return fout;
}

ofstream Of_Construct(const string& pre_string, const stringstream& base, bool output){
	stringstream filename;
	filename << pre_string << base.str();

	if (output) cout << filename.str() <<endl;

	ofstream fout(filename.str().c_str());

	return fout;
}

ofstream Of_Construct(const string& pre_string, const stringstream& base, 
const string& post_string, bool output){

	stringstream filename;
	filename << pre_string << base.str() << post_string;

	if (output) cout << filename.str() <<endl;

	ofstream fout(filename.str().c_str());

	return fout;
}

void Write_File(ofstream& fout, const vector<double>& data, int width){
	for (int i=0; i<data.size(); i++){
		fout << setw(width) << data[i];
	}
	fout<<endl;
}

void Write_File(ofstream& fout, int para, const vector<double>& data, int width){
	fout << setw(width) << para;
	for (int i=0; i<data.size(); i++){
		fout << setw(width) << data[i];
	}
	fout<<endl;
}