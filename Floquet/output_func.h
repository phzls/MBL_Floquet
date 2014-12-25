#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

/**
This file contains some functions which are used for writing out files.
**/

// Given a base filename and string after it. The bool determines whether the filename 
// is outputted
ofstream Of_Construct(const stringstream&, const string&, bool);

// Given a base filename and string before it. The bool determines whether the filename 
// is outputted
ofstream Of_Construct(const string&, const stringstream&, bool);

// Given a base filename and strings both before and after it. The bool determines whether  
// the filename is outputted
ofstream Of_Construct(const string&, const stringstream&, const string&, bool);

// Write out the vector in a line to a file. The integer gives the space between two numbers
void Write_File(ofstream&, const vector<double>&, int);

// Write out the vector in a line to a file. The first integer is some parameter, and the
// last integer gives the space between two numbers
void Write_File(ofstream&, int, const vector<double>&, int);
