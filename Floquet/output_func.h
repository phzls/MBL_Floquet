#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

/**
 ** This file contains some functions which are used for writing out files.
 **/

/* 
 * Given a base filename and string after it. The bool determines whether the filename 
 * is outputted
 */
void Of_Construct(ofstream&, const stringstream&, const string&, bool);

/* 
 * Given a base filename and string before it. The bool determines whether the filename 
 * is outputted
 */
void Of_Construct(ofstream&, const string&, const stringstream&, bool);

/* 
 * Given a base filename and strings both before and after it. The bool determines whether  
 * the filename is outputted
 */
void Of_Construct(ofstream&, const string&, const stringstream&, const string&, bool);

/* 
 * Write out the vector in a line to a file. The integer gives the space between two numbers
 * in the output files.
 */
void Write_File(ofstream&, const vector<double>&, int);

/* 
 * Write out the vector in a line to a file. The double is some parameter, and the
 * integer gives the space between two numbers in the output files.
 */
void Write_File(ofstream&, double, const vector<double>&, int);
