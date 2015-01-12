#include <algorithm>
#include <vector>
#include <utility>

using namespace std;

/**
 ** This is the header file for a variaty of comparison functions used in sorting functions
 **/

bool Vec_Double_Sort(double, double); 

template <class T>
bool Vec_Pair_Double_First_Sort(pair<double, T> i, pair<double, T> j){
	return (i.first < j.first);
}