#include <algorithm>
#include <vector>
#include <utility>

using namespace std;

/**
 ** This implements variaty of comparison functions used in sorting functions
 **/

bool Vec_Double_Sort(double i, double j) {return (i<j);}

bool Vec_Pair_Double_Int_First_Sort(pair<double, int> i, pair<double, int> j){
	return (i.first < j.first);
}