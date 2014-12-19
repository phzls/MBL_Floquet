#include "flo_evol_model.h"

using namespace std;

void RandomFloData::Initialize(){
	
}

void FloEvolRandom::Repr_Init_(){
	repr_ << "Random_Floquet_L=" << size_ << ",J=" << param_.J <<",tau_"
		  << param_.tau;
}

void FloEvolRandom::Evol_Construct(){

}