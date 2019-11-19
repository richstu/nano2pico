#include "ttz_producer.hpp"

#include "utilities.hpp"
#include <algorithm>

using namespace std;

TTZVarProducer::TTZVarProducer(int year_){
    year = year_;
}

TTZVarProducer::~TTZVarProducer(){
}

void TTZVarProducer::WriteTTZVars(pico_tree &pico){
	//write signal_lepton_pt
	for (unsigned int el_idx = 0; el_idx < pico.out_el_pt().size(); el_idx++) {
		if (pico.out_el_sig()[el_idx]) {
			pico.out_signal_lepton_pt().push_back(pico.out_el_pt()[el_idx]);
		}
	}
	for (unsigned int mu_idx = 0; mu_idx < pico.out_mu_pt().size(); mu_idx++) {
		if (pico.out_mu_sig()[mu_idx]) {
			pico.out_signal_lepton_pt().push_back(pico.out_mu_pt()[mu_idx]);
		}
	}
	sort(pico.out_signal_lepton_pt().begin(), pico.out_signal_lepton_pt().end(), greater<float>());	
	return;
}

