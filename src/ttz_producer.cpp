#include "ttz_producer.hpp"

#include "utilities.hpp"
#include "TLorentzVector.h"

using namespace std;

TTZVarProducer::TTZVarProducer(int year_){
    year = year_;
}

TTZVarProducer::~TTZVarProducer(){
}

float get_deltaphi(float phi1, float phi2) {
	//function to get Delta Phi between two particles
	return TMath::Min(TMath::Min(static_cast<Float_t>(TMath::Abs(phi2-phi1)),static_cast<Float_t>(TMath::Abs(phi2+2*3.1415-phi1))),static_cast<Float_t>(TMath::Abs(phi2-2*3.1415-phi1)));
}

float get_deltar(float eta1, float phi1, float eta2, float phi2) {
	//function that returns Delta R given eta and phi of particles 1 and 2
	return TMath::Sqrt(TMath::Power(eta2-eta1,2)+TMath::Power(TMath::Min(TMath::Min(static_cast<Float_t>(TMath::Abs(phi2-phi1)),static_cast<Float_t>(TMath::Abs(phi2+2*3.1415-phi1))),static_cast<Float_t>(TMath::Abs(phi2-2*3.1415-phi1))),2));
}

void TTZVarProducer::WriteTTZVars(pico_tree &pico){
	//calculate min_dilep_dr
	std::vector<float> signal_lepton_eta;
	std::vector<float> signal_lepton_phi;
	for (unsigned int el_idx = 0; el_idx < pico.out_el_sig().size(); el_idx++) {
		if (pico.out_el_sig()[el_idx]) {
			signal_lepton_eta.push_back(pico.out_el_eta()[el_idx]);
			signal_lepton_phi.push_back(pico.out_el_phi()[el_idx]);
		}
	}
	for (unsigned int mu_idx = 0; mu_idx < pico.out_mu_sig().size(); mu_idx++) {
		if (pico.out_mu_sig()[mu_idx]) {
			signal_lepton_eta.push_back(pico.out_mu_eta()[mu_idx]);
			signal_lepton_phi.push_back(pico.out_mu_phi()[mu_idx]);
		}
	}
	float min_dilep_dr = 999;
	for (unsigned int lep_idx_1 = 0; lep_idx_1 < signal_lepton_eta.size(); lep_idx_1++) {
		for (unsigned int lep_idx_2 = 0; lep_idx_2 < lep_idx_1; lep_idx_2++) {
			float ll_dr = get_deltar(signal_lepton_eta[lep_idx_1],signal_lepton_phi[lep_idx_1],signal_lepton_eta[lep_idx_2],signal_lepton_phi[lep_idx_2]);
			min_dilep_dr = min_dilep_dr < ll_dr ? min_dilep_dr : ll_dr;
		}
	}
	pico.out_min_dilep_dr() = min_dilep_dr;

	//calculate z_m
	float z_m = -999;
	int zcandidate_idx = -1;
	for (unsigned int dilep_idx = 0; dilep_idx < pico.out_ll_m().size(); dilep_idx++) {
		if (pico.out_ll_charge()[dilep_idx]==0) { //OSSF
			if (TMath::Abs(pico.out_ll_m()[dilep_idx]-91) < TMath::Abs(z_m-91)) {
				z_m = pico.out_ll_m()[dilep_idx];
				zcandidate_idx = dilep_idx;	
			}
		}
	}
	pico.out_z_idx() = zcandidate_idx;

	//calculate lll_m
	if (pico.out_nlep() == 3) {
		TLorentzVector temp_p4, sum_p4;
		sum_p4.SetPtEtaPhiM(0,0,0,0);
		for (unsigned int el_idx = 0; el_idx < pico.out_el_sig().size(); el_idx++) {
			if (pico.out_el_sig()[el_idx]) {
				temp_p4.SetPtEtaPhiM(pico.out_el_pt()[el_idx],pico.out_el_eta()[el_idx],pico.out_el_phi()[el_idx],0.000511);
				sum_p4 = sum_p4 + temp_p4;
			}
		}
		for (unsigned int mu_idx = 0; mu_idx < pico.out_mu_sig().size(); mu_idx++) {
			if (pico.out_mu_sig()[mu_idx]) {
				temp_p4.SetPtEtaPhiM(pico.out_mu_pt()[mu_idx],pico.out_mu_eta()[mu_idx],pico.out_mu_phi()[mu_idx],0.106);
				sum_p4 = sum_p4 + temp_p4;
			}
		}
		pico.out_lll_m() = sum_p4.M();
	}
	else {
		pico.out_lll_m() = -999;
	}

	//calculate l3_mt
	float l3_mt = -999;
	if (!(zcandidate_idx == -1)) {
		//loop over leptons, return M_T for signal lepton that is not z-candidate
		for (unsigned int el_idx = 0; el_idx < pico.out_el_pt().size(); el_idx++) {
			if (pico.out_el_sig()[el_idx]) {
				if (pico.out_ll_lepid()[zcandidate_idx] == 11 && (pico.out_ll_i1()[zcandidate_idx] == static_cast<int>(el_idx) || pico.out_ll_i2()[zcandidate_idx] == static_cast<int>(el_idx))) {
					continue;
				}
				//this is the lepton of interest
				l3_mt = TMath::Sqrt(2.0*pico.out_el_pt()[el_idx]*pico.out_met()*(1.0-TMath::Cos(get_deltaphi(pico.out_el_phi()[el_idx],pico.out_met_phi()))));
			}
		}
		for (unsigned int mu_idx = 0; mu_idx < pico.out_mu_pt().size(); mu_idx++) {
			if (pico.out_mu_sig()[mu_idx]) {
				if (pico.out_ll_lepid()[zcandidate_idx] == 13 && (pico.out_ll_i1()[zcandidate_idx] == static_cast<int>(mu_idx) || pico.out_ll_i2()[zcandidate_idx] == static_cast<int>(mu_idx))) {
					continue;
				}
				//this is the lepton of interest
				l3_mt = TMath::Sqrt(2.0*pico.out_mu_pt()[mu_idx]*pico.out_met()*(1.0-TMath::Cos(get_deltaphi(pico.out_mu_phi()[mu_idx],pico.out_met_phi()))));
			}
		}
	}
	pico.out_l3_mt() = l3_mt;


  return;
}

