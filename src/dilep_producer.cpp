#include "dilep_producer.hpp"

#include "utilities.hpp"

#include "TLorentzVector.h"

using namespace std;

DileptonProducer::DileptonProducer(int year_){
    year = year_;
}

DileptonProducer::~DileptonProducer(){
}


void DileptonProducer::WriteDielectrons(nano_tree &nano, pico_tree &pico, std::vector<int> sig_el_nano_idx){

  if (pico.out_nel()<2) return;

  for (unsigned i(0); i<sig_el_nano_idx.size(); i++){
    TLorentzVector el1; 
    el1.SetPtEtaPhiM(nano.Electron_pt()[i], nano.Electron_eta()[i], nano.Electron_phi()[i], nano.Electron_mass()[i]);
    for (unsigned j(i+1); j<sig_el_nano_idx.size(); j++){
      TLorentzVector diel; 
      diel.SetPtEtaPhiM(nano.Electron_pt()[j], nano.Electron_eta()[j], nano.Electron_phi()[j], nano.Electron_mass()[j]); 
      diel += el1;
      pico.out_elel_pt().push_back(diel.Pt());
      pico.out_elel_eta().push_back(diel.Eta());
      pico.out_elel_phi().push_back(diel.Phi());
      pico.out_elel_m().push_back(diel.M());
    }
  }
  
  return;
}


void DileptonProducer::WriteDimuons(nano_tree &nano, pico_tree &pico, std::vector<int> sig_mu_nano_idx){

  if (pico.out_nmu()<2) return;

  for (unsigned i(0); i<sig_mu_nano_idx.size(); i++){
    TLorentzVector mu1; 
    mu1.SetPtEtaPhiM(nano.Muon_pt()[i], nano.Muon_eta()[i], nano.Muon_phi()[i], nano.Muon_mass()[i]);
    for (unsigned j(i+1); j<sig_mu_nano_idx.size(); j++){
      TLorentzVector dimu; 
      dimu.SetPtEtaPhiM(nano.Muon_pt()[j], nano.Muon_eta()[j], nano.Muon_phi()[j], nano.Muon_mass()[j]); 
      dimu += mu1;
      pico.out_mumu_pt().push_back(dimu.Pt());
      pico.out_mumu_eta().push_back(dimu.Eta());
      pico.out_mumu_phi().push_back(dimu.Phi());
      pico.out_mumu_m().push_back(dimu.M());
    }
  }

  return;
}