#include "dilep_producer.hpp"

#include "utilities.hpp"

#include "TLorentzVector.h"

using namespace std;

DileptonProducer::DileptonProducer(int year_){
    year = year_;
}

DileptonProducer::~DileptonProducer(){
}

void DileptonProducer::WriteDileptons(nano_tree &nano, pico_tree &pico, 
                                      std::vector<int> sig_el_nano_idx, std::vector<int> sig_mu_nano_idx) {
  if (pico.out_nmu()<2 && pico.out_nel()<2) return;
  if (pico.out_nmu()>=2)
    for(size_t i(0); i < sig_mu_nano_idx.size(); i++) 
      for(size_t j(i+1); j < sig_mu_nano_idx.size(); j++) {
        int imu1 = sig_mu_nano_idx.at(i);
        int imu2 = sig_mu_nano_idx.at(j);
        if(pico.out_mu_charge()[imu1] + pico.out_mu_charge()[imu2] == 0) {
          TLorentzVector mu1, mu2, dimu;
          mu1.SetPtEtaPhiM(nano.Muon_pt()[imu1], nano.Muon_eta()[imu1],
                           nano.Muon_phi()[imu1],nano.Muon_mass()[imu1]);
          mu2.SetPtEtaPhiM(nano.Muon_pt()[imu2], nano.Muon_eta()[imu2],
                           nano.Muon_phi()[imu2],nano.Muon_mass()[imu2]);
          dimu = mu1 + mu2;
          pico.out_ll_pt().push_back(dimu.Pt());
          pico.out_ll_eta().push_back(dimu.Eta());
          pico.out_ll_phi().push_back(dimu.Phi());
          pico.out_ll_m().push_back(dimu.M());
          pico.out_ll_dr().push_back(mu1.DeltaR(mu2));
          pico.out_ll_dphi().push_back(fabs(mu1.DeltaPhi(mu2)));
          pico.out_ll_deta().push_back(fabs(mu1.Eta()-mu2.Eta()));
          pico.out_ll_lepid().push_back(13);
          pico.out_ll_i1().push_back(imu1);
          pico.out_ll_i2().push_back(imu2);
        }
      }
  if (pico.out_nel()>=2)
    for(size_t i(0); i < sig_el_nano_idx.size(); i++) 
      for(size_t j(i+1); j < sig_el_nano_idx.size(); j++) {
        int iel1 = sig_el_nano_idx.at(i);
        int iel2 = sig_el_nano_idx.at(j);
        if(pico.out_el_charge()[iel1] + pico.out_el_charge()[iel2] == 0) {
          TLorentzVector el1, el2, diel;
          el1.SetPtEtaPhiM(nano.Electron_pt()[iel1] ,nano.Electron_eta()[iel1],
                           nano.Electron_phi()[iel1],nano.Electron_mass()[iel1]);
          el2.SetPtEtaPhiM(nano.Electron_pt()[iel2], nano.Electron_eta()[iel2],
                           nano.Electron_phi()[iel2],nano.Electron_mass()[iel2]);
          diel = el1 + el2;
          pico.out_ll_pt().push_back(diel.Pt());
          pico.out_ll_eta().push_back(diel.Eta());
          pico.out_ll_phi().push_back(diel.Phi());
          pico.out_ll_m().push_back(diel.M());
          pico.out_ll_dr().push_back(el1.DeltaR(el2));
          pico.out_ll_dphi().push_back(fabs(el1.DeltaPhi(el2)));
          pico.out_ll_deta().push_back(fabs(el1.Eta()-el2.Eta()));
          pico.out_ll_lepid().push_back(11);
          pico.out_ll_i1().push_back(iel1);
          pico.out_ll_i2().push_back(iel2);
        }
      }
  return;
}
