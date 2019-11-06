#include "zgamma_producer.hpp"

#include "utilities.hpp"
#include "TLorentzVector.h"

using namespace std;

ZGammaVarProducer::ZGammaVarProducer(int year_){
    year = year_;
}

ZGammaVarProducer::~ZGammaVarProducer(){
}

void ZGammaVarProducer::WriteZGammaVars(pico_tree &pico){
  if(pico.out_njet() > 1) {
    TLorentzVector j1, j2, dijet;
    j1.SetPtEtaPhiM(pico.out_jet_pt()[0], pico.out_jet_eta()[0], 
                    pico.out_jet_phi()[0],pico.out_jet_m()[0]);
    j2.SetPtEtaPhiM(pico.out_jet_pt()[1], pico.out_jet_eta()[1], 
                    pico.out_jet_phi()[1],pico.out_jet_m()[1]);
    dijet = j1 + j2;
    pico.out_dijet_pt()   = dijet.Pt();
    pico.out_dijet_eta()  = dijet.Eta();
    pico.out_dijet_phi()  = dijet.Phi();
    pico.out_dijet_m()    = dijet.M();
    pico.out_dijet_dr()   = j1.DeltaR(j2);
    pico.out_dijet_dphi() = j1.DeltaPhi(j2);
    pico.out_dijet_deta() = fabs(j1.Eta() - j2.Eta());
  }
  if (pico.out_ll_pt().size() == 0 || pico.out_nphoton() == 0) return;
  for(size_t ill(0); ill < pico.out_ll_pt().size(); ill++) 
    for(int igamma(0); igamma < pico.out_nphoton(); igamma++) {
      TLorentzVector dilep, photon, llg;
      dilep.SetPtEtaPhiM(pico.out_ll_pt()[ill], pico.out_ll_eta()[ill],
                         pico.out_ll_phi()[ill],pico.out_ll_m()[ill]);
      photon.SetPtEtaPhiM(pico.out_photon_pt()[igamma], pico.out_photon_eta()[igamma],
                          pico.out_photon_phi()[igamma],0);
      llg = dilep + photon;
      pico.out_llphoton_pt().push_back(llg.Pt());
      pico.out_llphoton_eta().push_back(llg.Eta());
      pico.out_llphoton_phi().push_back(llg.Phi());
      pico.out_llphoton_m().push_back(llg.M());
      pico.out_llphoton_dr().push_back(dilep.DeltaR(photon));
      pico.out_llphoton_dphi().push_back(dilep.DeltaPhi(photon));
      pico.out_llphoton_deta().push_back(fabs(dilep.Eta() - photon.Eta()));
      pico.out_llphoton_iph().push_back(igamma);
      pico.out_llphoton_ill().push_back(ill);
      TLorentzVector lminus, lplus;
      if(pico.out_ll_lepid()[ill] == 11) {
        int iel1 = pico.out_ll_i1()[ill];
        int iel2 = pico.out_ll_i2()[ill];
        if(pico.out_el_charge()[iel1] < 0) {
          lminus.SetPtEtaPhiM(pico.out_el_pt()[iel1], pico.out_el_eta()[iel1],
                              pico.out_el_phi()[iel1],0.0005);
          lplus.SetPtEtaPhiM(pico.out_el_pt()[iel2], pico.out_el_eta()[iel2],
                             pico.out_el_phi()[iel2],0.0005);
        }
        else {
          lplus.SetPtEtaPhiM(pico.out_el_pt()[iel1], pico.out_el_eta()[iel1],
                             pico.out_el_phi()[iel1],0.0005);
          lminus.SetPtEtaPhiM(pico.out_el_pt()[iel2], pico.out_el_eta()[iel2],
                              pico.out_el_phi()[iel2],0.0005);
        }
      }
      else {
        int imu1 = pico.out_ll_i1()[ill];
        int imu2 = pico.out_ll_i2()[ill];
        if(pico.out_mu_charge()[imu1] < 0) {
          lminus.SetPtEtaPhiM(pico.out_mu_pt()[imu1], pico.out_mu_eta()[imu1],
                              pico.out_mu_phi()[imu1],0.105);
          lplus.SetPtEtaPhiM(pico.out_mu_pt()[imu2], pico.out_mu_eta()[imu2],
                             pico.out_mu_phi()[imu2],0.105);
        }
        else {
          lplus.SetPtEtaPhiM(pico.out_mu_pt()[imu1], pico.out_mu_eta()[imu1],
                             pico.out_mu_phi()[imu1],0.105);
          lminus.SetPtEtaPhiM(pico.out_mu_pt()[imu2], pico.out_mu_eta()[imu2],
                              pico.out_mu_phi()[imu2],0.105);
        }
      }
      pico.out_llphoton_costhj().push_back(cosThetaJeff(lminus,lplus,photon));
  }
  return;
}

