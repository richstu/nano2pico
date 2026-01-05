#include "dilep_producer.hpp"
#include "utilities.hpp"

#include "TLorentzVector.h"

using namespace std;

DileptonProducer::DileptonProducer(int year_){
    year = year_;
}

DileptonProducer::~DileptonProducer(){

}

void DileptonProducer::WriteDileptons(pico_tree &pico, bool is_signal) {
  //note el_producer and mu_producer must be run first
  pico.out_nll() = 0;
  if (pico.out_nmu() < 2 && pico.out_nel() < 2 && !is_signal) return;
  double mindm(99999), zmass(91.1876);
  TLorentzVector l1err, l2err;
  double ptl1err, ptl2err;
  double dml1, dml2;
  int dilep_charge;

  // loop over systematic variations
  // order: nominal, elscaleup, elscaledn, elresup, elresdn, muscaleup, 
  // muscaledn, muresup, muresdn
  unsigned n_variations = 1;
  if (is_signal)
    n_variations = 9;
  pico.out_ll_pt().resize(n_variations, 0.0);
  pico.out_ll_eta().resize(n_variations, 0.0);
  pico.out_ll_phi().resize(n_variations, 0.0);
  pico.out_ll_m().resize(n_variations, 0.0);
  pico.out_ll_l1_masserr().resize(n_variations, 0.0);
  pico.out_ll_l2_masserr().resize(n_variations, 0.0);
  pico.out_ll_dr().resize(n_variations, 0.0);
  pico.out_ll_dphi().resize(n_variations, 0.0);
  pico.out_ll_deta().resize(n_variations, 0.0);
  pico.out_ll_lepid().resize(n_variations, -1);
  pico.out_ll_i1().resize(n_variations, -1);
  pico.out_ll_i2().resize(n_variations, -1);
  pico.out_ll_l1_pt().resize(n_variations, 0.0);
  pico.out_ll_l2_pt().resize(n_variations, 0.0);
  for (unsigned ivar = 0; ivar < n_variations; ivar++) {
    mindm = 99999.0;
    for (unsigned imu1 = 0; imu1 < pico.out_mu_pt().size(); imu1++) {
      for (unsigned imu2 = imu1+1; imu2 < pico.out_mu_pt().size(); imu2++) {
        bool mu1_sig = pico.out_mu_sig()[imu1];
        bool mu2_sig = pico.out_mu_sig()[imu2];
        float mu1_pt = pico.out_mu_pt()[imu1];
        float mu2_pt = pico.out_mu_pt()[imu2];
        if (ivar == 5) {
          mu1_sig = pico.out_sys_mu_sig_scaleup()[imu1];
          mu2_sig = pico.out_sys_mu_sig_scaleup()[imu2];
          mu1_pt = pico.out_sys_mu_pt_scaleup()[imu1];
          mu2_pt = pico.out_sys_mu_pt_scaleup()[imu2];
        }
        else if (ivar == 6) {
          mu1_sig = pico.out_sys_mu_sig_scaledn()[imu1];
          mu2_sig = pico.out_sys_mu_sig_scaledn()[imu2];
          mu1_pt = pico.out_sys_mu_pt_scaledn()[imu1];
          mu2_pt = pico.out_sys_mu_pt_scaledn()[imu2];
        }
        else if (ivar == 7) {
          mu1_sig = pico.out_sys_mu_sig_resup()[imu1];
          mu2_sig = pico.out_sys_mu_sig_resup()[imu2];
          mu1_pt = pico.out_sys_mu_pt_resup()[imu1];
          mu2_pt = pico.out_sys_mu_pt_resup()[imu2];
        }
        else if (ivar == 8) {
          mu1_sig = pico.out_sys_mu_sig_resdn()[imu1];
          mu2_sig = pico.out_sys_mu_sig_resdn()[imu2];
          mu1_pt = pico.out_sys_mu_pt_resdn()[imu1];
          mu2_pt = pico.out_sys_mu_pt_resdn()[imu2];
        }
        if (!(mu1_sig && mu2_sig)) continue;

        dilep_charge = pico.out_mu_charge()[imu1]+pico.out_mu_charge()[imu2];
        if (dilep_charge != 0) {
          continue;
        }
        if (ivar == 0) {
          pico.out_nll() = 1;
        }


        TLorentzVector mu1, mu2, dimu;
        mu1.SetPtEtaPhiM(mu1_pt, pico.out_mu_eta()[imu1],
                         pico.out_mu_phi()[imu1], 0.10566);
        mu2.SetPtEtaPhiM(mu2_pt, pico.out_mu_eta()[imu2],
                         pico.out_mu_phi()[imu2], 0.10566);
        dimu = mu1 + mu2;
        if (abs(dimu.M() - zmass) < mindm) { 
          mindm = abs(dimu.M() - zmass);

          ptl1err = pico.out_mu_ptErr()[imu1];
          ptl2err = pico.out_mu_ptErr()[imu2];
          l1err.SetPtEtaPhiM(mu1_pt + ptl1err, pico.out_mu_eta()[imu1], 
                             pico.out_mu_phi()[imu1], 0.10566);
          l2err.SetPtEtaPhiM(mu2_pt + ptl2err, pico.out_mu_eta()[imu2], 
                             pico.out_mu_phi()[imu2], 0.10566);
          dml1 = (l1err + mu2).M() - dimu.M();
          dml2 = (mu1 + l2err).M() - dimu.M();

          pico.out_ll_pt()[ivar] = dimu.Pt();
          pico.out_ll_eta()[ivar] = dimu.Eta();
          pico.out_ll_phi()[ivar] = dimu.Phi();
          pico.out_ll_m()[ivar] = dimu.M();
          pico.out_ll_dr()[ivar] = mu1.DeltaR(mu2);
          pico.out_ll_dphi()[ivar] = fabs(mu1.DeltaPhi(mu2));
          pico.out_ll_deta()[ivar] = fabs(mu1.Eta()-mu2.Eta());
          pico.out_ll_lepid()[ivar] = 13;
          pico.out_ll_i1()[ivar] = imu1;
          pico.out_ll_i2()[ivar] = imu2;
          pico.out_ll_l1_masserr()[ivar] = dml1;
          pico.out_ll_l2_masserr()[ivar] = dml2;
          pico.out_ll_l1_pt()[ivar] = mu1_pt;
          pico.out_ll_l2_pt()[ivar] = mu2_pt;

        } // best Z candidate
      } // loop over muon 2
    } // loop over muon 1
    for (unsigned iel1 = 0; iel1 < pico.out_el_pt().size(); iel1++) {
      for (unsigned iel2 = iel1+1; iel2 < pico.out_el_pt().size(); iel2++) {
        bool el1_sig = pico.out_el_sig()[iel1];
        bool el2_sig = pico.out_el_sig()[iel2];
        float el1_pt = pico.out_el_pt()[iel1];
        float el2_pt = pico.out_el_pt()[iel2];
        if (ivar == 1) {
          el1_sig = pico.out_sys_el_sig_scaleup()[iel1];
          el2_sig = pico.out_sys_el_sig_scaleup()[iel2];
          el1_pt = pico.out_sys_el_pt_scaleup()[iel1];
          el2_pt = pico.out_sys_el_pt_scaleup()[iel2];
        }
        else if (ivar == 2) {
          el1_sig = pico.out_sys_el_sig_scaledn()[iel1];
          el2_sig = pico.out_sys_el_sig_scaledn()[iel2];
          el1_pt = pico.out_sys_el_pt_scaledn()[iel1];
          el2_pt = pico.out_sys_el_pt_scaledn()[iel2];
        }
        else if (ivar == 3) {
          el1_sig = pico.out_sys_el_sig_resup()[iel1];
          el2_sig = pico.out_sys_el_sig_resup()[iel2];
          el1_pt = pico.out_sys_el_pt_resup()[iel1];
          el2_pt = pico.out_sys_el_pt_resup()[iel2];
        }
        else if (ivar == 4) {
          el1_sig = pico.out_sys_el_sig_resdn()[iel1];
          el2_sig = pico.out_sys_el_sig_resdn()[iel2];
          el1_pt = pico.out_sys_el_pt_resdn()[iel1];
          el2_pt = pico.out_sys_el_pt_resdn()[iel2];
        }
        if (!(el1_sig && el2_sig)) continue;

        dilep_charge = pico.out_el_charge()[iel1]+pico.out_el_charge()[iel2];
        if (dilep_charge != 0) {
          continue;
        }
        if (ivar == 0) {
          pico.out_nll() = 1;
        }

        TLorentzVector el1, el2, diel;
        el1.SetPtEtaPhiM(el1_pt, pico.out_el_eta()[iel1],
                         pico.out_el_phi()[iel1], 0.000511);
        el2.SetPtEtaPhiM(el2_pt, pico.out_el_eta()[iel2],
                         pico.out_el_phi()[iel2], 0.000511);
        diel = el1 + el2;
        if (abs(diel.M() - zmass) < mindm) { 
          mindm = abs(diel.M() - zmass);

          ptl1err = pico.out_el_energyErr()[iel1] * el1.Pt() / el1.P();
          ptl2err = pico.out_el_energyErr()[iel2] * el2.Pt() / el2.P();
          l1err.SetPtEtaPhiM(el1_pt + ptl1err, pico.out_el_eta()[iel1], 
                             pico.out_el_phi()[iel1], 0.000511);
          l2err.SetPtEtaPhiM(el2_pt + ptl2err, pico.out_el_eta()[iel2], 
                             pico.out_el_phi()[iel2], 0.000511);
          dml1 = (l1err + el2).M() - diel.M();
          dml2 = (el1 + l2err).M() - diel.M();

          pico.out_ll_pt()[ivar] = diel.Pt();
          pico.out_ll_eta()[ivar] = diel.Eta();
          pico.out_ll_phi()[ivar] = diel.Phi();
          pico.out_ll_m()[ivar] = diel.M();
          pico.out_ll_dr()[ivar] = el1.DeltaR(el2);
          pico.out_ll_dphi()[ivar] = fabs(el1.DeltaPhi(el2));
          pico.out_ll_deta()[ivar] = fabs(el1.Eta()-el2.Eta());
          pico.out_ll_lepid()[ivar] = 11;
          pico.out_ll_i1()[ivar] = iel1;
          pico.out_ll_i2()[ivar] = iel2;
          pico.out_ll_l1_masserr()[ivar] = dml1;
          pico.out_ll_l2_masserr()[ivar] = dml2;
          pico.out_ll_l1_pt()[ivar] = el1_pt;
          pico.out_ll_l2_pt()[ivar] = el2_pt;

        } // best Z candidate
      } // loop over electron 2
    } //loop over electron 1
  } // loop over variations

}
