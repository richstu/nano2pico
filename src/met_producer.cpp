#include "met_producer.hpp"

#include <algorithm>
#include <iomanip> 

#include "TLorentzVector.h"

#include "utilities.hpp"

using namespace std;

MetProducer::MetProducer(int year_, bool isData_, bool verbose_){
  year = year_;
  isData = isData_;
  verbose = verbose_;
}

MetProducer::~MetProducer(){
}

void MetProducer::WriteMet(nano_tree &nano, pico_tree &pico, bool isFastsim, bool isSignal, bool isUL){
  float MET_pt, MET_phi;
  getMETWithJEC(nano, year, isFastsim, MET_pt, MET_phi, isUL);
  // Copy MET and ME ISR directly from NanoAOD
  pico.out_met()         = MET_pt;
  pico.out_met_phi()     = MET_phi;
  pico.out_met_calo()    = nano.CaloMET_pt();
  pico.out_met_tru()     = nano.GenMET_pt();
  pico.out_met_tru_phi() = nano.GenMET_phi();
  pico.out_ht_isr_me()   = nano.LHE_HTIncoming();

  //if (!isData_) {
  //currently apply only to fastsim
  if (isSignal) {
    pico.out_sys_met().resize(4,-999.0);
    pico.out_sys_met_phi().resize(4,-999.0);
    if (year==2017 && isFastsim && !isUL) {
      pico.out_sys_met()[0] = nano.METFixEE2017_T1Smear_pt_jerUp();
      pico.out_sys_met()[1] = nano.METFixEE2017_T1Smear_pt_jerDown();
      pico.out_sys_met()[2] = nano.METFixEE2017_T1_pt_jesTotalUp();
      pico.out_sys_met()[3] = nano.METFixEE2017_T1_pt_jesTotalDown();
      pico.out_sys_met_phi()[0] = nano.METFixEE2017_T1Smear_phi_jerUp();
      pico.out_sys_met_phi()[1] = nano.METFixEE2017_T1Smear_phi_jerDown();
      pico.out_sys_met_phi()[2] = nano.METFixEE2017_T1_phi_jesTotalUp();
      pico.out_sys_met_phi()[3] = nano.METFixEE2017_T1_phi_jesTotalDown();
    }
    else {
      pico.out_sys_met()[0] = nano.MET_T1Smear_pt_jerUp();
      pico.out_sys_met()[1] = nano.MET_T1Smear_pt_jerDown();
      pico.out_sys_met()[2] = nano.MET_T1_pt_jesTotalUp();
      pico.out_sys_met()[3] = nano.MET_T1_pt_jesTotalDown();
      pico.out_sys_met_phi()[0] = nano.MET_T1Smear_phi_jerUp();
      pico.out_sys_met_phi()[1] = nano.MET_T1Smear_phi_jerDown();
      pico.out_sys_met_phi()[2] = nano.MET_T1_phi_jesTotalUp();
      pico.out_sys_met_phi()[3] = nano.MET_T1_phi_jesTotalDown();
    }
  }
}
