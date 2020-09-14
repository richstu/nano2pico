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

void MetProducer::WriteMet(nano_tree &nano, pico_tree &pico){
  // Copy MET and ME ISR directly from NanoAOD
  if (year == 2017) {
    pico.out_met()         = nano.METFixEE2017_pt();
    pico.out_met_phi()     = nano.METFixEE2017_phi();
  } else {
    pico.out_met()         = nano.MET_pt();
    pico.out_met_phi()     = nano.MET_phi();
  }
  pico.out_met_calo()    = nano.CaloMET_pt();
  pico.out_met_tru()     = nano.GenMET_pt();
  pico.out_met_tru_phi() = nano.GenMET_phi();
  pico.out_ht_isr_me()   = nano.LHE_HTIncoming();
}
