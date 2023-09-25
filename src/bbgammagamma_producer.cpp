#include "bbgammagamma_producer.hpp"

#include "utilities.hpp"
#include "TLorentzVector.h"

using namespace std;

BBGammaGammaVarProducer::BBGammaGammaVarProducer(int year_){
    year = year_;
}

BBGammaGammaVarProducer::~BBGammaGammaVarProducer(){
}

void BBGammaGammaVarProducer::WriteBBGammaGammaVars(pico_tree &pico){

  if (pico.out_nphotonphoton()==0 || pico.out_nbb()==0) return;

  pico.out_nbbphotonphoton() = 0;
  if (pico.out_nbb()>0) {
    TLorentzVector bb_p4, gammagamma_p4;
    bb_p4.SetPtEtaPhiM(pico.out_bb_pt()[0], pico.out_bb_eta()[0], pico.out_bb_phi()[0], pico.out_bb_m()[0]);
    gammagamma_p4.SetPtEtaPhiM(pico.out_photonphoton_pt()[0], pico.out_photonphoton_eta()[0], pico.out_photonphoton_phi()[0], pico.out_photonphoton_m()[0]);
    pico.out_nbbphotonphoton()++;
    pico.out_bbphotonphoton_dm().push_back(fabs(bb_p4.M()-gammagamma_p4.M()));
    pico.out_bbphotonphoton_am().push_back((bb_p4.M()+gammagamma_p4.M())/2);
    pico.out_bbphotonphoton_drmax().push_back(max(pico.out_bb_dr()[0], pico.out_photonphoton_dr()[0]));
  }
  pico.out_nbbphotonphoton_df() = 0;
  if (pico.out_nbb_df()>0) {
    TLorentzVector bb_df_p4, gammagamma_p4;
    bb_df_p4.SetPtEtaPhiM(pico.out_bb_df_pt()[0], pico.out_bb_df_eta()[0], pico.out_bb_df_phi()[0], pico.out_bb_df_m()[0]);
    gammagamma_p4.SetPtEtaPhiM(pico.out_photonphoton_pt()[0], pico.out_photonphoton_eta()[0], pico.out_photonphoton_phi()[0], pico.out_photonphoton_m()[0]);
    pico.out_nbbphotonphoton_df()++;
    pico.out_bbphotonphoton_df_dm().push_back(fabs(bb_df_p4.M()-gammagamma_p4.M()));
    pico.out_bbphotonphoton_df_am().push_back((bb_df_p4.M()+gammagamma_p4.M())/2);
    pico.out_bbphotonphoton_df_drmax().push_back(max(pico.out_bb_df_dr()[0], pico.out_photonphoton_dr()[0]));
  }

  

}
