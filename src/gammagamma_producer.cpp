#include "gammagamma_producer.hpp"

#include "utilities.hpp"
#include "TLorentzVector.h"

using namespace std;

GammaGammaVarProducer::GammaGammaVarProducer(int year_){
    year = year_;
}

GammaGammaVarProducer::~GammaGammaVarProducer(){
}

void GammaGammaVarProducer::WriteGammaGammaVars(pico_tree &pico){
  pico.out_nphotonphoton() = 0;
  if (pico.out_nphoton() <= 1) return;
  TLorentzVector leadPhoton, subPhoton, gammagamma;
  leadPhoton.SetPtEtaPhiM(pico.out_photon_pt()[0], pico.out_photon_eta()[0], pico.out_photon_phi()[0], 0);
  subPhoton.SetPtEtaPhiM(pico.out_photon_pt()[1], pico.out_photon_eta()[1], pico.out_photon_phi()[1], 0);
  gammagamma = leadPhoton + subPhoton;
  pico.out_nphotonphoton()++;
  pico.out_photonphoton_pt().push_back(gammagamma.Pt());
  pico.out_photonphoton_eta().push_back(gammagamma.Eta());
  pico.out_photonphoton_phi().push_back(gammagamma.Phi());
  pico.out_photonphoton_m().push_back(gammagamma.M());
  pico.out_photonphoton_dr().push_back(leadPhoton.DeltaR(subPhoton));
  pico.out_photonphoton_dphi().push_back(leadPhoton.DeltaPhi(subPhoton));
  pico.out_photonphoton_deta().push_back(fabs(leadPhoton.Eta() - subPhoton.Eta()));
  pico.out_photonphoton_ileadph().push_back(0);
  pico.out_photonphoton_isubph().push_back(1);
}
