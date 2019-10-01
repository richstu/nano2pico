#include "jet_producer.hpp"

#include "utilities.hpp"

using namespace std;

JetProducer::JetProducer(int year_){
    year = year_;
}

JetProducer::~JetProducer(){
}

void JetProducer::WriteJets(nano_tree &nano, pico_tree &pico, 
                            vector<int> sig_el_nano_idx, 
                            vector<int> sig_mu_nano_idx){
  pico.out_njet() = 0; 
  for(int ijet(0); ijet<nano.nJet(); ++ijet){
    if (nano.Jet_pt()[ijet] <= JetPtCut) continue;
    if (fabs(nano.Jet_eta()[ijet]) > JetEtaCut) continue;

    pico.out_jet_pt().push_back(nano.Jet_pt()[ijet]);
    pico.out_jet_eta().push_back(nano.Jet_eta()[ijet]);
    pico.out_jet_phi().push_back(nano.Jet_phi()[ijet]);
    pico.out_jet_m().push_back(nano.Jet_mass()[ijet]);
    pico.out_jet_deepcsv().push_back(nano.Jet_btagDeepB()[ijet]);
    pico.out_jet_deepflav().push_back(nano.Jet_btagDeepFlavB()[ijet]);
    pico.out_jet_hflavor().push_back(nano.Jet_hadronFlavour()[ijet]);
    pico.out_jet_pflavor().push_back(nano.Jet_partonFlavour()[ijet]);

    // check overlap with signal leptons
    pico.out_jet_islep().push_back(false);
    for(auto &iel: sig_el_nano_idx){
      if (nano.Electron_isPFcand()[iel]) {
        if (iel==nano.Jet_electronIdx1()[ijet] || iel==nano.Jet_electronIdx2()[ijet]) 
          pico.out_jet_islep().back() = true;
      } else {
        if (dR(nano.Electron_eta()[iel], nano.Jet_eta()[ijet], nano.Electron_phi()[iel], nano.Jet_phi()[ijet]))
          pico.out_jet_islep().back() = true;
      }
    }
    for(auto &imu: sig_mu_nano_idx){
      if (nano.Muon_isPFcand()[imu]) {
        if (imu==nano.Jet_muonIdx1()[ijet] || imu==nano.Jet_muonIdx2()[ijet]) 
          pico.out_jet_islep().back() = true;
      } else {
        if (dR(nano.Muon_eta()[imu], nano.Jet_eta()[ijet], nano.Muon_phi()[imu], nano.Jet_phi()[ijet]))
          pico.out_jet_islep().back() = true;
      }
    }
    
    //this will be filled in hig_producer or later
    pico.out_jet_h1d().push_back(false);
    pico.out_jet_h2d().push_back(false);

    pico.out_njet()++;
  }
  return;
}
