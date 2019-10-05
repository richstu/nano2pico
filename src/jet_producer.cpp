#include "jet_producer.hpp"

#include "utilities.hpp"

using namespace std;

JetProducer::JetProducer(int year_){
    year = year_;
}

JetProducer::~JetProducer(){
}

void JetProducer::WriteJets(nano_tree &nano, pico_tree &pico, vector<int> jet_islep_nano_idx){
  pico.out_njet() = 0; 
  for(int ijet(0); ijet<nano.nJet(); ++ijet){
    if (nano.Jet_pt()[ijet] <= JetPtCut) continue;
    if (fabs(nano.Jet_eta()[ijet]) > JetEtaCut) continue;

    if (nano.Jet_jetId()[ijet] < 1) continue; //require just loosest possible ID for now

    pico.out_jet_pt().push_back(nano.Jet_pt()[ijet]);
    pico.out_jet_eta().push_back(nano.Jet_eta()[ijet]);
    pico.out_jet_phi().push_back(nano.Jet_phi()[ijet]);
    pico.out_jet_m().push_back(nano.Jet_mass()[ijet]);
    pico.out_jet_deepcsv().push_back(nano.Jet_btagDeepB()[ijet]);
    pico.out_jet_deepflav().push_back(nano.Jet_btagDeepFlavB()[ijet]);
    pico.out_jet_hflavor().push_back(nano.Jet_hadronFlavour()[ijet]);
    pico.out_jet_pflavor().push_back(nano.Jet_partonFlavour()[ijet]);

    // check overlap with signal leptons
    bool islep = find(jet_islep_nano_idx.begin(), jet_islep_nano_idx.end(), ijet) != jet_islep_nano_idx.end();
    pico.out_jet_islep().push_back(islep);

    //this will be filled in hig_producer or later
    pico.out_jet_h1d().push_back(false);
    pico.out_jet_h2d().push_back(false);

    pico.out_njet()++;
  }
  return;
}
