#include "jet_producer.hpp"

#include "utilities.hpp"
#include <algorithm>

using namespace std;

JetProducer::JetProducer(int year_){
  year = year_;
}

JetProducer::~JetProducer(){
}

void JetProducer::WriteJets(nano_tree &nano, pico_tree &pico, vector<int> jet_islep_nano_idx,
                            const vector<float> &btag_wpts, const vector<float> &btag_df_wpts){
  pico.out_njet() = 0; 
  pico.out_nbl() = 0; pico.out_nbm() = 0; pico.out_nbt() = 0; 
  pico.out_nbdfl() = 0; pico.out_nbdfm() = 0; pico.out_nbdft() = 0; 
  pico.out_pass_jets() = true;
  pico.out_pass_fsjets() = true;
  for(int ijet(0); ijet<nano.nJet(); ++ijet){
    // check overlap with signal leptons
    bool islep = find(jet_islep_nano_idx.begin(), jet_islep_nano_idx.end(), ijet) != jet_islep_nano_idx.end();

    // Fastsim: veto if certain central jets have no matching GenJet as per SUSY recommendation:
    // https://twiki.cern.ch/twiki/bin/view/CMS/SUSRecommendations18#Cleaning_up_of_fastsim_jets_from
    if(nano.Jet_pt()[ijet] > 20 && fabs(nano.Jet_eta()[ijet])>2.5 && 
      nano.Jet_chHEF()[ijet] < 0.1 && nano.Jet_genJetIdx()[ijet]<0) 
      pico.out_pass_fsjets() = false;

    // Fullsim: require just loosest possible ID for now (for all jets, not just central!)
    if (nano.Jet_pt()[ijet] > JetPtCut && nano.Jet_jetId()[ijet] < 1) 
      pico.out_pass_jets() = false;

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
    pico.out_jet_islep().push_back(islep);
    pico.out_jet_id().push_back(nano.Jet_jetId()[ijet]);

    //this will be filled in hig_producer or later
    pico.out_jet_h1d().push_back(false);
    pico.out_jet_h2d().push_back(false);

    if (!islep) {
      pico.out_njet()++;
      if (nano.Jet_btagDeepB()[ijet] > btag_wpts[0]) pico.out_nbl()++; 
      if (nano.Jet_btagDeepB()[ijet] > btag_wpts[1]) pico.out_nbm()++; 
      if (nano.Jet_btagDeepB()[ijet] > btag_wpts[2]) pico.out_nbt()++;
      if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[0]) pico.out_nbdfl()++; 
      if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[1]) pico.out_nbdfm()++; 
      if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[2]) pico.out_nbdft()++; 

    }
  }
  return;
}
