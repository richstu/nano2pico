#include "jet_producer.hpp"

#include <algorithm>

#include "TLorentzVector.h"

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

#include "utilities.hpp"

using namespace std;

JetProducer::JetProducer(int year_){
  year = year_;
}

JetProducer::~JetProducer(){
}

vector<int> JetProducer::WriteJets(nano_tree &nano, pico_tree &pico,
                                   vector<int> jet_islep_nano_idx, vector<int> jet_isphoton_nano_idx,
                                   const vector<float> &btag_wpts, const vector<float> &btag_df_wpts, 
                                   bool isZgamma){
  vector<int> sig_jet_nano_idx;
  pico.out_njet() = 0; pico.out_ht() = 0; 
  pico.out_nbl() = 0; pico.out_nbm() = 0; pico.out_nbt() = 0; 
  pico.out_nbdfl() = 0; pico.out_nbdfm() = 0; pico.out_nbdft() = 0; 
  pico.out_pass_jets() = true;
  pico.out_pass_fsjets() = true;
  for(int ijet(0); ijet<nano.nJet(); ++ijet){
    // check overlap with signal leptons (or photons)
    bool islep = find(jet_islep_nano_idx.begin(), jet_islep_nano_idx.end(), ijet) != jet_islep_nano_idx.end();
    bool isphoton = find(jet_isphoton_nano_idx.begin(), jet_isphoton_nano_idx.end(), ijet) != jet_isphoton_nano_idx.end();
    // Fastsim: veto if certain central jets have no matching GenJet as per SUSY recommendation:
    // https://twiki.cern.ch/twiki/bin/view/CMS/SUSRecommendations18#Cleaning_up_of_fastsim_jets_from
    if(nano.Jet_pt()[ijet] > 20 && fabs(nano.Jet_eta()[ijet])>2.5 && 
      nano.Jet_chHEF()[ijet] < 0.1 && nano.Jet_genJetIdx()[ijet]<0) 
      pico.out_pass_fsjets() = false;
    // Fullsim: require just loosest possible ID for now (for all jets, not just central!)
    if (nano.Jet_pt()[ijet] > JetPtCut && nano.Jet_jetId()[ijet] < 1) 
      pico.out_pass_jets() = false;

    if (nano.Jet_pt()[ijet] <= JetPtCut) continue;
    // Use a different eta cut for Zgamma production
    if (fabs(nano.Jet_eta()[ijet])  > (isZgamma ? ZgJetEtaCut:JetEtaCut)) continue;
    pico.out_jet_pt().push_back(nano.Jet_pt()[ijet]);
    pico.out_jet_eta().push_back(nano.Jet_eta()[ijet]);
    pico.out_jet_phi().push_back(nano.Jet_phi()[ijet]);
    pico.out_jet_m().push_back(nano.Jet_mass()[ijet]);
    pico.out_jet_deepcsv().push_back(nano.Jet_btagDeepB()[ijet]);
    pico.out_jet_deepflav().push_back(nano.Jet_btagDeepFlavB()[ijet]);
    pico.out_jet_qgl().push_back(nano.Jet_qgl()[ijet]);
    pico.out_jet_hflavor().push_back(nano.Jet_hadronFlavour()[ijet]);
    pico.out_jet_pflavor().push_back(nano.Jet_partonFlavour()[ijet]);
    pico.out_jet_islep().push_back(islep);
    pico.out_jet_isphoton().push_back(isphoton);
    pico.out_jet_id().push_back(nano.Jet_jetId()[ijet]);
    pico.out_jet_met_dphi().push_back(DeltaPhi(nano.Jet_phi()[ijet], nano.MET_phi()));
    
    // will be overwritten with the overlapping fat jet index, if such exists, in WriteFatJets
    pico.out_jet_fjet_idx().push_back(-999);

    //the jets for the higgs pair with smallest dm will be set to true in hig_producer
    pico.out_jet_h1d().push_back(false);
    pico.out_jet_h2d().push_back(false);

    if (!islep && !isphoton) {
      sig_jet_nano_idx.push_back(ijet);
      pico.out_njet()++;
      pico.out_ht() += nano.Jet_pt()[ijet];
      if (nano.Jet_btagDeepB()[ijet] > btag_wpts[0]) pico.out_nbl()++; 
      if (nano.Jet_btagDeepB()[ijet] > btag_wpts[1]) pico.out_nbm()++; 
      if (nano.Jet_btagDeepB()[ijet] > btag_wpts[2]) pico.out_nbt()++;
      if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[0]) pico.out_nbdfl()++; 
      if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[1]) pico.out_nbdfm()++; 
      if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[2]) pico.out_nbdft()++; 
    }
  } // end jet loop

  return sig_jet_nano_idx;
}

void JetProducer::WriteFatJets(nano_tree &nano, pico_tree &pico){

  pico.out_nfjet() = 0; 
  for(int ifjet(0); ifjet<nano.nFatJet(); ++ifjet){
    // if (nano.FatJet_pt()[ifjet] <= FatJetPtCut) continue;
    if (fabs(nano.FatJet_eta()[ifjet]) > JetEtaCut) continue;

    pico.out_fjet_pt().push_back(nano.FatJet_pt()[ifjet]);
    pico.out_fjet_eta().push_back(nano.FatJet_eta()[ifjet]);
    pico.out_fjet_phi().push_back(nano.FatJet_phi()[ifjet]);
    pico.out_fjet_m().push_back(nano.FatJet_mass()[ifjet]);
    pico.out_fjet_msoftdrop().push_back(nano.FatJet_msoftdrop()[ifjet]);
    // Mass-decorrelated Deep Double B, H->bb vs QCD discriminator, endorsed by BTV
    pico.out_fjet_deep_md_hbb_btv().push_back(nano.FatJet_btagDDBvL()[ifjet]);
    pico.out_fjet_mva_hbb_btv().push_back(nano.FatJet_btagHbb()[ifjet]);
    // Mass-decorrelated DeepAk8, H->bb vs QCD discriminator, endorsed by JME
    pico.out_fjet_deep_md_hbb_jme().push_back(nano.FatJet_deepTagMD_HbbvsQCD()[ifjet]);

    for(unsigned ijet(0); ijet<pico.out_jet_pt().size(); ++ijet){
      if (dR(pico.out_jet_eta()[ijet], nano.FatJet_eta()[ifjet], pico.out_jet_phi()[ijet], nano.FatJet_phi()[ifjet])<0.8)
        pico.out_jet_fjet_idx()[ijet] = ifjet;
    }

    pico.out_nfjet()++;
  }

  //calculate sum of R 1.4 jet masses
  vector<fastjet::PseudoJet> skinny_jets(0);
  //loop over Ak4 jets to add to skinny_jets vector
  for (unsigned int jet_idx(0); jet_idx < pico.out_jet_pt().size(); jet_idx++) {
	//currently includes leptons and photons in clustering
	if (pico.out_jet_pt()[jet_idx]> JetPtCut || pico.out_jet_islep()[jet_idx] || pico.out_jet_isphoton()[jet_idx]) {
		TLorentzVector jet_4p;
		jet_4p.SetPtEtaPhiM(pico.out_jet_pt()[jet_idx],pico.out_jet_eta()[jet_idx],pico.out_jet_phi()[jet_idx],pico.out_jet_m()[jet_idx]);
		const fastjet::PseudoJet this_pseudo_jet(jet_4p.Px(), jet_4p.Py(), jet_4p.Pz(), jet_4p.E());
		skinny_jets.push_back(this_pseudo_jet);
	}
  } //loop over Ak4 Jets
  //cluster into R 1.4 jets with fastjet
  fastjet::JetDefinition jet_definition(fastjet::antikt_algorithm, 1.4);
  fastjet::ClusterSequence cluster_sequence(skinny_jets, jet_definition);
  vector<fastjet::PseudoJet> fat_jets = cluster_sequence.inclusive_jets();
  pico.out_mj14() = 0;
  for (unsigned int fat_jet_idx = 0; fat_jet_idx < fat_jets.size(); fat_jet_idx++) {
	  pico.out_mj14() += fat_jets[fat_jet_idx].m();
  } //loop over Ak14 jets

  return;
}

void JetProducer::WriteJetSystemPt(nano_tree &nano, pico_tree &pico, 
                                   vector<int> &sig_jet_nano_idx, const float &btag_wpt) {

  TLorentzVector jetsys_v4, jetsys_nob_v4;
  int njet_nob(0);
  for (auto &idx: sig_jet_nano_idx) {
    TLorentzVector ijet_v4;
    ijet_v4.SetPtEtaPhiM(nano.Jet_pt()[idx], nano.Jet_eta()[idx], nano.Jet_phi()[idx], nano.Jet_mass()[idx]);
    jetsys_v4 += ijet_v4;

    if (nano.Jet_btagDeepB()[idx] <= btag_wpt){
      njet_nob++;
      jetsys_nob_v4 += ijet_v4;
    }
  }

  if (sig_jet_nano_idx.size()>0) {
    pico.out_jetsys_pt() = jetsys_v4.Pt();
    pico.out_jetsys_eta() = jetsys_v4.Eta();
    pico.out_jetsys_phi() = jetsys_v4.Phi();
    pico.out_jetsys_m() = jetsys_v4.M();
    if (njet_nob>0) {
      pico.out_jetsys_nob_pt() = jetsys_nob_v4.Pt();
      pico.out_jetsys_nob_eta() = jetsys_nob_v4.Eta();
      pico.out_jetsys_nob_phi() = jetsys_nob_v4.Phi();
      pico.out_jetsys_nob_m() = jetsys_nob_v4.M();
    }
  }
  return;
}
