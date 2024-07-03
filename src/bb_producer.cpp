#include "bb_producer.hpp"

#include "utilities.hpp"
#include "TLorentzVector.h"
using namespace std;

BBVarProducer::BBVarProducer(int year_){
    year = year_;
}

BBVarProducer::~BBVarProducer(){
}

void BBVarProducer::WriteBBVars(pico_tree &pico, bool doDeepFlav){
  
  /* Approach: Choosing two highest discriminator b jets

  // get jet 4-vectors ordered by decreasing b-tag discriminator value,
  // also saving their original index in the pico.out_jet* vectors
  vector<pair<int, float>>  ordered_by_discr;
  for (unsigned ijet(0); ijet<pico.out_jet_pt().size(); ijet++) {
    if (pico.out_jet_isgood()[ijet]) {
      float discr = doDeepFlav ? pico.out_jet_deepflav()[ijet] : pico.out_jet_deepcsv()[ijet];
      ordered_by_discr.push_back(make_pair(ijet, discr));
    }
  }
  
  // enough jets to make higgs?
  //if (ordered_by_discr.size()>=2) {
  if (pico.out_jet_pt().size()>=2) {
 
    sort(ordered_by_discr.begin(), ordered_by_discr.end(), 
          [](const pair<int, float> &a, const pair<int, float> &b) -> bool {
            return a.second > b.second;
          });

    // Make higgs candidate with best b tags
    vector<TLorentzVector> jets_lv;
    for (unsigned ijet(0); ijet<2; ijet++) {
      TLorentzVector lv;
      lv.SetPtEtaPhiM(pico.out_jet_pt()[ordered_by_discr[ijet].first], 
                      pico.out_jet_eta()[ordered_by_discr[ijet].first], 
                      pico.out_jet_phi()[ordered_by_discr[ijet].first], 
                      pico.out_jet_m()[ordered_by_discr[ijet].first]);
      jets_lv.push_back(lv);
    }


    // Make higgs candidate with highest b pt 
    vector<TLorentzVector> jets_lv;
    for (unsigned ijet(0); ijet<2; ijet++) {
      TLorentzVector lv;
      lv.SetPtEtaPhiM(pico.out_jet_pt()[ijet], 
                      pico.out_jet_eta()[ijet], 
                      pico.out_jet_phi()[ijet], 
                      pico.out_jet_m()[ijet]);
      jets_lv.push_back(lv);
    }

    */
      

  // Approach: Choosing two highest pt b jets

  vector<pair<int, float>>  ordered_by_pt;
  for (unsigned ijet(0); ijet<pico.out_jet_pt().size(); ijet++) {
    if (pico.out_jet_isgood()[ijet]) ordered_by_pt.push_back(make_pair(ijet, pico.out_jet_pt()[ijet]));
  }
  
  // enough jets to make higgs?
  if (ordered_by_pt.size()>=2) {
 
    sort(ordered_by_pt.begin(), ordered_by_pt.end(), 
          [](const pair<int, float> &a, const pair<int, float> &b) -> bool {
            return a.second > b.second;
          });

    vector<TLorentzVector> jets_lv;
    for (unsigned ijet(0); ijet<2; ijet++) {
      TLorentzVector lv;
      lv.SetPtEtaPhiM(pico.out_jet_pt()[ordered_by_pt[ijet].first], 
                      pico.out_jet_eta()[ordered_by_pt[ijet].first], 
                      pico.out_jet_phi()[ordered_by_pt[ijet].first], 
                      pico.out_jet_m()[ordered_by_pt[ijet].first]);
      jets_lv.push_back(lv);
    }
    
    TLorentzVector higgs = jets_lv[0] + jets_lv[1];

    if (doDeepFlav) {
      pico.out_nbb_df() = 0;
      pico.out_nbb_df()++;
      pico.out_bb_df_pt().push_back(higgs.Pt());
      pico.out_bb_df_eta().push_back(higgs.Eta());
      pico.out_bb_df_phi().push_back(higgs.Phi());
      pico.out_bb_df_m().push_back(higgs.M());
      pico.out_bb_df_dr().push_back(jets_lv[0].DeltaR(jets_lv[1]));
      pico.out_bb_df_dphi().push_back(jets_lv[0].DeltaPhi(jets_lv[1]));
      pico.out_bb_df_deta().push_back(fabs(jets_lv[0].Eta()-jets_lv[1].Eta()));
      //pico.out_bb_df_ileadscorejet().push_back(ordered_by_pt[0].first);
      //pico.out_bb_df_isubscorejet().push_back(ordered_by_pt[1].first);
    } else {
      pico.out_nbb() = 0;
      pico.out_nbb()++;
      pico.out_bb_pt().push_back(higgs.Pt());
      pico.out_bb_eta().push_back(higgs.Eta());
      pico.out_bb_phi().push_back(higgs.Phi());
      pico.out_bb_m().push_back(higgs.M());
      pico.out_bb_dr().push_back(jets_lv[0].DeltaR(jets_lv[1]));
      pico.out_bb_dphi().push_back(jets_lv[0].DeltaPhi(jets_lv[1]));
      pico.out_bb_deta().push_back(fabs(jets_lv[0].Eta()-jets_lv[1].Eta()));
      pico.out_bb_ileadscorejet().push_back(ordered_by_pt[0].first); //here they are lead and sublead pt, rather than lead score jet.
      pico.out_bb_isubscorejet().push_back(ordered_by_pt[1].first); //here they are lead and sublead pt, rather than lead score jet.
    }
  }
}