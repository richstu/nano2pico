#include "hig_producer.hpp"

#include <iomanip> 
#include <algorithm>
#include <functional>

#include "TLorentzVector.h"

#include "utilities.hpp"

using namespace std;

HigVarProducer::HigVarProducer(int year_){
  year = year_;
}

HigVarProducer::~HigVarProducer(){
}

void HigVarProducer::WriteHigVars(pico_tree &pico, bool doDeepFlav, bool isSignal,
                                  vector<HiggsConstructionVariables> sys_higvars, float nanoaod_version){

  // get jet 4-vectors ordered by decreasing b-tag discriminator value,
  // also saving their original index in the pico.out_jet* vectors
  vector<pair<int, float>>  ordered_by_discr;
  vector<pair<int, float>>  ordered_by_pt; //* new addition
  
  for (unsigned ijet(0); ijet<pico.out_jet_pt().size(); ijet++) {
    if (pico.out_jet_isgood()[ijet]) {
      float discr = -999;
      if (nanoaod_version+0.01 < 11.9) {
        discr = doDeepFlav ? pico.out_jet_deepflav()[ijet] : pico.out_jet_deepcsv()[ijet];
      } else {
        discr = doDeepFlav ? pico.out_jet_deepflav()[ijet] : pico.out_jet_btagpnetb()[ijet];
      }
      ordered_by_discr.push_back(make_pair(ijet, discr));
      ordered_by_pt.push_back(make_pair(ijet, pico.out_jet_pt()[ijet])); //* new addition
    }
  }

  if (ordered_by_pt.size()>=1) {
    // ideally this should be in jetmet_producer.cpp but it will be moved once I understand what is going on in there
    sort(ordered_by_pt.begin(), ordered_by_pt.end(), //* new addition
          [](const pair<int, float> &a, const pair<int, float> &b) -> bool {
            return a.second > b.second;
          }); 
    for (unsigned ijet(0); ijet<pico.out_jet_pt().size(); ijet++){
      pico.out_jet_ordered_pt_indices().push_back(ordered_by_pt[ijet].first);
    }    
  }

  // enough jets to make two higgses?
  if (ordered_by_discr.size()>=4) {

    sort(ordered_by_discr.begin(), ordered_by_discr.end(), 
          [](const pair<int, float> &a, const pair<int, float> &b) -> bool {
            return a.second > b.second;
          });
   
    vector<TLorentzVector> jets_lv;
    for (unsigned ijet(0); ijet<4; ijet++) {
      TLorentzVector lv;
      lv.SetPtEtaPhiM(pico.out_jet_pt()[ordered_by_discr[ijet].first], 
                      pico.out_jet_eta()[ordered_by_discr[ijet].first], 
                      pico.out_jet_phi()[ordered_by_discr[ijet].first], 
                      pico.out_jet_m()[ordered_by_discr[ijet].first]);
      jets_lv.push_back(lv);
    }

    //make possible combinations
    vector<vector<unsigned>>  hcombs = {{0,1,2,3},{0,2,1,3},{0,3,1,2}};
    vector<float> dm, am, drmax;
    unsigned icomb_min_dm(-999);
    for (unsigned ic(0); ic<hcombs.size(); ic++){
      TLorentzVector h1 = jets_lv[hcombs[ic][0]] + jets_lv[hcombs[ic][1]]; 
      TLorentzVector h2 = jets_lv[hcombs[ic][2]] + jets_lv[hcombs[ic][3]]; 

      float idm = fabs(h1.M()-h2.M());
      float iam = (h1.M()+h2.M())/2.;
      float idr1 = dR(jets_lv[hcombs[ic][0]].Eta(), jets_lv[hcombs[ic][1]].Eta(), 
                      jets_lv[hcombs[ic][0]].Phi(), jets_lv[hcombs[ic][1]].Phi());
      float idr2 = dR(jets_lv[hcombs[ic][2]].Eta(), jets_lv[hcombs[ic][3]].Eta(), 
                      jets_lv[hcombs[ic][2]].Phi(), jets_lv[hcombs[ic][3]].Phi());
      float idrmax = idr1 > idr2 ? idr1 : idr2;

      // insert them in order of increasing hig_dm
      unsigned pos = 0;
      if (doDeepFlav) {
        for (unsigned i(0); i< pico.out_hig_df_cand_dm().size(); i++){
          if (idm > pico.out_hig_df_cand_dm()[i]) pos = i+1;
        }
        pico.out_hig_df_cand_dm().insert(pico.out_hig_df_cand_dm().begin()+pos, idm);
        pico.out_hig_df_cand_am().insert(pico.out_hig_df_cand_am().begin()+pos, iam);
        pico.out_hig_df_cand_drmax().insert(pico.out_hig_df_cand_drmax().begin()+pos, idrmax);
        if (pos==0) icomb_min_dm = ic; //* new addition       
      } else {
        for (unsigned i(0); i< pico.out_hig_cand_dm().size(); i++){
          if (idm > pico.out_hig_cand_dm()[i]) pos = i+1;
        }
        pico.out_hig_cand_dm().insert(pico.out_hig_cand_dm().begin()+pos, idm);
        pico.out_hig_cand_am().insert(pico.out_hig_cand_am().begin()+pos, iam);
        pico.out_hig_cand_drmax().insert(pico.out_hig_cand_drmax().begin()+pos, idrmax);

        //save index of jet combination with smallest dm in order to mark the utilized jets (DeepCSV only)
        if (pos==0) icomb_min_dm = ic;
      }
    }

    if (!doDeepFlav){
      // set the jet h1d/h2d variables indicating that the jet was used in the hig pair with smallest dm
      pico.out_jet_h1d()[ordered_by_discr[hcombs[icomb_min_dm][0]].first] = true;
      pico.out_jet_h1d()[ordered_by_discr[hcombs[icomb_min_dm][1]].first] = true;
      pico.out_jet_h2d()[ordered_by_discr[hcombs[icomb_min_dm][2]].first] = true;
      pico.out_jet_h2d()[ordered_by_discr[hcombs[icomb_min_dm][3]].first] = true;
      pico.out_jet_h1_indices().push_back(ordered_by_discr[hcombs[icomb_min_dm][0]].first);
      pico.out_jet_h1_indices().push_back(ordered_by_discr[hcombs[icomb_min_dm][1]].first);
      pico.out_jet_h2_indices().push_back(ordered_by_discr[hcombs[icomb_min_dm][2]].first);
      pico.out_jet_h2_indices().push_back(ordered_by_discr[hcombs[icomb_min_dm][3]].first);
    } else{ //* new addition (save indices for deepflav)
      pico.out_jet_h1_df_indices().push_back(ordered_by_discr[hcombs[icomb_min_dm][0]].first);
      pico.out_jet_h1_df_indices().push_back(ordered_by_discr[hcombs[icomb_min_dm][1]].first);
      pico.out_jet_h2_df_indices().push_back(ordered_by_discr[hcombs[icomb_min_dm][2]].first);
      pico.out_jet_h2_df_indices().push_back(ordered_by_discr[hcombs[icomb_min_dm][3]].first);
    }

  } //if at least 4 good jets

  //make systematic variations. for now leave NanoAODv12 and NanoAODv9 (No JEC corrections for fastsim v9 signal -- Higgsino)
  if (isSignal && !doDeepFlav && (nanoaod_version+0.01)<8) {
    pico.out_sys_hig_cand_dm().resize(4, -999.0);
    pico.out_sys_hig_cand_am().resize(4, -999.0);
    pico.out_sys_hig_cand_drmax().resize(4, -999.0);

    for (unsigned ijec(0); ijec < 4; ijec++) {
      // enough jets to make two higgses?
      if (sys_higvars[ijec].jet_lv.size()>=4) {

        vector<pair<TLorentzVector, float>> sys_ord_jet;
        for (unsigned ijet(0); ijet<sys_higvars[ijec].jet_lv.size(); ijet++) {
          sys_ord_jet.push_back(make_pair(sys_higvars[ijec].jet_lv[ijet], 
          sys_higvars[ijec].jet_deepcsv[ijet]));
        }

        sort(sys_ord_jet.begin(), sys_ord_jet.end(), 
             [](const pair<TLorentzVector, float> &a, const pair<TLorentzVector, float> &b) -> bool {
               return a.second > b.second;
             });

        //make possible combinations
        vector<vector<unsigned>>  hcombs = {{0,1,2,3},{0,2,1,3},{0,3,1,2}};
        vector<float> dm, am, drmax;
        float min_dm(999);
        for (unsigned ic(0); ic<hcombs.size(); ic++){
          TLorentzVector h1 = sys_ord_jet[hcombs[ic][0]].first + sys_ord_jet[hcombs[ic][1]].first;
          TLorentzVector h2 = sys_ord_jet[hcombs[ic][2]].first + sys_ord_jet[hcombs[ic][3]].first;
          float idm = fabs(h1.M()-h2.M());
          if (idm < min_dm) {
            min_dm = idm;
            pico.out_sys_hig_cand_dm()[ijec] = idm;
            pico.out_sys_hig_cand_am()[ijec] = (h1.M()+h2.M())/2.;
            float idr1 = dR(sys_ord_jet[hcombs[ic][0]].first.Eta(), 
                            sys_ord_jet[hcombs[ic][1]].first.Eta(), 
                            sys_ord_jet[hcombs[ic][0]].first.Phi(), 
                            sys_ord_jet[hcombs[ic][1]].first.Phi());
            float idr2 = dR(sys_ord_jet[hcombs[ic][2]].first.Eta(), 
                            sys_ord_jet[hcombs[ic][3]].first.Eta(), 
                            sys_ord_jet[hcombs[ic][2]].first.Phi(), 
                            sys_ord_jet[hcombs[ic][3]].first.Phi());
            pico.out_sys_hig_cand_drmax()[ijec] = idr1 > idr2 ? idr1 : idr2;
          } //if current dm is the lowest
        } //loop over Higgs pairings
      } //if at least 4 good jets
    } //loop over JEC/JER variations
  } //if (isSignal && !doDeepFlav)

  return;
}
