/*
 * Utility to calculate trigger scale factors for H->Zgamma analysis
 */

#include "trigger_weighter.hpp"

#include "correction.hpp"
#include "pico_tree.hpp"
#include "utilities.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>


TriggerWeighter::TriggerWeighter(int year, bool preVFP) {
  std::string in_file_path;
  std::string in_file_ello;
  std::string in_file_elup;
  std::string in_file_elsi;
  std::string in_file_mulo;
  std::string in_file_muup;
  std::string in_file_musi;
  if (year==2016 && preVFP) {
    in_file_path = "data/zgamma/2016preVFP_UL/";
    in_file_ello = "trigeff_ele12";
    in_file_elup = "trigeff_ele23";
    in_file_elsi = "trigeff_ele27";
    in_file_mulo = "trigeff_mu8";
    in_file_muup = "trigeff_mu17";
    in_file_musi = "trigeff_mu24";
  } 
  else if (year==2016 && !preVFP) {
    in_file_path = "data/zgamma/2016postVFP_UL/";
    in_file_ello = "trigeff_ele12";
    in_file_elup = "trigeff_ele23";
    in_file_elsi = "trigeff_ele27";
    in_file_mulo = "trigeff_mu8";
    in_file_muup = "trigeff_mu17";
    in_file_musi = "trigeff_mu24";
  } 
  else if (year==2017) {
    in_file_path = "data/zgamma/2017_UL/";
    in_file_ello = "trigeff_ele12";
    in_file_elup = "trigeff_ele23";
    in_file_elsi = "trigeff_ele32";
    in_file_mulo = "trigeff_mu8";
    in_file_muup = "trigeff_mu17";
    in_file_musi = "trigeff_mu27";
  } 
  else { //2018
    in_file_path = "data/zgamma/2018_UL/";
    in_file_ello = "trigeff_ele12";
    in_file_elup = "trigeff_ele23";
    in_file_elsi = "trigeff_ele32";
    in_file_mulo = "trigeff_mu8";
    in_file_muup = "trigeff_mu17";
    in_file_musi = "trigeff_mu24";
  } 
  cs_ello_ = correction::CorrectionSet::from_file(in_file_path+in_file_ello+".json");
  cs_elup_ = correction::CorrectionSet::from_file(in_file_path+in_file_elup+".json");
  cs_elsi_ = correction::CorrectionSet::from_file(in_file_path+in_file_elsi+".json");
  cs_mulo_ = correction::CorrectionSet::from_file(in_file_path+in_file_mulo+".json");
  cs_muup_ = correction::CorrectionSet::from_file(in_file_path+in_file_muup+".json");
  cs_musi_ = correction::CorrectionSet::from_file(in_file_path+in_file_musi+".json");
  map_dielectron_lowerleg_ = cs_ello_->at(in_file_ello);
  map_dielectron_upperleg_ = cs_elup_->at(in_file_elup);
  map_single_electron_ = cs_elsi_->at(in_file_elsi);
  map_dimuon_lowerleg_ = cs_mulo_->at(in_file_mulo);
  map_dimuon_upperleg_ = cs_muup_->at(in_file_muup);
  map_single_muon_ = cs_musi_->at(in_file_musi);
}


std::vector<float> TriggerWeighter::GetSF(pico_tree &pico) {
  //this is just a wrapper around the other GetSF that extracts relevant info
  //from pico tree

  std::vector<float> electron_pt;
  std::vector<float> electron_eta;
  std::vector<float> muon_pt;
  std::vector<float> muon_eta;
  for (unsigned int iel = 0; iel < pico.out_el_sig().size(); iel++) {
    if (pico.out_el_sig()[iel]) {
      electron_pt.push_back(pico.out_el_pt()[iel]);
      electron_eta.push_back(pico.out_el_eta()[iel]);
    }
  }
  for (unsigned int imu = 0; imu < pico.out_mu_sig().size(); imu++) {
    if (pico.out_mu_sig()[imu]) {
      muon_pt.push_back(pico.out_mu_pt()[imu]);
      muon_eta.push_back(pico.out_mu_eta()[imu]);
    }
  }

  return GetSF(electron_pt, muon_pt, electron_eta, muon_eta, 
      pico.out_trig_single_el(), pico.out_trig_single_mu(), 
      pico.out_trig_double_el(), pico.out_trig_double_mu());
}


std::vector<float> TriggerWeighter::GetSF(std::vector<float> electron_pt, 
    std::vector<float> muon_pt, std::vector<float> electron_eta, 
    std::vector<float> muon_eta, bool pass_singleel, bool pass_singlemu, 
    bool pass_diel, bool pass_dimu) {
  //note that this only weights leptons that pass the signal criteria
  //i.e. trigger efficiencies will remain uncorrected for leptons failing
  //signal criteria

  if (muon_pt.size()==0 && electron_pt.size()==0) return {1.0,1.0,1.0};

  //get data/mc probability/uncertainty from appropriate functions
  std::vector<float> data_prob = GetTotalProbability(electron_pt, muon_pt, 
    electron_eta, muon_eta, pass_singleel, pass_singlemu, pass_diel, 
    pass_dimu, true);
  std::vector<float> mc_prob = GetTotalProbability(electron_pt, muon_pt, 
    electron_eta, muon_eta, pass_singleel, pass_singlemu, pass_diel, 
    pass_dimu, false);

  //calculate SFs
  if (mc_prob[0] < 0.001f) mc_prob[0] = 0.0;
  float sf = 1.0;
  float unc = 0.0;
  propagate_uncertainty_ratio(data_prob[0], data_prob[1], mc_prob[0], mc_prob[1], sf, unc);
  float sf_up = sf+unc;
  float sf_dn = fmax(sf-unc,0);

  ////deal with signal leptons with low probabilities
  //if (mc_prob[0]<0.001f || data_prob[0]<0.001f) sf = 1.0;
  //if (mc_prob[1]<0.001f || data_prob[1]<0.001f) sf_up = 1.0;
  //if (mc_prob[2]<0.001f || data_prob[2]<0.001f) sf_dn = 1.0;

  return {sf, sf_up, sf_dn};
}


std::vector<float> TriggerWeighter::GetTotalProbability(
    std::vector<float> electron_pt, std::vector<float> muon_pt, 
    std::vector<float> electron_eta, std::vector<float> muon_eta, 
    bool pass_singleel, bool pass_singlemu, bool pass_diel, bool pass_dimu, 
    bool is_data) {

  //get lepton probability/uncertainty from appropriate functions
  std::vector<float> electron_prob = GetFlavorProbability(electron_pt, electron_eta,
    pass_singleel, pass_diel, is_data, true);
  std::vector<float> muon_prob = GetFlavorProbability(muon_pt, muon_eta,
    pass_singlemu, pass_dimu, is_data, false);

  //calculate probability and uncertainty
  //assume electron and muon triggers independent
  float prob = electron_prob[0]+muon_prob[0]-electron_prob[0]*muon_prob[0];
  float unc = hypotf(electron_prob[1],hypotf(muon_prob[1],hypotf(electron_prob[0]*muon_prob[1],muon_prob[0]*electron_prob[1])));
  if (isnan(prob) || isnan(unc)) {
    DBG("Error: NaN propagated in trigger weighter.\n");
    exit(1);
  }
  return {prob, unc};
}


std::vector<float> TriggerWeighter::GetFlavorProbability(
    std::vector<float> lepton_pt, std::vector<float> lepton_eta, 
    bool pass_singlelep, bool pass_dilep, bool is_data, bool is_electron) {

  //The algorithm for calculating probability is somewhat cumbersome due to the
  //correlation between different triggers. The idea is to loop over every
  //possible category where each lepton is assigned one of four exclusive classes
  //  fail all triggers
  //  pass dilepton lower leg but fails upper leg
  //  pass dilepton upper leg but fails single lepton
  //  pass single lepton trigger
  //which should be a complete and exclusive categorization since each of the
  //criteria are increasingly tight. The probability for each of these categories
  //can be directly calculated and the sum of all categories consistent with
  //the event trigger decision gives the total probability.
  
  if (lepton_pt.size()==0)
    return {0.0, 0.0};

  std::vector<LeptonHLTStatus> lepton_status;
  float tot_prob = 0.0;
  float tot_unc = 0.0;
  for (unsigned int ilep = 0; ilep < lepton_pt.size(); ilep++) 
    lepton_status.push_back(LeptonHLTStatus::fail_all);
  //loop over every category
  bool finished = false;
  while (!finished) {

    //check if this exclusive category contributes to the event trigger category
    //if not, skip
    int n_pass_lower = 0;
    int n_pass_upper = 0;
    int n_pass_single = 0;
    for (unsigned int ilep = 0; ilep < lepton_status.size(); ilep++) {
      if (lepton_status[ilep]>=LeptonHLTStatus::pass_lowerdilep) n_pass_lower++;
      if (lepton_status[ilep]>=LeptonHLTStatus::pass_upperdilep) n_pass_upper++;
      if (lepton_status[ilep]==LeptonHLTStatus::pass_singlelep) n_pass_single++;
    }
    bool relevant_cat = true;
    if ((n_pass_single>0)!=pass_singlelep) relevant_cat = false;
    if ((n_pass_upper>0&&n_pass_lower>1)!=pass_dilep) relevant_cat = false;

    if (relevant_cat) {
      //the probability for an exclusive category is just the product of the 
      //probabilities for each lepton in the category
      float cat_prob = 1.0;
      float cat_unc = 0.0;
      for (unsigned int ilep = 0; ilep < lepton_status.size(); ilep++) {
        float lep_prob = 0.0;
        float lep_unc = 0.0;

        //probability to fail all is 1 - probability to pass dilep lower leg
        if (lepton_status[ilep]==LeptonHLTStatus::fail_all) {
          std::vector<float> prob_lower = GetLeptonProbability(lepton_pt[ilep],
              lepton_eta[ilep],is_data,is_electron,
              LeptonHLTStatus::pass_lowerdilep);
          lep_prob = 1.0-prob_lower[0];
          lep_unc = prob_lower[1];
        }

        //probability to pass only lower leg is prob(lower leg)-prob(upper leg)
        else if (lepton_status[ilep]==LeptonHLTStatus::pass_lowerdilep) {
          std::vector<float> prob_lower = GetLeptonProbability(lepton_pt[ilep],
              lepton_eta[ilep],is_data,is_electron,
              LeptonHLTStatus::pass_lowerdilep);
          std::vector<float> prob_upper = GetLeptonProbability(lepton_pt[ilep],
              lepton_eta[ilep],is_data,is_electron,
              LeptonHLTStatus::pass_upperdilep);
          lep_prob = prob_lower[0]-prob_upper[0];
          lep_unc = hypotf(prob_lower[1], prob_upper[1]);
        }

        //probability to pass upper leg but fail single lep is 
        //prob(upper leg)-prob(single lep)
        else if (lepton_status[ilep]==LeptonHLTStatus::pass_upperdilep) {
          std::vector<float> prob_upper = GetLeptonProbability(lepton_pt[ilep],
              lepton_eta[ilep],is_data,is_electron,
              LeptonHLTStatus::pass_upperdilep);
          std::vector<float> prob_single = GetLeptonProbability(lepton_pt[ilep],
              lepton_eta[ilep],is_data,is_electron,
              LeptonHLTStatus::pass_singlelep);
          lep_prob = prob_upper[0]-prob_single[0];
          lep_unc = hypotf(prob_upper[1], prob_single[1]);
        }

        //probability to pass single lep trigger
        else {
          std::vector<float> prob_single = GetLeptonProbability(lepton_pt[ilep],
              lepton_eta[ilep],is_data,is_electron,
              LeptonHLTStatus::pass_singlelep);
          lep_prob = prob_single[0];
          lep_unc = prob_single[1];
        }

        //deal with rounding errors/low stats categories
        if (lep_prob < 0) lep_prob = 0;

        //combine lepton probability into total
        cat_unc = hypotf(lep_unc*cat_prob,cat_unc*lep_prob);
        cat_prob *= lep_prob;
      }

      //add to total probability
      tot_prob += cat_prob;
      tot_unc = hypotf(tot_unc,cat_unc);
    }

    //iterate to next category
    for (unsigned int ilep = 0; ilep < lepton_status.size(); ilep++) {
      if (lepton_status[ilep] == LeptonHLTStatus::fail_all) {
        lepton_status[ilep] = LeptonHLTStatus::pass_lowerdilep;
        break;
      }
      else if (lepton_status[ilep] == LeptonHLTStatus::pass_lowerdilep) {
        lepton_status[ilep] = LeptonHLTStatus::pass_upperdilep;
        break;
      }
      else if (lepton_status[ilep] == LeptonHLTStatus::pass_upperdilep) {
        lepton_status[ilep] = LeptonHLTStatus::pass_singlelep;
        break;
      }
      else {
        lepton_status[ilep] = LeptonHLTStatus::fail_all;
        if (ilep == lepton_status.size()-1)
          finished = true;
      }
    }
  }

  return {tot_prob, tot_unc};
}


std::vector<float> TriggerWeighter::GetLeptonProbability(float lepton_pt, float lepton_eta,
    bool is_data, bool is_electron, LeptonHLTStatus trigger_leg) {

  //select map and key based on arguments
  correction::Correction::Ref* prob_map;
  if (is_electron) {
    if (trigger_leg == LeptonHLTStatus::pass_lowerdilep) {
      prob_map = &map_dielectron_lowerleg_;
    }
    else if (trigger_leg == LeptonHLTStatus::pass_upperdilep) {
      prob_map = &map_dielectron_upperleg_;
    }
    else {
      prob_map = &map_single_electron_;
    }
  }
  else {
    if (trigger_leg == LeptonHLTStatus::pass_lowerdilep) {
      prob_map = &map_dimuon_lowerleg_;
    }
    else if (trigger_leg == LeptonHLTStatus::pass_upperdilep) {
      prob_map = &map_dimuon_upperleg_;
    }
    else {
      prob_map = &map_single_muon_;
    }
  }
  std::string eff_name = "effdata";
  std::string unc_name = "systdata";
  if (!is_data) {
    eff_name = "effmc";
    unc_name = "systmc";
  }

  //evaluate
  float prob = (*prob_map)->evaluate({eff_name, std::abs(lepton_eta), lepton_pt});
  float uncr = (*prob_map)->evaluate({unc_name, std::abs(lepton_eta), lepton_pt});
  return {prob, uncr};
}



