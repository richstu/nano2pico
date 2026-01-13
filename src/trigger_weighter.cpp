/*
 * Utility to calculate trigger scale factors for H->Zgamma analysis
 */

#include "trigger_weighter.hpp"

#include "correction.h"
#include "pico_tree.hpp"
#include "utilities.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::isnan;

float bound(float value, float upper, float lower) {
  if (value>upper)
    return upper;
  if (value<lower)
    return lower;
  return value;
}

float safe_div(float num, float den) {
  if (den > 0.0f)
    return num/den;
  return 1.0f;
}

TriggerWeighter::TriggerWeighter(string year) {
  string in_file_path;
  string in_file_ello;
  string in_file_elup;
  string in_file_elsi;
  string in_file_mulo;
  string in_file_muup;
  string in_file_musi;
  if (year=="2016APV") {
    in_file_path = "data/zgamma/2016preVFP_UL/";
    in_file_ello = "hzg_eltrig12_2016APV_efficiencies";
    in_file_elup = "hzg_eltrig23_2016APV_efficiencies";
    in_file_elsi = "hzg_eltrig27_2016APV_efficiencies";
    in_file_mulo = "hzg_mutrig8_2016APV_efficiencies";
    in_file_muup = "hzg_mutrig17_2016APV_efficiencies";
    in_file_musi = "hzg_mutrig24_2016APV_efficiencies";
  } 
  else if (year=="2016") {
    in_file_path = "data/zgamma/2016postVFP_UL/";
    in_file_ello = "hzg_eltrig12_2016_efficiencies";
    in_file_elup = "hzg_eltrig23_2016_efficiencies";
    in_file_elsi = "hzg_eltrig27_2016_efficiencies";
    in_file_mulo = "hzg_mutrig8_2016_efficiencies";
    in_file_muup = "hzg_mutrig17_2016_efficiencies";
    in_file_musi = "hzg_mutrig24_2016_efficiencies";
  } 
  else if (year=="2017") {
    in_file_path = "data/zgamma/2017_UL/";
    in_file_ello = "hzg_eltrig12_2017_efficiencies";
    in_file_elup = "hzg_eltrig23_2017_efficiencies";
    in_file_elsi = "hzg_eltrig32_2017_efficiencies";
    in_file_mulo = "hzg_mutrig8_2017_efficiencies";
    in_file_muup = "hzg_mutrig17_2017_efficiencies";
    in_file_musi = "hzg_mutrig27_2017_efficiencies";
  } 
  else if (year=="2018") {
    in_file_path = "data/zgamma/2018_UL/";
    in_file_ello = "hzg_eltrig12_2018_efficiencies";
    in_file_elup = "hzg_eltrig23_2018_efficiencies";
    in_file_elsi = "hzg_eltrig32_2018_efficiencies";
    in_file_mulo = "hzg_mutrig8_2018_efficiencies";
    in_file_muup = "hzg_mutrig17_2018_efficiencies";
    in_file_musi = "hzg_mutrig24_2018_efficiencies";
  }
  else if (year=="2022") {
    in_file_path = "data/zgamma/2022/";
    in_file_ello = "hzg_eltrig12_2022_efficiencies";
    in_file_elup = "hzg_eltrig23_2022_efficiencies";
    in_file_elsi = "hzg_eltrig30_2022_efficiencies";
    in_file_mulo = "hzg_mutrig8_2022_efficiencies";
    in_file_muup = "hzg_mutrig17_2022_efficiencies";
    in_file_musi = "hzg_mutrig24_2022_efficiencies";
  }
  else if (year=="2022EE") {
    in_file_path = "data/zgamma/2022EE/";
    in_file_ello = "hzg_eltrig12_2022EE_efficiencies";
    in_file_elup = "hzg_eltrig23_2022EE_efficiencies";
    in_file_elsi = "hzg_eltrig30_2022EE_efficiencies";
    in_file_mulo = "hzg_mutrig8_2022EE_efficiencies";
    in_file_muup = "hzg_mutrig17_2022EE_efficiencies";
    in_file_musi = "hzg_mutrig24_2022EE_efficiencies";
  }
  else if (year=="2023") {
    in_file_path = "data/zgamma/2023/";
    in_file_ello = "hzg_eltrig12_2023_efficiencies";
    in_file_elup = "hzg_eltrig23_2023_efficiencies";
    in_file_elsi = "hzg_eltrig30_2023_efficiencies";
    in_file_mulo = "hzg_mutrig8_2023_efficiencies";
    in_file_muup = "hzg_mutrig17_2023_efficiencies";
    in_file_musi = "hzg_mutrig24_2023_efficiencies";
  }
  else if (year=="2023BPix") {
    in_file_path = "data/zgamma/2023BPix/";
    in_file_ello = "hzg_eltrig12_2023BPix_efficiencies";
    in_file_elup = "hzg_eltrig23_2023BPix_efficiencies";
    in_file_elsi = "hzg_eltrig30_2023BPix_efficiencies";
    in_file_mulo = "hzg_mutrig8_2023BPix_efficiencies";
    in_file_muup = "hzg_mutrig17_2023BPix_efficiencies";
    in_file_musi = "hzg_mutrig24_2023BPix_efficiencies";
    cs_ello_hole_ = correction::CorrectionSet::from_file(in_file_path
        +"hzg_eltrig12_2023BPixHole_efficiencies.json");
    cs_elup_hole_ = correction::CorrectionSet::from_file(in_file_path
        +"hzg_eltrig23_2023BPixHole_efficiencies.json");
    cs_elsi_hole_ = correction::CorrectionSet::from_file(in_file_path
        +"hzg_eltrig30_2023BPixHole_efficiencies.json");
  }
  else { //2018
    cout << "WARNING: No trigger weights, defaulting to 2018" << endl;
    in_file_path = "data/zgamma/2018_UL/";
    in_file_ello = "hzg_eltrig12_2018_efficiencies";
    in_file_elup = "hzg_eltrig23_2018_efficiencies";
    in_file_elsi = "hzg_eltrig32_2018_efficiencies";
    in_file_mulo = "hzg_mutrig8_2018_efficiencies";
    in_file_muup = "hzg_mutrig17_2018_efficiencies";
    in_file_musi = "hzg_mutrig24_2018_efficiencies";
  } 
  post_bpix_ = false;
  cs_ello_ = correction::CorrectionSet::from_file(in_file_path+in_file_ello
                                                  +".json");
  cs_elup_ = correction::CorrectionSet::from_file(in_file_path+in_file_elup
                                                  +".json");
  cs_elsi_ = correction::CorrectionSet::from_file(in_file_path+in_file_elsi
                                                  +".json");
  cs_mulo_ = correction::CorrectionSet::from_file(in_file_path+in_file_mulo
                                                  +".json");
  cs_muup_ = correction::CorrectionSet::from_file(in_file_path+in_file_muup
                                                  +".json");
  cs_musi_ = correction::CorrectionSet::from_file(in_file_path+in_file_musi
                                                  +".json");
  map_diel_lower_dataeff_ = cs_ello_->at("effdata");
  map_diel_lower_dataunc_ = cs_ello_->at("systdata");
  map_diel_lower_mceff_ = cs_ello_->at("effmc");
  map_diel_lower_mcunc_ = cs_ello_->at("systmc");
  map_diel_upper_dataeff_ = cs_elup_->at("effdata");
  map_diel_upper_dataunc_ = cs_elup_->at("systdata");
  map_diel_upper_mceff_ = cs_elup_->at("effmc");
  map_diel_upper_mcunc_ = cs_elup_->at("systmc");
  map_singleel_dataeff_ = cs_elsi_->at("effdata");
  map_singleel_dataunc_ = cs_elsi_->at("systdata");
  map_singleel_mceff_ = cs_elsi_->at("effmc");
  map_singleel_mcunc_ = cs_elsi_->at("systmc");
  map_dimu_lower_dataeff_ = cs_mulo_->at("effdata");
  map_dimu_lower_dataunc_ = cs_mulo_->at("systdata");
  map_dimu_lower_mceff_ = cs_mulo_->at("effmc");
  map_dimu_lower_mcunc_ = cs_mulo_->at("systmc");
  map_dimu_upper_dataeff_ = cs_muup_->at("effdata");
  map_dimu_upper_dataunc_ = cs_muup_->at("systdata");
  map_dimu_upper_mceff_ = cs_muup_->at("effmc");
  map_dimu_upper_mcunc_ = cs_muup_->at("systmc");
  map_singlemu_dataeff_ = cs_musi_->at("effdata");
  map_singlemu_dataunc_ = cs_musi_->at("systdata");
  map_singlemu_mceff_ = cs_musi_->at("effmc");
  map_singlemu_mcunc_ = cs_musi_->at("systmc");
  if (year == "2023BPix") {
    map_hole_diel_lower_dataeff_ = cs_ello_hole_->at("effdata");
    map_hole_diel_lower_dataunc_ = cs_ello_hole_->at("systdata");
    map_hole_diel_lower_mceff_   = cs_ello_hole_->at("effmc");
    map_hole_diel_lower_mcunc_   = cs_ello_hole_->at("systmc");
    map_hole_diel_upper_dataeff_ = cs_elup_hole_->at("effdata");
    map_hole_diel_upper_dataunc_ = cs_elup_hole_->at("systdata");
    map_hole_diel_upper_mceff_   = cs_elup_hole_->at("effmc");
    map_hole_diel_upper_mcunc_   = cs_elup_hole_->at("systmc");
    map_hole_singleel_dataeff_   = cs_elsi_hole_->at("effdata");
    map_hole_singleel_dataunc_   = cs_elsi_hole_->at("systdata");
    map_hole_singleel_mceff_     = cs_elsi_hole_->at("effmc");
    map_hole_singleel_mcunc_     = cs_elsi_hole_->at("systmc");
    post_bpix_ = true;
  }
}


void TriggerWeighter::GetSF(pico_tree &pico) {
  //this is just a wrapper around the other GetSF that extracts relevant info
  //from pico tree

  vector<float> electron_pt;
  vector<float> electron_eta;
  vector<float> electron_phi;
  vector<float> muon_pt;
  vector<float> muon_eta;
  for (unsigned int iel = 0; iel < pico.out_el_sig().size(); iel++) {
    if (pico.out_el_sig()[iel]) {
      electron_pt.push_back(pico.out_el_pt()[iel]);
      electron_eta.push_back(pico.out_el_eta()[iel]);
      electron_phi.push_back(pico.out_el_phi()[iel]);
    }
  }
  for (unsigned int imu = 0; imu < pico.out_mu_sig().size(); imu++) {
    if (pico.out_mu_sig()[imu]) {
      muon_pt.push_back(pico.out_mu_pt()[imu]);
      muon_eta.push_back(pico.out_mu_eta()[imu]);
    }
  }

  vector<float> sfs =  GetSF(electron_pt, muon_pt, electron_eta, muon_eta, 
      electron_phi, pico.out_trig_single_el(), pico.out_trig_single_mu(), 
      pico.out_trig_double_el(), pico.out_trig_double_mu());
  pico.out_w_trig() = sfs[0];
  pico.out_sys_trig().resize(2,1.0f); //unused in ZGamma
  pico.out_sys_trig_el().resize(2,1.0f);
  pico.out_sys_trig_mu().resize(2,1.0f);
  pico.out_sys_trig_el()[0] = sfs[1];
  pico.out_sys_trig_el()[1] = sfs[2];
  pico.out_sys_trig_mu()[0] = sfs[3];
  pico.out_sys_trig_mu()[1] = sfs[4];
}

vector<float> TriggerWeighter::GetSF(
    vector<float> electron_pt, vector<float> muon_pt, 
    vector<float> electron_eta, vector<float> muon_eta, 
    vector<float> electron_phi, bool pass_singleel, bool pass_singlemu, 
    bool pass_diel, bool pass_dimu) {

  //note that this only weights leptons that pass the signal criteria
  //i.e. trigger efficiencies will remain uncorrected for leptons failing
  //signal criteria
  //get lepton probability/uncertainty from appropriate functions
  vector<float> el_prob_data, mu_prob_data, el_prob_mc, mu_prob_mc;
  vector<float> muon_phi(muon_pt.size(), 0.0);
  if (electron_pt.size()==0) {
    el_prob_data.resize(3,1.0f);
    el_prob_mc.resize(3,1.0f);
  }
  else {
    el_prob_data = GetFlavorProbability(electron_pt, 
        electron_eta, electron_phi, pass_singleel, pass_diel, true, true);
    el_prob_mc = GetFlavorProbability(electron_pt, 
        electron_eta, electron_phi, pass_singleel, pass_diel, false, true);
  }
  if (muon_pt.size()==0) {
    mu_prob_data.resize(3,1.0f);
    mu_prob_mc.resize(3,1.0f);
  }
  else {
    mu_prob_data = GetFlavorProbability(muon_pt, muon_eta, 
        muon_phi, pass_singlemu, pass_dimu, true, false);
    mu_prob_mc = GetFlavorProbability(muon_pt, muon_eta, 
        muon_phi, pass_singlemu, pass_dimu, false, false);
  }

  //calculate SFs
  //assume worst case for variations (data overestimated, MC under or v.v.)
  float sflep_el_nm = safe_div(el_prob_data[0], el_prob_mc[0]);
  float sflep_el_up = safe_div(el_prob_data[1], el_prob_mc[2]);
  float sflep_el_dn = safe_div(el_prob_data[2], el_prob_mc[1]);
  float sflep_mu_nm = safe_div(mu_prob_data[0], mu_prob_mc[0]);
  float sflep_mu_up = safe_div(mu_prob_data[1], mu_prob_mc[2]);
  float sflep_mu_dn = safe_div(mu_prob_data[2], mu_prob_mc[1]);
  float sf = sflep_el_nm*sflep_mu_nm;
  float sf_elup = sflep_el_up*sflep_mu_nm;
  float sf_eldn = sflep_el_dn*sflep_mu_nm;
  float sf_muup = sflep_el_nm*sflep_mu_up;
  float sf_mudn = sflep_el_nm*sflep_mu_dn;
  if (isnan(sf) || isnan(sf_elup) || isnan(sf_eldn) || isnan(sf_muup) 
      || isnan(sf_mudn)) {
    DBG("Error: NaN propagated in trigger weighter.\n");
    exit(1);
  }
  sf = bound(sf, 10.0f, 0.0f);
  sf_elup = bound(sf_elup, 10.0f, 0.0f);
  sf_eldn = bound(sf_eldn, 10.0f, 0.0f);
  sf_muup = bound(sf_muup, 10.0f, 0.0f);
  sf_mudn = bound(sf_mudn, 10.0f, 0.0f);

  return {sf, sf_elup, sf_eldn, sf_muup, sf_mudn};
}


vector<float> TriggerWeighter::GetFlavorProbability(
    vector<float> lepton_pt, vector<float> lepton_eta, 
    vector<float> lepton_phi, bool pass_singlelep, 
    bool pass_dilep, bool is_data, bool is_electron) {

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

  vector<LeptonHLTStatus> lepton_status;
  float tot_prob = 0.0;
  float tot_prob_up = 0.0;
  float tot_prob_dn = 0.0;
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
      float cat_prob_up = 1.0;
      float cat_prob_dn = 1.0;
      for (unsigned int ilep = 0; ilep < lepton_status.size(); ilep++) {
        float lep_prob = 0.0;
        float lep_prob_up = 0.0;
        float lep_prob_dn = 0.0;

        //probability to fail all is 1 - probability to pass dilep lower leg
        if (lepton_status[ilep]==LeptonHLTStatus::fail_all) {
          vector<float> prob_lower = GetLeptonProbability(lepton_pt[ilep],
              lepton_eta[ilep],lepton_phi[ilep],is_data,is_electron,
              LeptonHLTStatus::pass_lowerdilep);
          lep_prob = 1.0-prob_lower[0];
          lep_prob_up = 1.0-prob_lower[1];
          lep_prob_dn = 1.0-prob_lower[2];
        }

        //probability to pass only lower leg is prob(lower leg)-prob(upper leg)
        else if (lepton_status[ilep]==LeptonHLTStatus::pass_lowerdilep) {
          vector<float> prob_lower = GetLeptonProbability(lepton_pt[ilep],
              lepton_eta[ilep],lepton_phi[ilep],is_data,is_electron,
              LeptonHLTStatus::pass_lowerdilep);
          vector<float> prob_upper = GetLeptonProbability(lepton_pt[ilep],
              lepton_eta[ilep],lepton_phi[ilep],is_data,is_electron,
              LeptonHLTStatus::pass_upperdilep);
          lep_prob = prob_lower[0]-prob_upper[0];
          lep_prob_up = prob_lower[1]-prob_upper[1];
          lep_prob_dn = prob_lower[2]-prob_upper[2];
        }

        //probability to pass upper leg but fail single lep is 
        //prob(upper leg)-prob(single lep)
        else if (lepton_status[ilep]==LeptonHLTStatus::pass_upperdilep) {
          vector<float> prob_upper = GetLeptonProbability(lepton_pt[ilep],
              lepton_eta[ilep],lepton_phi[ilep],is_data,is_electron,
              LeptonHLTStatus::pass_upperdilep);
          vector<float> prob_single = GetLeptonProbability(lepton_pt[ilep],
              lepton_eta[ilep],lepton_phi[ilep],is_data,is_electron,
              LeptonHLTStatus::pass_singlelep);
          lep_prob = prob_upper[0]-prob_single[0];
          lep_prob_up = prob_upper[1]-prob_single[1];
          lep_prob_dn = prob_upper[2]-prob_single[2];
        }

        //probability to pass single lep trigger
        else {
          vector<float> prob_single = GetLeptonProbability(lepton_pt[ilep],
              lepton_eta[ilep],lepton_phi[ilep],is_data,is_electron,
              LeptonHLTStatus::pass_singlelep);
          lep_prob = prob_single[0];
          lep_prob_up = prob_single[1];
          lep_prob_dn = prob_single[2];
        }

        //deal with rounding errors/low stats categories
        if (lep_prob < 0) lep_prob = 0;
        if (lep_prob_up < 0) lep_prob_up = 0;
        if (lep_prob_dn < 0) lep_prob_dn = 0;

        //combine lepton probability into total
        cat_prob *= lep_prob;
        cat_prob_up *= lep_prob_up;
        cat_prob_dn *= lep_prob_dn;
      }

      //add to total probability
      tot_prob += cat_prob;
      tot_prob_up += cat_prob_up;
      tot_prob_dn += cat_prob_dn;
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

  return {tot_prob, tot_prob_up, tot_prob_dn};
}


vector<float> TriggerWeighter::GetLeptonProbability(float lepton_pt, 
    float lepton_eta, float lepton_phi, bool is_data, bool is_electron, 
    LeptonHLTStatus trigger_leg) {

  //select map and key based on arguments
  correction::Correction::Ref* prob_map;
  correction::Correction::Ref* unc_map;
  if (is_electron) {
    bool in_bpix_region = (lepton_eta > -1.566 && lepton_eta < 0.0 
                           && lepton_phi > -1.2 && lepton_phi < -0.8);
    if (trigger_leg == LeptonHLTStatus::pass_lowerdilep) {
      if (is_data) {
        if (in_bpix_region && post_bpix_) {
          prob_map = &map_hole_diel_lower_dataeff_;
          unc_map = &map_hole_diel_lower_dataunc_;
        }
        else {
          prob_map = &map_diel_lower_dataeff_;
          unc_map = &map_diel_lower_dataunc_;
        }
      }
      else {
        if (in_bpix_region && post_bpix_) {
          prob_map = &map_hole_diel_lower_mceff_;
          unc_map = &map_hole_diel_lower_mcunc_;
        }
        else {
          prob_map = &map_diel_lower_mceff_;
          unc_map = &map_diel_lower_mcunc_;
        }
      }
    }
    else if (trigger_leg == LeptonHLTStatus::pass_upperdilep) {
      if (is_data) {
        if (in_bpix_region && post_bpix_) {
          prob_map = &map_hole_diel_upper_dataeff_;
          unc_map = &map_hole_diel_upper_dataunc_;
        }
        else {
          prob_map = &map_diel_upper_dataeff_;
          unc_map = &map_diel_upper_dataunc_;
        }
      }
      else {
        if (in_bpix_region && post_bpix_) {
          prob_map = &map_hole_diel_upper_mceff_;
          unc_map = &map_hole_diel_upper_mcunc_;
        }
        else {
          prob_map = &map_diel_upper_mceff_;
          unc_map = &map_diel_upper_mcunc_;
        }
      }
    }
    else {
      if (is_data) {
        if (in_bpix_region && post_bpix_) {
          prob_map = &map_hole_singleel_dataeff_;
          unc_map = &map_hole_singleel_dataunc_;
        }
        else {
          prob_map = &map_singleel_dataeff_;
          unc_map = &map_singleel_dataunc_;
        }
      }
      else {
        if (in_bpix_region && post_bpix_) {
          prob_map = &map_hole_singleel_mceff_;
          unc_map = &map_hole_singleel_mcunc_;
        }
        else {
          prob_map = &map_singleel_mceff_;
          unc_map = &map_singleel_mcunc_;
        }
      }
    }
  }
  else {
    if (trigger_leg == LeptonHLTStatus::pass_lowerdilep) {
      if (is_data) {
        prob_map = &map_dimu_lower_dataeff_;
        unc_map = &map_dimu_lower_dataunc_;
      }
      else {
        prob_map = &map_dimu_lower_mceff_;
        unc_map = &map_dimu_lower_mcunc_;
      }
    }
    else if (trigger_leg == LeptonHLTStatus::pass_upperdilep) {
      if (is_data) {
        prob_map = &map_dimu_upper_dataeff_;
        unc_map = &map_dimu_upper_dataunc_;
      }
      else {
        prob_map = &map_dimu_upper_mceff_;
        unc_map = &map_dimu_upper_mcunc_;
      }
    }
    else {
      if (is_data) {
        prob_map = &map_singlemu_dataeff_;
        unc_map = &map_singlemu_dataunc_;
      }
      else {
        prob_map = &map_singlemu_mceff_;
        unc_map = &map_singlemu_mcunc_;
      }
    }
  }
  float map_eta = lepton_eta;
  if (!is_electron) map_eta = fabs(lepton_eta);

  //evaluate
  float prob = (*prob_map)->evaluate({lepton_pt, map_eta});
  float uncr = (*unc_map)->evaluate({lepton_pt, map_eta});
  return {prob, prob+uncr, prob-uncr};
}



