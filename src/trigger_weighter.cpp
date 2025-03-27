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

using std::cout;
using std::endl;
using std::string;
using std::vector;

float bound(float value, float upper, float lower) {
  if (value>upper)
    return upper;
  if (value<lower)
    return lower;
  return value;
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


vector<float> TriggerWeighter::GetSF(pico_tree &pico) {
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

  return GetSF(electron_pt, muon_pt, electron_eta, muon_eta, electron_phi,
      pico.out_trig_single_el(), pico.out_trig_single_mu(), 
      pico.out_trig_double_el(), pico.out_trig_double_mu());
}


vector<float> TriggerWeighter::GetSF(std::vector<float> electron_pt, 
    vector<float> muon_pt, std::vector<float> electron_eta, 
    std::vector<float> muon_eta, std::vector<float> electron_phi,
    bool pass_singleel, bool pass_singlemu, 
    bool pass_diel, bool pass_dimu) {
  //note that this only weights leptons that pass the signal criteria
  //i.e. trigger efficiencies will remain uncorrected for leptons failing
  //signal criteria

  if (muon_pt.size()==0 && electron_pt.size()==0) return {1.0,1.0,1.0};

  //get data/mc probability/uncertainty from appropriate functions
  vector<float> data_prob = GetTotalProbability(electron_pt, muon_pt, 
    electron_eta, muon_eta, electron_phi, pass_singleel, pass_singlemu, 
    pass_diel, 
    pass_dimu, true);
  vector<float> mc_prob = GetTotalProbability(electron_pt, muon_pt, 
    electron_eta, muon_eta, electron_phi, pass_singleel, pass_singlemu, 
    pass_diel, pass_dimu, false);

  //calculate SFs
  if (mc_prob[0] < 0.001f) mc_prob[0] = 0.0;
  float sf = 1.0f;
  float unc = 0.0f;
  propagate_uncertainty_ratio(data_prob[0], data_prob[1], mc_prob[0], 
                              mc_prob[1], sf, unc);
  bool pass_trig = (pass_singleel || pass_singlemu || pass_diel || pass_dimu);
  float sf_up = sf+unc;
  float sf_dn = sf-unc;
  if (!pass_trig) {
    sf_up = sf-unc;
    sf_dn = sf+unc;
  }
  sf = bound(sf,5.0,0.0);
  sf_up = bound(sf_up,5.0,0.0);
  sf_dn = bound(sf_dn,5.0,0.0);

  return {sf, sf_up, sf_dn};
}


vector<float> TriggerWeighter::GetTotalProbability(
    vector<float> electron_pt, vector<float> muon_pt, 
    vector<float> electron_eta, vector<float> muon_eta, 
    vector<float> electron_phi,
    bool pass_singleel, bool pass_singlemu, bool pass_diel, bool pass_dimu, 
    bool is_data) {

  //get lepton probability/uncertainty from appropriate functions
  vector<float> muon_phi(muon_pt.size(), 0.0);
  vector<float> electron_prob = GetFlavorProbability(electron_pt, electron_eta,
    electron_phi, pass_singleel, pass_diel, is_data, true);
  vector<float> muon_prob = GetFlavorProbability(muon_pt, muon_eta, muon_phi,
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
          vector<float> prob_lower = GetLeptonProbability(lepton_pt[ilep],
              lepton_eta[ilep],lepton_phi[ilep],is_data,is_electron,
              LeptonHLTStatus::pass_lowerdilep);
          lep_prob = 1.0-prob_lower[0];
          lep_unc = prob_lower[1];
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
          lep_unc = hypotf(prob_lower[1], prob_upper[1]);
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
          lep_unc = hypotf(prob_upper[1], prob_single[1]);
        }

        //probability to pass single lep trigger
        else {
          vector<float> prob_single = GetLeptonProbability(lepton_pt[ilep],
              lepton_eta[ilep],lepton_phi[ilep],is_data,is_electron,
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
  return {prob, uncr};
}



