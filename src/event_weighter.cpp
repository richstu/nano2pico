#include "event_weighter.hpp"

#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "TFile.h"
#include "TGraphAsymmErrors.h"

#include "utilities.hpp"

using namespace std;

EventWeighter::EventWeighter(int year, bool preVFP){
  if (year==2016 && preVFP) {
    in_file_electron_        = "data/zgamma/2016preVFP_UL/electron_WPL.json";
    in_file_photon_          = "data/zgamma/2016preVFP_UL/photon.json";
    in_file_muon_            = "data/zgamma/2016preVFP_UL/muon_Z.json";
    in_file_muon_lowpt_reco_ = "data/zgamma/2016preVFP_UL/muon_Jpsi_reco.json";
    in_file_muon_lowpt_id_   = "data/zgamma/2016preVFP_UL/muon_Jpsi_id.json";
    in_file_muon_mceff_      = "data/zgamma/2016preVFP_UL/muon_mceff.json";
    in_file_pu_              = "data/zgamma/2016preVFP_UL/puWeights.json";
    in_file_btag_            = "data/zgamma/2016preVFP_UL/btagging.json";
    key_                     = "2016preVFP";
    puName_                  = "Collisions16_UltraLegacy_goldenJSON";
  } else if (year==2016) {
    in_file_electron_        = "data/zgamma/2016postVFP_UL/electron_WPL.json";
    in_file_photon_          = "data/zgamma/2016postVFP_UL/photon.json";
    in_file_muon_            = "data/zgamma/2016postVFP_UL/muon_Z.json";
    in_file_muon_lowpt_reco_ = "data/zgamma/2016postVFP_UL/muon_Jpsi_reco.json";
    in_file_muon_lowpt_id_   = "data/zgamma/2016postVFP_UL/muon_Jpsi_id.json";
    in_file_muon_mceff_      = "data/zgamma/2016postVFP_UL/muon_mceff.json";
    in_file_pu_              = "data/zgamma/2016postVFP_UL/puWeights.json";
    in_file_btag_            = "data/zgamma/2016postVFP_UL/btagging.json";
    key_                     = "2016postVFP";
    puName_                  = "Collisions16_UltraLegacy_goldenJSON";
  } else if (year==2017) {
    in_file_electron_        = "data/zgamma/2017_UL/electron_WPL.json";
    in_file_photon_          = "data/zgamma/2017_UL/photon.json";
    in_file_muon_            = "data/zgamma/2017_UL/muon_Z.json";
    in_file_muon_lowpt_reco_ = "data/zgamma/2017_UL/muon_Jpsi_reco.json";
    in_file_muon_lowpt_id_   = "data/zgamma/2017_UL/muon_Jpsi_id.json";
    in_file_muon_mceff_      = "data/zgamma/2017_UL/muon_mceff.json";
    in_file_pu_              = "data/zgamma/2017_UL/puWeights.json";
    in_file_btag_            = "data/zgamma/2017_UL/btagging.json";
    key_                     = "2017";
    puName_                  = "Collisions17_UltraLegacy_goldenJSON";
  } else if (year==2018) {
    in_file_electron_        = "data/zgamma/2018_UL/electron_WPL.json";
    in_file_photon_          = "data/zgamma/2018_UL/photon.json";
    in_file_muon_            = "data/zgamma/2018_UL/muon_Z.json";
    in_file_muon_lowpt_reco_ = "data/zgamma/2018_UL/muon_Jpsi_reco.json";
    in_file_muon_lowpt_id_   = "data/zgamma/2018_UL/muon_Jpsi_id.json";
    in_file_muon_mceff_      = "data/zgamma/2018_UL/muon_mceff.json";
    in_file_pu_              = "data/zgamma/2018_UL/puWeights.json";
    in_file_btag_            = "data/zgamma/2018_UL/btagging.json";
    key_                     = "2018";
    puName_                  = "Collisions18_UltraLegacy_goldenJSON";
  } else if (year==2022){
    std::cout<<"Using 2018 JSONs by default for now in event_weighter.cpp"<<endl;
    in_file_electron_        = "data/zgamma/2018_UL/electron_WPL.json";
    in_file_photon_          = "data/zgamma/2018_UL/photon.json";
    in_file_muon_            = "data/zgamma/2018_UL/muon_Z.json";
    in_file_muon_lowpt_reco_ = "data/zgamma/2018_UL/muon_Jpsi_reco.json";
    in_file_muon_lowpt_id_   = "data/zgamma/2018_UL/muon_Jpsi_id.json";
    in_file_muon_mceff_      = "data/zgamma/2018_UL/muon_mceff.json";
    in_file_pu_              = "data/zgamma/2018_UL/puWeights.json";
    in_file_btag_            = "data/zgamma/2018_UL/btagging.json";
    key_                     = "2018";
    puName_                  = "Collisions18_UltraLegacy_goldenJSON";
  } else if (year==2023){
    std::cout<<"Using 2018 JSONs by default for now in event_weighter.cpp"<<endl;
    in_file_electron_        = "data/zgamma/2018_UL/electron_WPL.json";
    in_file_photon_          = "data/zgamma/2018_UL/photon.json";
    in_file_muon_            = "data/zgamma/2018_UL/muon_Z.json";
    in_file_muon_lowpt_reco_ = "data/zgamma/2018_UL/muon_Jpsi_reco.json";
    in_file_muon_lowpt_id_   = "data/zgamma/2018_UL/muon_Jpsi_id.json";
    in_file_muon_mceff_      = "data/zgamma/2018_UL/muon_mceff.json";
    in_file_pu_              = "data/zgamma/2018_UL/puWeights.json";
    in_file_btag_            = "data/zgamma/2018_UL/btagging.json";
    key_                     = "2018";
    puName_                  = "Collisions18_UltraLegacy_goldenJSON";
  } else {
    std::cout<<"Year has not been implemented in event_weighter"<<endl;
  }
  cs_electron_         = correction::CorrectionSet::from_file(in_file_electron_);
  cs_photon_           = correction::CorrectionSet::from_file(in_file_photon_);
  cs_muon_             = correction::CorrectionSet::from_file(in_file_muon_);
  cs_muon_lowpt_reco_  = correction::CorrectionSet::from_file(in_file_muon_lowpt_reco_);
  cs_muon_lowpt_id_    = correction::CorrectionSet::from_file(in_file_muon_lowpt_id_);
  cs_muon_mceff_       = correction::CorrectionSet::from_file(in_file_muon_mceff_);
  cs_pileup_           = correction::CorrectionSet::from_file(in_file_pu_);
  cs_btag_             = correction::CorrectionSet::from_file(in_file_btag_);
  map_electron_        = cs_electron_->at("ElectronWPL");
  map_photon_id_       = cs_photon_->at("UL-Photon-ID-SF");
  map_photon_csev_     = cs_photon_->at("UL-Photon-CSEV-SF");
  map_muon_looseid_    = cs_muon_->at("NUM_LooseID_DEN_genTracks");
  map_muon_highptid_   = cs_muon_->at("NUM_HighPtID_DEN_genTracks");
  map_muon_iso_        = cs_muon_->at("NUM_LooseRelIso_DEN_LooseID");
  map_muon_lowpt_reco_ = cs_muon_lowpt_reco_->at("NUM_TrackerMuons_DEN_genTracks");
  map_muon_lowpt_id_   = cs_muon_lowpt_id_->at("NUM_LooseID_DEN_TrackerMuons");
  map_muon_mceff_      = cs_muon_mceff_->at("Muon_LooseID_MCeff");
  map_btag_            = cs_btag_->at("deepCSV_mujets");
  // DeepJet can be used instead of DeepCSV
  // map_btag_           = cs_btag_->at("deepJet_mujets");
  map_pileup_          = cs_pileup_->at(puName_);
}

// Electron MVA ID Scale Factors
// note: electron prodcer and gen particle producer should already have been run
void EventWeighter::ElectronIDSF(pico_tree &pico, float &w_el_id, std::vector<float> &sys_lep){
  float sf_tot = 1.0;
  float sf_tot_up = 1.0;
  float sf_tot_dn = 1.0;
  //SF logic: for each true electron, weight by prob(MC)/prob(data) which is
  //prob(reco and pass id) or prob(fail)=1-prob(reco and pass id)
  for (unsigned imc = 0; imc < pico.out_mc_id().size(); imc++) {
    if (abs(pico.out_mc_id().at(imc))==11 && ((pico.out_mc_statusflag().at(imc) & 0x2000)!=0)) {
      //is electron and last copy
      if ((pico.out_mc_pt().at(imc)<10) || (fabs(pico.out_mc_eta().at(imc))>2.5)) continue;
      bool pass_id = false;
      float reco_pt = -999;
      float reco_eta = -999;
      float min_dr = 999;
      for(unsigned iel = 0; iel < pico.out_el_sig().size(); ++iel){
        if (pico.out_el_sig().at(iel)) {
          float dr = dR(pico.out_mc_eta().at(imc),pico.out_el_eta().at(iel),
                        pico.out_mc_phi().at(imc),pico.out_el_phi().at(iel));
          if (dr < 0.2 && dr < min_dr) {
            reco_pt = pico.out_el_pt().at(iel);
            reco_eta = pico.out_el_eta().at(iel);
            min_dr = dr;
            pass_id = true;
          }
        }
      }
      if (reco_pt < 0) {
        reco_pt = pico.out_mc_pt().at(imc);
        reco_eta = pico.out_mc_eta().at(imc);
      }
      float mc_eff = map_electron_->evaluate({"effmc", reco_eta, reco_pt});
      float data_eff = map_electron_->evaluate({"effdata", reco_eta, reco_pt});
      float mc_unc = map_electron_->evaluate({"systmc", reco_eta, reco_pt});
      float data_unc = map_electron_->evaluate({"systdata", reco_eta, reco_pt});
      float data_eff_up = data_eff+data_unc;
      float data_eff_dn = data_eff-data_unc;
      float mc_eff_up = mc_eff+mc_unc;
      float mc_eff_dn = mc_eff-mc_unc;
      if (data_eff_up > 1.0) data_eff_up = 1.0;
      if (data_eff_dn < 0.0) data_eff_dn = 0.0;
      if (mc_eff_up > 1.0) mc_eff_up = 1.0;
      if (mc_eff_dn < 0.0) mc_eff_dn = 0.0;
      //for variations consider "worst case": data overestimated and mc 
      //underestimated or vice-versa
      float sf = data_eff/mc_eff;
      float sf_up = data_eff_up/mc_eff_dn;
      float sf_dn = data_eff_dn/mc_eff_up;
      if (!pass_id) {
        sf = (1.0-data_eff)/(1.0-mc_eff);
        sf_up = (1.0-data_eff_up)/(1.0-mc_eff_dn);
        sf_dn = (1.0-data_eff_dn)/(1.0-mc_eff_up);
      }
      if (isnan(sf)||isinf(sf)) sf = 1.0;
      if (isnan(sf_up)||isinf(sf_up)) sf_up = 1.0;
      if (isnan(sf_dn)||isinf(sf_dn)) sf_dn = 1.0;
      sf_tot *= sf;
      sf_tot_up *= sf_up;
      sf_tot_dn *= sf_dn;
    }
  }
  w_el_id = sf_tot;
  sys_lep[0] *= sf_tot_up;
  sys_lep[1] *= sf_tot_dn;
}

// Photon MVA ID Scale Factors
void EventWeighter::PhotonIDSF(pico_tree &pico, float &w_photon_id){
  double sf_tot = 1.0;
  for(size_t i = 0; i < pico.out_photon_pt().size(); ++i){
    if(pico.out_photon_sig().at(i) && pico.out_photon_pt().at(i) > 20.0){
      auto sf = map_photon_id_->evaluate({key_, "sf", "wp90", std::abs(pico.out_photon_eta().at(i)), pico.out_photon_pt().at(i)});
      sf_tot *= sf;
    }
  }
  w_photon_id = sf_tot;
}

// Photon CSEV Scale Factors
void EventWeighter::PhotonCSEVSF(pico_tree &pico, float &w_photon_csev){
  double sf_tot = 1.0;
  string category = "";
  for(size_t i = 0; i < pico.out_photon_pt().size(); ++i){
    if(pico.out_photon_sig().at(i)){
      if (std::abs(pico.out_photon_eta().at(i)) < 1.5){
        category = "EBInc";
        if (std::abs(pico.out_photon_r9().at(i)) > 0.94) {
          category = "EBHighR9";
        } else {
          category = "EBLowR9";
        }
      } else {
        category = "EEInc";
        if (std::abs(pico.out_photon_r9().at(i)) > 0.94) {
          category = "EEHighR9";
        } else {
          category = "EELowR9";
        }
      }
      auto sf = map_photon_csev_->evaluate({key_, "sf", "MVA", category});
      sf_tot *= sf;
    }
  }
  w_photon_csev = sf_tot;
}

// Total Muon Scale Factors
void EventWeighter::MuonTotalSF(pico_tree &pico, float &w_muon_tot, std::vector<float> &sys_lep){
  float sf_tot = 1.0;
  float sf_tot_up = 1.0;
  float sf_tot_dn = 1.0;
  //SF logic: for each true electron, weight by prob(MC)/prob(data) which is
  //prob(reco and pass id) or prob(fail)=1-prob(reco and pass id)
  for (unsigned imc = 0; imc < pico.out_mc_id().size(); imc++) {
    if (abs(pico.out_mc_id().at(imc))==13 && ((pico.out_mc_statusflag().at(imc) & 0x2000)!=0)) {
      //is muon and last copy
      if (pico.out_mc_pt().at(imc)<5 || fabs(pico.out_mc_eta().at(imc))>2.4) continue;
      bool pass_id = false;
      float reco_pt = -999;
      float reco_eta = -999;
      float min_dr = 999;
      for(unsigned imu = 0; imu < pico.out_mu_sig().size(); ++imu){
        if (pico.out_mu_sig().at(imu)) {
          float dr = dR(pico.out_mc_eta().at(imc),pico.out_mu_eta().at(imu),
                        pico.out_mc_phi().at(imc),pico.out_mu_phi().at(imu));
          if (dr < 0.1 && dr < min_dr) {
            reco_pt = pico.out_mu_pt().at(imu);
            reco_eta = pico.out_mu_eta().at(imu);
            min_dr = dr;
            pass_id = true;
          }
        }
      }
      if (reco_pt < 0) {
        reco_pt = pico.out_mc_pt().at(imc);
        reco_eta = pico.out_mc_eta().at(imc);
      }
      float sf = 1.0;
      float sf_up = 1.0;
      float sf_dn = 1.0;
      if (reco_pt < 15) {
        float reco_sf = map_muon_lowpt_reco_->evaluate({"sf", std::abs(reco_eta), reco_pt});
        float reco_sf_unc = map_muon_lowpt_reco_->evaluate({"unc", std::abs(reco_eta), reco_pt});
        float id_sf = map_muon_lowpt_id_->evaluate({"sf", std::abs(reco_eta), reco_pt});
        float id_sf_unc = map_muon_lowpt_id_->evaluate({"unc", std::abs(reco_eta), reco_pt});
        sf = reco_sf*id_sf;
        sf_up = (reco_sf+reco_sf_unc)*(id_sf+id_sf_unc);
        sf_dn = (reco_sf-reco_sf_unc)*(id_sf-id_sf_unc);
        if (sf_dn < 0) sf_dn = 0;
      }
      else if (reco_pt < 200) {
        float id_sf = map_muon_looseid_->evaluate({key_ + "_UL", std::abs(reco_eta), reco_pt, "sf"});
        float id_sf_up = map_muon_looseid_->evaluate({key_ + "_UL", std::abs(reco_eta), reco_pt, "systup"});
        float id_sf_dn = map_muon_looseid_->evaluate({key_ + "_UL", std::abs(reco_eta), reco_pt, "systdown"});
        float iso_sf = map_muon_iso_->evaluate({key_ + "_UL", std::abs(reco_eta), reco_pt, "sf"});
        float iso_sf_up = map_muon_iso_->evaluate({key_ + "_UL", std::abs(reco_eta), reco_pt, "systup"}); 
        float iso_sf_dn = map_muon_iso_->evaluate({key_ + "_UL", std::abs(reco_eta), reco_pt, "systdown"}); 
        sf = id_sf*iso_sf;
        sf_up = id_sf_up*iso_sf_up;
        sf_dn = id_sf_dn*iso_sf_dn;
      }
      else {
        float id_sf = map_muon_highptid_->evaluate({key_ + "_UL", std::abs(reco_eta), reco_pt, "sf"});
        float id_sf_up = map_muon_highptid_->evaluate({key_ + "_UL", std::abs(reco_eta), reco_pt, "systup"});
        float id_sf_dn = map_muon_highptid_->evaluate({key_ + "_UL", std::abs(reco_eta), reco_pt, "systdown"});
        float iso_sf = map_muon_iso_->evaluate({key_ + "_UL", std::abs(reco_eta), reco_pt, "sf"});
        float iso_sf_up = map_muon_iso_->evaluate({key_ + "_UL", std::abs(reco_eta), reco_pt, "systup"}); 
        float iso_sf_dn = map_muon_iso_->evaluate({key_ + "_UL", std::abs(reco_eta), reco_pt, "systdown"}); 
        sf = id_sf*iso_sf;
        sf_up = id_sf_up*iso_sf_up;
        sf_dn = id_sf_dn*iso_sf_dn;
      }
      float mc_eff = map_muon_mceff_->evaluate({"effmc", std::abs(reco_eta), reco_pt});
      float mc_syst = map_muon_mceff_->evaluate({"systmc", std::abs(reco_eta), reco_pt});
      float mc_eff_up = mc_eff+mc_syst;
      float mc_eff_dn = mc_eff-mc_syst;
      if (mc_eff_up > 1.0) mc_eff_up = 1.0;
      if (mc_eff_dn < 0.0) mc_eff_dn = 0.0;
      float data_eff = sf*mc_eff;
      float data_eff_up = sf_up*mc_eff_up;
      float data_eff_dn = sf_dn*mc_eff_dn;
      //for variations consider "worst case": data overestimated and mc 
      //underestimated or vice-versa
      sf = data_eff/mc_eff;
      sf_up = data_eff_up/mc_eff_dn;
      sf_dn = data_eff_dn/mc_eff_up;
      if (!pass_id) {
        sf = (1.0-data_eff)/(1.0-mc_eff);
        sf_up = (1.0-data_eff_up)/(1.0-mc_eff_dn);
        sf_dn = (1.0-data_eff_dn)/(1.0-mc_eff_up);
      }
      if (isnan(sf)||isinf(sf)) sf = 1.0;
      if (isnan(sf_up)||isinf(sf_up)) sf_up = 1.0;
      if (isnan(sf_dn)||isinf(sf_dn)) sf_dn = 1.0;
      sf_tot *= sf;
      sf_tot_up *= sf_up;
      sf_tot_dn *= sf_dn;
    }
  }
  w_muon_tot = sf_tot;
  sys_lep[0] *= sf_tot_up;
  sys_lep[1] *= sf_tot_dn;
}

// Pileup Scale Factors
void EventWeighter::PileupSF(pico_tree &pico, float &w_pu, float &sys_pu_up, float &sys_pu_down){
  auto sf = map_pileup_->evaluate({float(pico.out_npu_tru_mean()), "nominal"});
  w_pu = sf;
  sys_pu_up = (map_pileup_->evaluate({float(pico.out_npu_tru_mean()), "up"}));
  sys_pu_down = (map_pileup_->evaluate({float(pico.out_npu_tru_mean()), "down"}));
}

// b-tagging Scale Factors
void EventWeighter::bTaggingSF(pico_tree &pico, float &w_btag){
  double sf_tot = 1.0;
  auto n_jets = pico.out_jet_isgood().size();
  for(size_t i = 0; i < n_jets; ++i){
    int hadronFlavour = abs(pico.out_jet_hflavor().at(i));
    if(pico.out_jet_isgood().at(i) && hadronFlavour==5 && std::abs(pico.out_jet_eta().at(i)) < 2.4){
      sf_tot *= map_btag_->evaluate({"central", "M", hadronFlavour, std::abs(pico.out_jet_eta().at(i)), pico.out_jet_pt().at(i)}); // M stands for Medium WP
    }
  }
  w_btag = sf_tot;
}
