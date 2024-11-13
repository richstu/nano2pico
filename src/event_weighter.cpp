#include "event_weighter.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "TFile.h"
#include "TGraphAsymmErrors.h"

#include "utilities.hpp"

using namespace std;

EventWeighter::EventWeighter(string year, const vector<float> &btag_wpts){
  if (year=="2016APV") {
    in_file_electron_        = "data/zgamma/2016preVFP_UL/hzg_elid_2016APV_scalefactors.json";
    in_file_photon_          = "data/zgamma/2016preVFP_UL/photon.json";
    in_file_photon_mceff_    = "data/zgamma/2016preVFP_UL/photon_csev_mceff.json";
    in_file_muon_            = "data/zgamma/2016postVFP_UL/muid_2016_2016APV.json";
    in_file_pu_              = "data/zgamma/2016preVFP_UL/puWeights.json";
    in_file_btag_            = "data/zgamma/2016preVFP_UL/btagging.json";
    in_file_btag_mceff_      = "data/zgamma/2016preVFP_UL/btag_mceff.json";
    key_                     = "2016preVFP";
    puName_                  = "Collisions16_UltraLegacy_goldenJSON";
  } else if (year=="2016") {
    in_file_electron_        = "data/zgamma/2016postVFP_UL/hzg_elid_2016_scalefactors.json";
    in_file_photon_          = "data/zgamma/2016postVFP_UL/photon.json";
    in_file_photon_mceff_    = "data/zgamma/2016postVFP_UL/photon_csev_mceff.json";
    in_file_muon_            = "data/zgamma/2016postVFP_UL/muid_2016_2016APV.json";
    in_file_pu_              = "data/zgamma/2016postVFP_UL/puWeights.json";
    in_file_btag_            = "data/zgamma/2016postVFP_UL/btagging.json";
    in_file_btag_mceff_      = "data/zgamma/2016postVFP_UL/btag_mceff.json";
    key_                     = "2016postVFP";
    puName_                  = "Collisions16_UltraLegacy_goldenJSON";
  } else if (year=="2017") {
    in_file_electron_        = "data/zgamma/2017_UL/hzg_elid_2017_scalefactors.json";
    in_file_photon_          = "data/zgamma/2017_UL/photon.json";
    in_file_photon_mceff_    = "data/zgamma/2017_UL/photon_csev_mceff.json";
    in_file_muon_            = "data/zgamma/2017_UL/muid_2017.json";
    in_file_pu_              = "data/zgamma/2017_UL/puWeights.json";
    in_file_btag_            = "data/zgamma/2017_UL/btagging.json";
    in_file_btag_mceff_      = "data/zgamma/2017_UL/btag_mceff.json";
    key_                     = "2017";
    puName_                  = "Collisions17_UltraLegacy_goldenJSON";
  } else if (year=="2018") {
    in_file_electron_        = "data/zgamma/2018_UL/hzg_elid_2018_scalefactors.json";
    in_file_photon_          = "data/zgamma/2018_UL/photon.json";
    in_file_photon_mceff_    = "data/zgamma/2018_UL/photon_csev_mceff.json";
    in_file_muon_            = "data/zgamma/2018_UL/muid_2018.json";
    in_file_pu_              = "data/zgamma/2018_UL/puWeights.json";
    in_file_btag_            = "data/zgamma/2018_UL/btagging.json";
    in_file_btag_mceff_      = "data/zgamma/2018_UL/btag_mceff.json";
    key_                     = "2018";
    puName_                  = "Collisions18_UltraLegacy_goldenJSON";
  } else if (year=="2022"){
    std::cout<<"Using 2018 JSONs by default for now in event_weighter.cpp"<<endl;
    in_file_electron_        = "data/zgamma/2018_UL/hzg_elid_2018_scalefactors.json";
    in_file_photon_          = "data/zgamma/2018_UL/photon.json";
    in_file_photon_mceff_    = "data/zgamma/2018_UL/photon_csev_mceff.json";
    in_file_muon_            = "data/zgamma/2022/muid_2022.json";
    in_file_pu_              = "data/zgamma/2022/puWeights.json";
    in_file_btag_            = "data/zgamma/2018_UL/btagging.json";
    in_file_btag_mceff_      = "data/zgamma/2018_UL/btag_mceff.json";
    key_                     = "2018";
    puName_                  = "Collisions2022_355100_357900_eraBCD_GoldenJson";
  } else if (year=="2022EE"){
    std::cout<<"Using 2018 JSONs by default for now in event_weighter.cpp"<<endl;
    in_file_electron_        = "data/zgamma/2018_UL/hzg_elid_2018_scalefactors.json";
    in_file_photon_          = "data/zgamma/2018_UL/photon.json";
    in_file_photon_mceff_    = "data/zgamma/2018_UL/photon_csev_mceff.json";
    in_file_muon_            = "data/zgamma/2022EE/muid_2022EE.json";
    in_file_pu_              = "data/zgamma/2022EE/puWeights.json";
    in_file_btag_            = "data/zgamma/2018_UL/btagging.json";
    in_file_btag_mceff_      = "data/zgamma/2018_UL/btag_mceff.json";
    key_                     = "2018";
    puName_                  = "Collisions2022_359022_362760_eraEFG_GoldenJson";
  } else if (year=="2023"){
    std::cout<<"Using 2018 JSONs by default for now in event_weighter.cpp"<<endl;
    in_file_electron_        = "data/zgamma/2018_UL/hzg_elid_2018_scalefactors.json";
    in_file_photon_          = "data/zgamma/2018_UL/photon.json";
    in_file_photon_mceff_    = "data/zgamma/2018_UL/photon_csev_mceff.json";
    in_file_muon_            = "data/zgamma/2023/hzg_muid_2023_scalefactors.json";
    in_file_pu_              = "data/zgamma/2023/puWeights.json";
    in_file_btag_            = "data/zgamma/2018_UL/btagging.json";
    in_file_btag_mceff_      = "data/zgamma/2018_UL/btag_mceff.json";
    key_                     = "2018";
    puName_                  = "Collisions2023_366403_369802_eraBC_GoldenJson";
  } else if (year=="2023BPix"){
    std::cout<<"Using 2018 JSONs by default for now in event_weighter.cpp"<<endl;
    in_file_electron_        = "data/zgamma/2018_UL/hzg_elid_2018_scalefactors.json";
    in_file_photon_          = "data/zgamma/2018_UL/photon.json";
    in_file_photon_mceff_    = "data/zgamma/2018_UL/photon_csev_mceff.json";
    in_file_muon_            = "data/zgamma/2023BPix/hzg_muid_2023BPix_scalefactors.json";
    in_file_pu_              = "data/zgamma/2023BPix/puWeights.json";
    in_file_btag_            = "data/zgamma/2018_UL/btagging.json";
    in_file_btag_mceff_      = "data/zgamma/2018_UL/btag_mceff.json";
    key_                     = "2018";
    puName_                  = "Collisions2023_369803_370790_eraD_GoldenJson";
  } else {
    std::cout<<"Year has not been implemented in event_weighter"<<endl;
  }
  cs_electron_              = correction::CorrectionSet::from_file(in_file_electron_);
  cs_photon_                = correction::CorrectionSet::from_file(in_file_photon_);
  cs_photon_mceff_          = correction::CorrectionSet::from_file(in_file_photon_mceff_);
  cs_muon_                  = correction::CorrectionSet::from_file(in_file_muon_);
  cs_pileup_                = correction::CorrectionSet::from_file(in_file_pu_);
  cs_btag_                  = correction::CorrectionSet::from_file(in_file_btag_);
  cs_btag_mceff_            = correction::CorrectionSet::from_file(in_file_btag_mceff_);
  map_photon_id_            = cs_photon_->at("UL-Photon-ID-SF");
  map_photon_csev_          = cs_photon_->at("UL-Photon-CSEV-SF");
  map_photon_csev_mceff_    = cs_photon_mceff_->at("Photon_CSEV_MCeff");
  map_electron_id_pass_     = cs_electron_->at("sf_pass");
  map_electron_id_pass_unc_ = cs_electron_->at("unc_pass");
  map_electron_id_fail_     = cs_electron_->at("sf_fail");
  map_electron_id_fail_unc_ = cs_electron_->at("unc_fail");
  map_muon_id_pass_         = cs_muon_->at("sf_pass");
  map_muon_id_pass_unc_     = cs_muon_->at("unc_pass");
  map_muon_id_fail_         = cs_muon_->at("sf_fail");
  map_muon_id_fail_unc_     = cs_muon_->at("unc_fail");
  map_btag_                 = cs_btag_->at("deepJet_mujets");
  map_udsgtag_              = cs_btag_->at("deepJet_incl");
  // DeepCSV can be used instead of DeepJet
  // map_btag_              = cs_btag_->at("deepCSV_mujets");
  map_pileup_               = cs_pileup_->at(puName_);
  year_                     = year;
  btag_wp_loose_            = btag_wpts[0];
  btag_wp_medium_           = btag_wpts[1];
  btag_wp_tight_            = btag_wpts[2];
}

// Electron MVA ID Scale Factors
// note: electron prodcer and gen particle producer should already have been run
void EventWeighter::ElectronSF(pico_tree &pico){
  float sf_tot = 1.0;
  float sf_tot_up = 1.0;
  float sf_tot_dn = 1.0;
  //SF logic: for each true electron, weight by prob(MC)/prob(data) which is
  //prob(reco and pass id) or prob(fail)=1-prob(reco and pass id)
  for (unsigned imc = 0; imc < pico.out_mc_id().size(); imc++) {
    if (abs(pico.out_mc_id().at(imc))==11 && 
        ((pico.out_mc_statusflag().at(imc) & 0x2000)!=0) &&
        ((pico.out_mc_statusflag().at(imc) & 0x1) != 0)) {
      //is electron and last copy and prompt
      if ((pico.out_mc_pt().at(imc)<7.0f) || (fabs(pico.out_mc_eta().at(imc))>2.5f)) continue;
      bool pass_id = false;
      float reco_pt = -999;
      float reco_eta = -999;
      float min_dr = 999;
      for(unsigned iel = 0; iel < pico.out_el_sig().size(); ++iel){
        if (pico.out_el_sig().at(iel)) {
          float dr = dR(pico.out_mc_eta().at(imc),pico.out_el_eta().at(iel),
                        pico.out_mc_phi().at(imc),pico.out_el_phi().at(iel));
          if (dr < 0.4f && dr < min_dr) {
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
      float sf = 1.0f;
      float sf_up = 1.0f;
      float sf_dn = 1.0f;
      if (pass_id) {
        sf = map_electron_id_pass_->evaluate({reco_pt,reco_eta});
        float unc = map_electron_id_pass_unc_->evaluate({reco_pt,reco_eta});
        sf_up = sf+unc;
        sf_dn = sf-unc;
      }
      else {
        sf = map_electron_id_fail_->evaluate({reco_pt,reco_eta});
        float unc = map_electron_id_fail_unc_->evaluate({reco_pt,reco_eta});
        sf_up = sf+unc;
        sf_dn = sf-unc;
      }
      sf_tot *= sf;
      sf_tot_up *= sf_up;
      sf_tot_dn *= fmax(sf_dn,0.0);
    }
  }
  pico.out_w_el() = sf_tot;
  pico.out_sys_el().resize(2,1.); 
  pico.out_sys_el()[0] = sf_tot_up;
  pico.out_sys_el()[1] = sf_tot_dn;
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
void EventWeighter::PhotonCSEVSF(pico_tree &pico, float &w_photon_csev, std::vector<float> &sys_photon_csev){
  double sf_tot = 1.0;
  double sf_tot_up = 1.0;
  double sf_tot_dn = 1.0;
  std::string category = "";
  for(size_t iph = 0; iph < pico.out_photon_pt().size(); ++iph){
    //apply to truth photons passing ID
    float pt = pico.out_photon_pt().at(iph);
    float mva = pico.out_photon_idmva().at(iph);
    bool pass_mva = ((pico.out_photon_isScEtaEB().at(iph) && mva>-0.4) ||
                     (pico.out_photon_isScEtaEE().at(iph) && mva>-0.58));
    float drmin = pico.out_photon_drmin().at(iph);
    bool eveto = pico.out_photon_elveto().at(iph);
    if (pt < 15.0f || !pass_mva || drmin < 0.4f) continue;
    if (pico.out_photon_pflavor().at(iph) != 1) continue;
    //find category
    if (pico.out_photon_isScEtaEB().at(iph)) {
      if (fabs(pico.out_photon_r9().at(iph)) > 0.94f) {
        category = "EBHighR9";
      } else {
        category = "EBLowR9";
      }
    } else {
      if (fabs(pico.out_photon_r9().at(iph)) > 0.94f) {
        category = "EEHighR9";
      } else {
        category = "EELowR9";
      }
    }
    float sf = map_photon_csev_->evaluate({key_, "sf", "MVA", category});
    float sf_up = map_photon_csev_->evaluate({key_, "sfup", "MVA", category});
    float sf_dn = map_photon_csev_->evaluate({key_, "sfdown", "MVA", category});
    float mc_eff = map_photon_csev_mceff_->evaluate({"effmc", category});
    float mc_syst_up = map_photon_csev_mceff_->evaluate({"systmc_up", category});
    float mc_syst_dn = map_photon_csev_mceff_->evaluate({"systmc_dn", category});
    float data_eff = 1.0;
    float data_syst_up = 0.0;
    float data_syst_dn = 0.0;
    float fail_sf = 1.0;
    float fail_unc_up = 0.0;
    float fail_unc_dn = 0.0;
    propagate_uncertainty_product(mc_eff, mc_syst_up, sf, fabs(sf_up-sf), data_eff, data_syst_up);
    propagate_uncertainty_product(mc_eff, mc_syst_dn, sf, fabs(sf-sf_dn), data_eff, data_syst_dn);
    propagate_uncertainty_ratio(1.0-data_eff, data_syst_up, 1.0-mc_eff, mc_syst_dn, fail_sf, fail_unc_up);
    propagate_uncertainty_ratio(1.0-data_eff, data_syst_dn, 1.0-mc_eff, mc_syst_up, fail_sf, fail_unc_dn);
    if (eveto) { //passes csev
      sf_tot *= sf;
      sf_tot_up *= sf_up;
      sf_tot_dn *= sf_dn;
    }
    else { //fails csev
      sf_tot *= fail_sf;
      sf_tot_up *= fmax(fail_sf-fail_unc_up,0);
      sf_tot_dn *= (fail_sf+fail_unc_dn);
    }
  }
  w_photon_csev = sf_tot;
  sys_photon_csev[0] *= sf_tot_up;
  sys_photon_csev[1] *= sf_tot_dn;
}

// Muon Scale Factors
void EventWeighter::MuonSF(pico_tree &pico){
  float sf_tot = 1.0f;
  float sf_tot_up = 1.0f;
  float sf_tot_dn = 1.0f;
  //SF logic: for each true electron, weight by prob(MC)/prob(data) which is
  //prob(reco and pass id) or prob(fail)=1-prob(reco and pass id)
  for (unsigned imc = 0; imc < pico.out_mc_id().size(); imc++) {
    if (abs(pico.out_mc_id().at(imc))==13 && ((pico.out_mc_statusflag().at(imc) & 0x2000)!=0)) {
      //is muon and last copy
      if (pico.out_mc_pt().at(imc)<5.0f || fabs(pico.out_mc_eta().at(imc))>2.4f) continue;
      bool pass_id = false;
      float reco_pt = -999.0f;
      float reco_eta = -999.0f;
      float min_dr = 999.0f;
      for(unsigned imu = 0; imu < pico.out_mu_sig().size(); ++imu){
        if (pico.out_mu_sig().at(imu)) {
          float dr = dR(pico.out_mc_eta().at(imc),pico.out_mu_eta().at(imu),
                        pico.out_mc_phi().at(imc),pico.out_mu_phi().at(imu));
          if (dr < 0.1f && dr < min_dr) {
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
      if (year_=="2023" || year_=="2023BPix") {
        reco_eta = fabs(reco_eta);
      }
      float sf = 1.0f;
      float sf_up = 1.0f;
      float sf_dn = 1.0f;
      if (pass_id) {
        sf = map_muon_id_pass_->evaluate({reco_pt,reco_eta});
        float unc = map_muon_id_pass_unc_->evaluate({reco_pt,reco_eta});
        sf_up = sf+unc;
        sf_dn = sf-unc;
      }
      else {
        sf = map_muon_id_fail_->evaluate({reco_pt,reco_eta});
        float unc = map_muon_id_fail_unc_->evaluate({reco_pt,reco_eta});
        sf_up = sf+unc;
        sf_dn = sf-unc;
      }
      sf_tot *= sf;
      sf_tot_up *= sf_up;
      sf_tot_dn *= fmax(sf_dn,0.0);
    }
  }
  pico.out_w_mu() = sf_tot;
  pico.out_sys_mu().resize(2,1.); 
  pico.out_sys_mu()[0] = sf_tot_up;
  pico.out_sys_mu()[1] = sf_tot_dn;
}

// Pileup Scale Factors
void EventWeighter::PileupSF(pico_tree &pico){
  pico.out_w_pu() = map_pileup_->evaluate({float(pico.out_npu_tru_mean()), "nominal"});
  pico.out_sys_pu().resize(2, 1.);
  pico.out_sys_pu()[0] = (map_pileup_->evaluate({float(pico.out_npu_tru_mean()), "up"}));
  pico.out_sys_pu()[1] = (map_pileup_->evaluate({float(pico.out_npu_tru_mean()), "down"}));
}

// b-tagging Scale Factors
void EventWeighter::bTaggingSF(pico_tree &pico){
  float sf_tot_nm = 1.0;
  float sf_tot_up_bc = 1.0;
  float sf_tot_dn_bc = 1.0;
  float sf_tot_up_udsg = 1.0;
  float sf_tot_dn_udsg = 1.0;
  for (unsigned ijet = 0; ijet < pico.out_jet_pt().size(); ijet++) {

    if(pico.out_jet_isgood().at(ijet) && std::abs(pico.out_jet_eta().at(ijet)) < 2.4){

      //get true flavor
      int jet_flavor = abs(pico.out_jet_hflavor().at(ijet));
      if (jet_flavor != 5 && jet_flavor != 4) jet_flavor = 0;
      correction::Correction::Ref *btag_map = &map_btag_;
      std::string mc_string = "Btag_b_WP";
      if (jet_flavor == 4) mc_string = "Btag_c_WP";
      if (jet_flavor == 0) {
        mc_string = "Btag_uds_WP";
        btag_map = &map_udsgtag_;
      }

      //calculate probability to be in exclusive categories (i.e. tight, 
      // medium-but-not-tight, loose-but-not-medium, not-loose) in data and MC.
      //This can be done by combining SFs (data/MC) with MC efficiencies
      //Then use prob(data)/prob(MC) to get final SFs
      float cat_data_eff(1.0), cat_data_eff_up(1.0), cat_data_eff_dn(1.0);
      float cat_mc_eff(1.0), cat_mc_eff_up(1.0), cat_mc_eff_dn(1.0);
      if (pico.out_jet_deepflav().at(ijet) > btag_wp_tight_) { 
        float t_sf = (*btag_map)->evaluate({"central", "T", jet_flavor, std::abs(pico.out_jet_eta().at(ijet)), pico.out_jet_pt().at(ijet)});
        float t_sf_up = (*btag_map)->evaluate({"up_correlated", "T", jet_flavor, std::abs(pico.out_jet_eta().at(ijet)), pico.out_jet_pt().at(ijet)});
        float t_sf_dn = (*btag_map)->evaluate({"down_correlated", "T", jet_flavor, std::abs(pico.out_jet_eta().at(ijet)), pico.out_jet_pt().at(ijet)});
        float t_mc_eff = cs_btag_mceff_->at((mc_string+"tight_MCeff").c_str())->evaluate({"effmc", pico.out_jet_eta().at(ijet), pico.out_jet_pt().at(ijet)});
        float t_mc_syst = cs_btag_mceff_->at((mc_string+"tight_MCeff").c_str())->evaluate({"systmc", pico.out_jet_eta().at(ijet), pico.out_jet_pt().at(ijet)});
        cat_mc_eff = t_mc_eff;
        cat_mc_eff_up = t_mc_eff+t_mc_syst;
        cat_mc_eff_dn = t_mc_eff+t_mc_syst;
        cat_data_eff = t_mc_eff*t_sf;
        cat_data_eff_up = t_mc_eff*t_sf_up;
        cat_data_eff_dn = t_mc_eff*t_sf_dn;
      }
      else if (pico.out_jet_deepflav().at(ijet) > btag_wp_medium_) {
        float t_sf = (*btag_map)->evaluate({"central", "T", jet_flavor, std::abs(pico.out_jet_eta().at(ijet)), pico.out_jet_pt().at(ijet)});
        float t_sf_up = (*btag_map)->evaluate({"up_correlated", "T", jet_flavor, std::abs(pico.out_jet_eta().at(ijet)), pico.out_jet_pt().at(ijet)});
        float t_sf_dn = (*btag_map)->evaluate({"down_correlated", "T", jet_flavor, std::abs(pico.out_jet_eta().at(ijet)), pico.out_jet_pt().at(ijet)});
        float t_mc_eff = cs_btag_mceff_->at((mc_string+"tight_MCeff").c_str())->evaluate({"effmc", pico.out_jet_eta().at(ijet), pico.out_jet_pt().at(ijet)});
        float t_mc_syst = cs_btag_mceff_->at((mc_string+"tight_MCeff").c_str())->evaluate({"systmc", pico.out_jet_eta().at(ijet), pico.out_jet_pt().at(ijet)});
        float m_sf = (*btag_map)->evaluate({"central", "M", jet_flavor, std::abs(pico.out_jet_eta().at(ijet)), pico.out_jet_pt().at(ijet)});
        float m_sf_up = (*btag_map)->evaluate({"up_correlated", "M", jet_flavor, std::abs(pico.out_jet_eta().at(ijet)), pico.out_jet_pt().at(ijet)});
        float m_sf_dn = (*btag_map)->evaluate({"down_correlated", "M", jet_flavor, std::abs(pico.out_jet_eta().at(ijet)), pico.out_jet_pt().at(ijet)});
        float m_mc_eff = cs_btag_mceff_->at((mc_string+"medium_MCeff").c_str())->evaluate({"effmc", pico.out_jet_eta().at(ijet), pico.out_jet_pt().at(ijet)});
        float m_mc_syst = cs_btag_mceff_->at((mc_string+"medium_MCeff").c_str())->evaluate({"systmc", pico.out_jet_eta().at(ijet), pico.out_jet_pt().at(ijet)});
        cat_mc_eff = m_mc_eff-t_mc_eff;
        cat_mc_eff_up = m_mc_eff+m_mc_syst-t_mc_eff-t_mc_syst;
        cat_mc_eff_dn = m_mc_eff-m_mc_syst-t_mc_eff+t_mc_syst;
        cat_data_eff = m_mc_eff*m_sf-t_mc_eff*t_sf;
        cat_data_eff_up = m_mc_eff*m_sf_up-t_mc_eff*t_sf_up;
        cat_data_eff_dn = m_mc_eff*m_sf_dn-t_mc_eff*t_sf_dn;
      }
      else if (pico.out_jet_deepflav().at(ijet) > btag_wp_loose_) {
        float m_sf = (*btag_map)->evaluate({"central", "M", jet_flavor, std::abs(pico.out_jet_eta().at(ijet)), pico.out_jet_pt().at(ijet)});
        float m_sf_up = (*btag_map)->evaluate({"up_correlated", "M", jet_flavor, std::abs(pico.out_jet_eta().at(ijet)), pico.out_jet_pt().at(ijet)});
        float m_sf_dn = (*btag_map)->evaluate({"down_correlated", "M", jet_flavor, std::abs(pico.out_jet_eta().at(ijet)), pico.out_jet_pt().at(ijet)});
        float m_mc_eff = cs_btag_mceff_->at((mc_string+"medium_MCeff").c_str())->evaluate({"effmc", pico.out_jet_eta().at(ijet), pico.out_jet_pt().at(ijet)});
        float m_mc_syst = cs_btag_mceff_->at((mc_string+"medium_MCeff").c_str())->evaluate({"systmc", pico.out_jet_eta().at(ijet), pico.out_jet_pt().at(ijet)});
        float l_sf = (*btag_map)->evaluate({"central", "L", jet_flavor, std::abs(pico.out_jet_eta().at(ijet)), pico.out_jet_pt().at(ijet)});
        float l_sf_up = (*btag_map)->evaluate({"up_correlated", "L", jet_flavor, std::abs(pico.out_jet_eta().at(ijet)), pico.out_jet_pt().at(ijet)});
        float l_sf_dn = (*btag_map)->evaluate({"down_correlated", "L", jet_flavor, std::abs(pico.out_jet_eta().at(ijet)), pico.out_jet_pt().at(ijet)});
        float l_mc_eff = cs_btag_mceff_->at((mc_string+"loose_MCeff").c_str())->evaluate({"effmc", pico.out_jet_eta().at(ijet), pico.out_jet_pt().at(ijet)});
        float l_mc_syst = cs_btag_mceff_->at((mc_string+"loose_MCeff").c_str())->evaluate({"systmc", pico.out_jet_eta().at(ijet), pico.out_jet_pt().at(ijet)});
        cat_mc_eff = l_mc_eff-m_mc_eff;
        cat_mc_eff_up = l_mc_eff+l_mc_syst-m_mc_eff-m_mc_syst;
        cat_mc_eff_dn = l_mc_eff-l_mc_syst-m_mc_eff+m_mc_syst;
        cat_data_eff = l_mc_eff*l_sf-m_mc_eff*m_sf;
        cat_data_eff_up = l_mc_eff*l_sf_up-m_mc_eff*m_sf_up;
        cat_data_eff_dn = l_mc_eff*l_sf_dn-m_mc_eff*m_sf_dn;
      }
      else {
        float l_sf = (*btag_map)->evaluate({"central", "L", jet_flavor, std::abs(pico.out_jet_eta().at(ijet)), pico.out_jet_pt().at(ijet)});
        float l_sf_up = (*btag_map)->evaluate({"up_correlated", "L", jet_flavor, std::abs(pico.out_jet_eta().at(ijet)), pico.out_jet_pt().at(ijet)});
        float l_sf_dn = (*btag_map)->evaluate({"down_correlated", "L", jet_flavor, std::abs(pico.out_jet_eta().at(ijet)), pico.out_jet_pt().at(ijet)});
        float l_mc_eff = cs_btag_mceff_->at((mc_string+"loose_MCeff").c_str())->evaluate({"effmc", pico.out_jet_eta().at(ijet), pico.out_jet_pt().at(ijet)});
        float l_mc_syst = cs_btag_mceff_->at((mc_string+"loose_MCeff").c_str())->evaluate({"systmc", pico.out_jet_eta().at(ijet), pico.out_jet_pt().at(ijet)});
        cat_mc_eff = 1.0-l_mc_eff;
        cat_mc_eff_up = 1.0-(l_mc_eff+l_mc_syst);
        cat_mc_eff_dn = 1.0-(l_mc_eff-l_mc_syst);
        cat_data_eff = 1.0-l_mc_eff*l_sf;
        cat_data_eff_up = 1.0-l_mc_eff*l_sf_up;
        cat_data_eff_dn = 1.0-l_mc_eff*l_sf_dn;
      }

      //total SF is product of per-jet SFs
      float sf_nm = cat_data_eff/cat_mc_eff;
      float sf_up = cat_data_eff_up/cat_mc_eff_dn;
      float sf_dn = cat_data_eff_dn/cat_mc_eff_up;
      if (isinf(sf_nm)||isnan(sf_nm)) sf_nm = 1.0;
      if (isinf(sf_up)||isnan(sf_up)) sf_up = 1.0;
      if (isinf(sf_dn)||isnan(sf_dn)) sf_dn = 1.0;
      sf_tot_nm *= sf_nm;
      if (jet_flavor != 0) { //bottom and charm
        sf_tot_up_bc *= sf_up;
        sf_tot_dn_bc *= sf_dn;
        sf_tot_up_udsg *= sf_nm;
        sf_tot_dn_udsg *= sf_nm;
      }
      else { //light flavor
        sf_tot_up_bc *= sf_nm;
        sf_tot_dn_bc *= sf_nm;
        sf_tot_up_udsg *= sf_up;
        sf_tot_dn_udsg *= sf_dn;
      }
    } //jet is good
  } //loop over jets

  pico.out_w_bhig_df() = sf_tot_nm;
  pico.out_sys_bchig().resize(2,1.); 
  pico.out_sys_bchig()[0] = sf_tot_up_bc;
  pico.out_sys_bchig()[1] = sf_tot_dn_bc;
  pico.out_sys_udsghig().resize(2,1.); 
  pico.out_sys_udsghig()[0] = sf_tot_up_udsg;
  pico.out_sys_udsghig()[1] = sf_tot_dn_udsg;
}
