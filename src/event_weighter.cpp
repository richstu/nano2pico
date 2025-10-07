#include "event_weighter.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "TFile.h"
#include "TGraphAsymmErrors.h"

#include "photon_shape_weighter.hpp"
#include "zgbkg_isr_weighter.hpp"
#include "utilities.hpp"

using namespace std;

EventWeighter::EventWeighter(string year, const vector<float> &btag_wpts){
  string photon_idmapname = "Photon-ID-SF";
  string photon_csevmapname = "Photon-CSEV-SF";
  string btag_lightname = "deepJet_incl";
  if (year=="2016APV") {
    in_file_electron_         = "data/zgamma/2016preVFP_UL/hzg_elid_2016APV_scalefactors.json";
    in_file_electron_reco_    = "data/zgamma/2016preVFP_UL/electron_recoSF2016preVFP.json";
    in_file_photon_           = "data/zgamma/2016preVFP_UL/photon.json";
    in_file_photon_low_       = "data/zgamma/2016preVFP_UL/hzg_phidel_2016APV_scalefactors.json";
    in_file_photon_mceff_     = "data/zgamma/2016preVFP_UL/photon_wp80mceff_2016APV.json";
    in_file_muon_             = "data/zgamma/2016postVFP_UL/muid_2016_2016APV.json";
    in_file_pu_               = "data/zgamma/2016preVFP_UL/puWeights.json";
    in_file_btag_             = "data/zgamma/2016preVFP_UL/btagging.json";
    in_file_btag_mceff_       = "data/zgamma/2016preVFP_UL/btag_mceff.json";
    in_file_electron_iso0p10_ = "data/zgamma/2016preVFP_UL/hzg_eliso0p1_2016APV_efficiencies.json";
    in_file_electron_iso0p15_ = "data/zgamma/2016preVFP_UL/hzg_eliso0p15_2016APV_efficiencies.json";
    in_file_muon_iso0p10_     = "data/zgamma/2016preVFP_UL/hzg_muiso0p1_2016APV_efficiencies.json";
    in_file_muon_iso0p15_     = "data/zgamma/2016preVFP_UL/hzg_muiso0p15_2016APV_efficiencies.json";
    in_file_ggf_nnlo_         = "data/zgamma/GluGluHToZG_NNLO_reweight_run2.json";
    key_                      = "2016preVFP";
    puName_                   = "Collisions16_UltraLegacy_goldenJSON";
    photon_idmapname          = "UL-Photon-ID-SF";
    photon_csevmapname        = "UL-Photon-CSEV-SF";
    ph_shape_weighter_        = make_unique<rw_mmp>();
    zgbkg_isr_weighter_       = make_unique<kinr2_weighter>();
  } else if (year=="2016") {
    in_file_electron_         = "data/zgamma/2016postVFP_UL/hzg_elid_2016_scalefactors.json";
    in_file_electron_reco_    = "data/zgamma/2016postVFP_UL/electron_recoSF2016postVFP.json";
    in_file_photon_           = "data/zgamma/2016postVFP_UL/photon.json";
    in_file_photon_low_       = "data/zgamma/2016postVFP_UL/hzg_phidvalidate_2016_scalefactors.json";
    in_file_photon_mceff_     = "data/zgamma/2016postVFP_UL/photon_wp80mceff_2016.json";
    in_file_muon_             = "data/zgamma/2016postVFP_UL/muid_2016_2016APV.json";
    in_file_pu_               = "data/zgamma/2016postVFP_UL/puWeights.json";
    in_file_btag_             = "data/zgamma/2016postVFP_UL/btagging.json";
    in_file_btag_mceff_       = "data/zgamma/2016postVFP_UL/btag_mceff.json";
    in_file_electron_iso0p10_ = "data/zgamma/2016postVFP_UL/hzg_eliso0p1_2016_efficiencies.json";
    in_file_electron_iso0p15_ = "data/zgamma/2016postVFP_UL/hzg_eliso0p15_2016_efficiencies.json";
    in_file_muon_iso0p10_     = "data/zgamma/2016postVFP_UL/hzg_muiso0p1_2016_efficiencies.json";
    in_file_muon_iso0p15_     = "data/zgamma/2016postVFP_UL/hzg_muiso0p15_2016_efficiencies.json";
    in_file_ggf_nnlo_         = "data/zgamma/GluGluHToZG_NNLO_reweight_run2.json";
    key_                      = "2016postVFP";
    puName_                   = "Collisions16_UltraLegacy_goldenJSON";
    photon_idmapname          = "UL-Photon-ID-SF";
    photon_csevmapname        = "UL-Photon-CSEV-SF";
    ph_shape_weighter_        = make_unique<rw_mmp>();
    zgbkg_isr_weighter_       = make_unique<kinr2_weighter>();
  } else if (year=="2017") {
    in_file_electron_         = "data/zgamma/2017_UL/hzg_elid_2017_scalefactors.json";
    in_file_electron_reco_    = "data/zgamma/2017_UL/electron_recoSF2017.json";
    in_file_photon_           = "data/zgamma/2017_UL/photon.json";
    in_file_photon_low_       = "data/zgamma/2017_UL/hzg_phidvalidate_2017_scalefactors.json";
    in_file_photon_mceff_     = "data/zgamma/2017_UL/photon_wp80mceff_2017.json";
    in_file_muon_             = "data/zgamma/2017_UL/muid_2017.json";
    in_file_pu_               = "data/zgamma/2017_UL/puWeights.json";
    in_file_btag_             = "data/zgamma/2017_UL/btagging.json";
    in_file_btag_mceff_       = "data/zgamma/2017_UL/btag_mceff.json";
    in_file_electron_iso0p10_ = "data/zgamma/2017_UL/hzg_eliso0p1_2017_efficiencies.json";
    in_file_electron_iso0p15_ = "data/zgamma/2017_UL/hzg_eliso0p15_2017_efficiencies.json";
    in_file_muon_iso0p10_     = "data/zgamma/2017_UL/hzg_muiso0p1_2017_efficiencies.json";
    in_file_muon_iso0p15_     = "data/zgamma/2017_UL/hzg_muiso0p15_2017_efficiencies.json";
    in_file_ggf_nnlo_         = "data/zgamma/GluGluHToZG_NNLO_reweight_run2.json";
    key_                      = "2017";
    puName_                   = "Collisions17_UltraLegacy_goldenJSON";
    photon_idmapname          = "UL-Photon-ID-SF";
    photon_csevmapname        = "UL-Photon-CSEV-SF";
    ph_shape_weighter_        = make_unique<rw_mmp>();
    zgbkg_isr_weighter_       = make_unique<kinr2_weighter>();
  } else if (year=="2018") {
    in_file_electron_         = "data/zgamma/2018_UL/hzg_elid_2018_scalefactors.json";
    in_file_electron_reco_    = "data/zgamma/2018_UL/electron_recoSF2018.json";
    in_file_photon_           = "data/zgamma/2018_UL/photon.json";
    in_file_photon_low_       = "data/zgamma/2018_UL/hzg_phidvalidate_2018_scalefactors.json";
    in_file_photon_mceff_     = "data/zgamma/2018_UL/photon_wp80mceff_2018.json";
    in_file_muon_             = "data/zgamma/2018_UL/muid_2018.json";
    in_file_pu_               = "data/zgamma/2018_UL/puWeights.json";
    in_file_btag_             = "data/zgamma/2018_UL/btagging.json";
    in_file_btag_mceff_       = "data/zgamma/2018_UL/btag_mceff.json";
    in_file_electron_iso0p10_ = "data/zgamma/2018_UL/hzg_eliso0p1_2018_efficiencies.json";
    in_file_electron_iso0p15_ = "data/zgamma/2018_UL/hzg_eliso0p15_2018_efficiencies.json";
    in_file_muon_iso0p10_     = "data/zgamma/2018_UL/hzg_muiso0p1_2018_efficiencies.json";
    in_file_muon_iso0p15_     = "data/zgamma/2018_UL/hzg_muiso0p15_2018_efficiencies.json";
    in_file_ggf_nnlo_         = "data/zgamma/GluGluHToZG_NNLO_reweight_run2.json";
    key_                      = "2018";
    puName_                   = "Collisions18_UltraLegacy_goldenJSON";
    photon_idmapname          = "UL-Photon-ID-SF";
    photon_csevmapname        = "UL-Photon-CSEV-SF";
    ph_shape_weighter_        = make_unique<rw_mmp>();
    zgbkg_isr_weighter_       = make_unique<kinr2_weighter>();
  } else if (year=="2022"){
    in_file_electron_         = "data/zgamma/2022/hzg_elid_2022_scalefactors.json";
    in_file_electron_reco_    = "data/zgamma/2022/electron_recoSF2022.json";
    in_file_photon_           = "data/zgamma/2022/photon.json";
    in_file_photon_low_       = "data/zgamma/2022/hzg_phidvalidate_2022_scalefactors.json";
    in_file_photon_mceff_     = "data/zgamma/2022/photon_wp80mceff_2022.json";
    in_file_muon_             = "data/zgamma/2022/muid_2022.json";
    in_file_pu_               = "data/zgamma/2022/puWeights.json";
    in_file_btag_             = "data/zgamma/2022/btagging.json";
    in_file_btag_mceff_       = "data/zgamma/2022/btag_mceff.json";
    in_file_electron_iso0p10_ = "data/zgamma/2022/hzg_eliso0p1_2022_efficiencies.json";
    in_file_electron_iso0p15_ = "data/zgamma/2022/hzg_eliso0p15_2022_efficiencies.json";
    in_file_muon_iso0p10_     = "data/zgamma/2022/hzg_muiso0p1_2022_efficiencies.json";
    in_file_muon_iso0p15_     = "data/zgamma/2022/hzg_muiso0p15_2022_efficiencies.json";
    in_file_ggf_nnlo_         = "data/zgamma/GluGluHToZG_NNLO_reweight_run3.json";
    key_                      = "2022Re-recoBCD";
    puName_                   = "Collisions2022_355100_357900_eraBCD_GoldenJson";
    btag_lightname            = "deepJet_light";
    ph_shape_weighter_        = make_unique<rw_mmp_r3>();
    zgbkg_isr_weighter_       = make_unique<kinr3_weighter>();
  } else if (year=="2022EE"){
    in_file_electron_         = "data/zgamma/2022EE/hzg_elid_2022EE_scalefactors.json";
    in_file_electron_reco_    = "data/zgamma/2022EE/electron_recoSF2022EE.json";
    in_file_photon_           = "data/zgamma/2022EE/photon.json";
    in_file_photon_low_       = "data/zgamma/2022EE/hzg_phidvalidate_2022EE_scalefactors.json";
    in_file_photon_mceff_     = "data/zgamma/2022EE/photon_wp80mceff_2022EE.json";
    in_file_muon_             = "data/zgamma/2022EE/muid_2022EE.json";
    in_file_pu_               = "data/zgamma/2022EE/puWeights.json";
    in_file_btag_             = "data/zgamma/2022EE/btagging.json";
    in_file_btag_mceff_       = "data/zgamma/2022EE/btag_mceff.json";
    in_file_electron_iso0p10_ = "data/zgamma/2022EE/hzg_eliso0p1_2022EE_efficiencies.json";
    in_file_electron_iso0p15_ = "data/zgamma/2022EE/hzg_eliso0p15_2022EE_efficiencies.json";
    in_file_muon_iso0p10_     = "data/zgamma/2022EE/hzg_muiso0p1_2022EE_efficiencies.json";
    in_file_muon_iso0p15_     = "data/zgamma/2022EE/hzg_muiso0p15_2022EE_efficiencies.json";
    in_file_ggf_nnlo_         = "data/zgamma/GluGluHToZG_NNLO_reweight_run3.json";
    key_                      = "2022Re-recoE+PromptFG";
    puName_                   = "Collisions2022_359022_362760_eraEFG_GoldenJson";
    btag_lightname            = "deepJet_light";
    ph_shape_weighter_        = make_unique<rw_mmp_r3>();
    zgbkg_isr_weighter_       = make_unique<kinr3_weighter>();
  } else if (year=="2023"){
    in_file_electron_         = "data/zgamma/2023/hzg_elid_2023_scalefactors.json";
    in_file_electron_reco_    = "data/zgamma/2023/electron_recoSF2023.json";
    in_file_photon_           = "data/zgamma/2023/photon.json";
    in_file_photon_low_       = "data/zgamma/2022EE/hzg_phidvalidate_2022EE_scalefactors.json";
    in_file_photon_mceff_     = "data/zgamma/2023/photon_wp80mceff_2023.json";
    in_file_muon_             = "data/zgamma/2023/hzg_muid_2023_scalefactors.json";
    in_file_pu_               = "data/zgamma/2023/puWeights.json";
    in_file_btag_             = "data/zgamma/2023/btagging.json";
    in_file_btag_mceff_       = "data/zgamma/2023/btag_mceff.json";
    in_file_electron_iso0p10_ = "data/zgamma/2023/hzg_eliso0p1_2023_efficiencies.json";
    in_file_electron_iso0p15_ = "data/zgamma/2023/hzg_eliso0p15_2023_efficiencies.json";
    in_file_muon_iso0p10_     = "data/zgamma/2023/hzg_muiso0p1_2023_efficiencies.json";
    in_file_muon_iso0p15_     = "data/zgamma/2023/hzg_muiso0p15_2023_efficiencies.json";
    in_file_ggf_nnlo_         = "data/zgamma/GluGluHToZG_NNLO_reweight_run3.json";
    key_                      = "2023PromptC";
    puName_                   = "Collisions2023_366403_369802_eraBC_GoldenJson";
    btag_lightname            = "deepJet_light";
    ph_shape_weighter_        = make_unique<rw_mmp_r3>();
    zgbkg_isr_weighter_       = make_unique<kinr3_weighter>();
  } else if (year=="2023BPix"){
    in_file_electron_         = "data/zgamma/2023BPix/hzg_elid_2023BPix_scalefactors.json";
    in_file_electron_reco_    = "data/zgamma/2023BPix/electron_recoSF2023BPix.json";
    in_file_photon_           = "data/zgamma/2023BPix/photon.json";
    in_file_photon_low_       = "data/zgamma/2022EE/hzg_phidvalidate_2022EE_scalefactors.json";
    in_file_photon_mceff_     = "data/zgamma/2023BPix/photon_wp80mceff_2023BPix.json";
    in_file_muon_             = "data/zgamma/2023BPix/hzg_muid_2023BPix_scalefactors.json";
    in_file_pu_               = "data/zgamma/2023BPix/puWeights.json";
    in_file_btag_             = "data/zgamma/2023BPix/btagging.json";
    in_file_btag_mceff_       = "data/zgamma/2023BPix/btag_mceff.json";
    in_file_electron_iso0p10_ = "data/zgamma/2023BPix/hzg_eliso0p1_2023BPix_efficiencies.json";
    in_file_electron_iso0p15_ = "data/zgamma/2023BPix/hzg_eliso0p15_2023BPix_efficiencies.json";
    in_file_muon_iso0p10_     = "data/zgamma/2023BPix/hzg_muiso0p1_2023BPix_efficiencies.json";
    in_file_muon_iso0p15_     = "data/zgamma/2023BPix/hzg_muiso0p15_2023BPix_efficiencies.json";
    in_file_ggf_nnlo_         = "data/zgamma/GluGluHToZG_NNLO_reweight_run3.json";
    key_                      = "2023PromptD";
    puName_                   = "Collisions2023_369803_370790_eraD_GoldenJson";
    btag_lightname            = "deepJet_light";
    cs_electron_bpixhole_     = correction::CorrectionSet::from_file(
        "data/zgamma/2023BPix/hzg_elid_2023BPixHole_scalefactors.json");
    cs_el_hole_iso0p10_       = correction::CorrectionSet::from_file(
        "data/zgamma/2023BPix/hzg_eliso0p1_2023BPixHole_efficiencies.json");
    cs_el_hole_iso0p15_       = correction::CorrectionSet::from_file(
        "data/zgamma/2023BPix/hzg_eliso0p15_2023BPixHole_efficiencies.json");
    ph_shape_weighter_        = make_unique<rw_mmp_r3>();
    zgbkg_isr_weighter_       = make_unique<kinr3_weighter>();
  } else {
    cout<<"Year has not been implemented in event_weighter"<<endl;
  }
  cs_electron_              = correction::CorrectionSet::from_file(in_file_electron_);
  cs_electron_reco_         = correction::CorrectionSet::from_file(in_file_electron_reco_);
  cs_photon_                = correction::CorrectionSet::from_file(in_file_photon_);
  cs_photon_low_            = correction::CorrectionSet::from_file(in_file_photon_low_);
  cs_photon_mceff_          = correction::CorrectionSet::from_file(in_file_photon_mceff_);
  cs_muon_                  = correction::CorrectionSet::from_file(in_file_muon_);
  cs_pileup_                = correction::CorrectionSet::from_file(in_file_pu_);
  cs_btag_                  = correction::CorrectionSet::from_file(in_file_btag_);
  cs_btag_mceff_            = correction::CorrectionSet::from_file(in_file_btag_mceff_);
  cs_el_iso0p10_            = correction::CorrectionSet::from_file(in_file_electron_iso0p10_);
  cs_el_iso0p15_            = correction::CorrectionSet::from_file(in_file_electron_iso0p15_);
  cs_mu_iso0p10_            = correction::CorrectionSet::from_file(in_file_muon_iso0p10_);
  cs_mu_iso0p15_            = correction::CorrectionSet::from_file(in_file_muon_iso0p15_);
  cs_fakephoton_            = correction::CorrectionSet::from_file("data/zgamma/fakephoton.json");
  cs_ggf_nnlo_              = correction::CorrectionSet::from_file(in_file_ggf_nnlo_);
  map_photon_id_            = cs_photon_->at(photon_idmapname);
  map_photon_csev_          = cs_photon_->at(photon_csevmapname);
  map_photon_mceff_         = cs_photon_mceff_->at("effmc");
  map_photon_mcunc_         = cs_photon_mceff_->at("systmc");
  map_electron_id_pass_     = cs_electron_->at("sf_pass");
  map_electron_id_pass_unc_ = cs_electron_->at("unc_pass");
  map_electron_id_fail_     = cs_electron_->at("sf_fail");
  map_electron_id_fail_unc_ = cs_electron_->at("unc_fail");
  map_electron_reco_pass_     = cs_electron_reco_->at("sf_pass");
  map_electron_reco_pass_unc_ = cs_electron_reco_->at("unc_pass");
  map_electron_reco_fail_     = cs_electron_reco_->at("sf_fail");
  map_electron_reco_fail_unc_ = cs_electron_reco_->at("unc_fail");
  map_photon_id_low_pass_     = cs_photon_low_->at("sf_pass");
  map_photon_id_low_pass_unc_ = cs_photon_low_->at("unc_pass");
  map_ggf_nnlo_             = cs_ggf_nnlo_->at("w_nnlo");
  post_bpix_                = false;
  if (year == "2023BPix") {
    map_electron_hole_id_pass_     = cs_electron_bpixhole_->at("sf_pass");
    map_electron_hole_id_pass_unc_ = cs_electron_bpixhole_->at("unc_pass");
    map_electron_hole_id_fail_     = cs_electron_bpixhole_->at("sf_fail");
    map_electron_hole_id_fail_unc_ = cs_electron_bpixhole_->at("unc_fail");
    post_bpix_                     = true;
  }
  map_muon_id_pass_         = cs_muon_->at("sf_pass");
  map_muon_id_pass_unc_     = cs_muon_->at("unc_pass");
  map_muon_id_fail_         = cs_muon_->at("sf_fail");
  map_muon_id_fail_unc_     = cs_muon_->at("unc_fail");
  map_btag_                 = cs_btag_->at("deepJet_comb");
  map_udsgtag_              = cs_btag_->at(btag_lightname);
  map_pileup_               = cs_pileup_->at(puName_);
  map_fakephoton_           = cs_fakephoton_->at("fakephoton_corrections");
  year_                     = year;
  btag_wp_loose_            = btag_wpts[0];
  btag_wp_medium_           = btag_wpts[1];
  btag_wp_tight_            = btag_wpts[2];
}

// Electron Reco+MVA ID Scale Factors
// note: electron and gen particle producer should already have been run
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
      if ((pico.out_mc_pt().at(imc)<7.0f) 
          || (fabs(pico.out_mc_eta().at(imc))>2.5f)) continue;
      bool pass_id = false;
      float reco_pt = -999;
      float reco_eta = -999;
      float reco_phi = -999;
      float min_dr = 999;
      for(unsigned iel = 0; iel < pico.out_el_sig().size(); ++iel){
        if (pico.out_el_sig().at(iel)) {
          float dr = dR(pico.out_mc_eta().at(imc),pico.out_el_eta().at(iel),
                        pico.out_mc_phi().at(imc),pico.out_el_phi().at(iel));
          if (dr < 0.4f && dr < min_dr) {
            reco_pt = pico.out_el_pt().at(iel);
            reco_eta = pico.out_el_eta().at(iel);
            reco_phi = pico.out_el_phi().at(iel);
            min_dr = dr;
            pass_id = true;
          }
        }
      }
      if (reco_pt < 0) {
        reco_pt = pico.out_mc_pt().at(imc);
        reco_eta = pico.out_mc_eta().at(imc);
        reco_phi = pico.out_mc_phi().at(imc);
      }
      float sf_reco = 1.0;
      float unc_reco = 1.0;
      float sf = 1.0;
      float unc = 1.0;
      float sf_up = 1.0;
      float sf_dn = 1.0;
      bool in_bpix_region = (reco_eta > -1.566 && reco_eta < 0.0
                             && reco_phi > -1.2 && reco_phi < -0.8);
      if (pass_id) {
        if (post_bpix_ && in_bpix_region) {
          sf = map_electron_hole_id_pass_->evaluate({reco_pt,reco_eta});
          unc = map_electron_hole_id_pass_unc_->evaluate({reco_pt,reco_eta});
          sf_reco = map_electron_reco_pass_->evaluate({reco_pt,reco_eta});
          unc_reco = map_electron_reco_pass_unc_->evaluate({reco_pt,reco_eta});
        }
        else {
          sf = map_electron_id_pass_->evaluate({reco_pt,reco_eta});
          unc = map_electron_id_pass_unc_->evaluate({reco_pt,reco_eta});
          sf_reco = map_electron_reco_pass_->evaluate({reco_pt,reco_eta});
          unc_reco = map_electron_reco_pass_unc_->evaluate({reco_pt,reco_eta});
        }
      }
      else {
        if (post_bpix_ && in_bpix_region) {
          sf = map_electron_hole_id_pass_->evaluate({reco_pt,reco_eta});
          unc = -1.0*map_electron_hole_id_pass_unc_->evaluate({reco_pt,
                                                               reco_eta});
          sf_reco = map_electron_reco_fail_->evaluate({reco_pt,reco_eta});
          unc_reco = -1.0*map_electron_reco_fail_unc_->evaluate({reco_pt,
                                                                 reco_eta});
        }
        else {
          sf = map_electron_id_fail_->evaluate({reco_pt,reco_eta});
          unc = -1.0*map_electron_id_fail_unc_->evaluate({reco_pt,reco_eta});
        }
      }
      sf_up = (sf+unc)*(sf_reco+unc_reco);
      sf_dn = (sf-unc)*(sf_reco-unc_reco);
      if (isinf(sf) || isnan(sf)) sf = 1.0;
      if (isinf(sf_reco) || isnan(sf_reco)) sf_reco = 1.0;
      if (isinf(sf_up) || isnan(sf_up)) sf_up = 1.0;
      if (isinf(sf_dn) || isnan(sf_dn)) sf_dn = 1.0;
      sf_tot *= sf*sf_reco;
      sf_tot_up *= sf_up;
      sf_tot_dn *= fmax(sf_dn,0.0);
    }
  }
  pico.out_w_el() = sf_tot;
  pico.out_sys_el().resize(2,1.); 
  pico.out_sys_el()[0] = sf_tot_up;
  pico.out_sys_el()[1] = sf_tot_dn;
}

// note: call after ElectronSF
void EventWeighter::ElectronMinisoSF(pico_tree &pico){
  float sf_tot = 1.0;
  float sf_tot_up = 1.0;
  float sf_tot_dn = 1.0;
  for (unsigned iel = 0; iel < pico.out_el_pt().size(); iel++) {
    if (!pico.out_el_sig().at(iel)) continue;
    float pt = pico.out_el_pt().at(iel);
    float eta = pico.out_el_eta().at(iel);
    float phi = pico.out_el_phi().at(iel);
    float miniso = pico.out_el_miniso().at(iel);
    float sf = 1.0;
    float sf_up = 1.0;
    float sf_dn = 1.0;
    bool in_bpix_region = (eta > -1.566 && eta < 0.0
                           && phi > -1.2 && phi < -0.8);
    float data_eff_0p10 = 1.0;
    float simu_eff_0p10 = 1.0;
    float data_eff_0p15 = 1.0;
    float simu_eff_0p15 = 1.0;
    float data_unc_0p10 = 1.0;
    float simu_unc_0p10 = 1.0;
    float data_unc_0p15 = 1.0;
    float simu_unc_0p15 = 1.0;
    if (post_bpix_ && in_bpix_region) {
      data_eff_0p10 = cs_el_hole_iso0p10_->at("effdata")->evaluate({pt,eta});
      simu_eff_0p10 = cs_el_hole_iso0p10_->at("effmc")->evaluate({pt,eta});
      data_eff_0p15 = cs_el_hole_iso0p15_->at("effdata")->evaluate({pt,eta});
      simu_eff_0p15 = cs_el_hole_iso0p15_->at("effmc")->evaluate({pt,eta});
      data_unc_0p10 = cs_el_hole_iso0p10_->at("systdata")->evaluate({pt,eta});
      simu_unc_0p10 = cs_el_hole_iso0p10_->at("systmc")->evaluate({pt,eta});
      data_unc_0p15 = cs_el_hole_iso0p15_->at("systdata")->evaluate({pt,eta});
      simu_unc_0p15 = cs_el_hole_iso0p15_->at("systmc")->evaluate({pt,eta});
    }
    else {
      data_eff_0p10 = cs_el_iso0p10_->at("effdata")->evaluate({pt,eta});
      simu_eff_0p10 = cs_el_iso0p10_->at("effmc")->evaluate({pt,eta});
      data_eff_0p15 = cs_el_iso0p15_->at("effdata")->evaluate({pt,eta});
      simu_eff_0p15 = cs_el_iso0p15_->at("effmc")->evaluate({pt,eta});
      data_unc_0p10 = cs_el_iso0p10_->at("systdata")->evaluate({pt,eta});
      simu_unc_0p10 = cs_el_iso0p10_->at("systmc")->evaluate({pt,eta});
      data_unc_0p15 = cs_el_iso0p15_->at("systdata")->evaluate({pt,eta});
      simu_unc_0p15 = cs_el_iso0p15_->at("systmc")->evaluate({pt,eta});
    }
    float data_eff = 1.0;
    float data_eff_up = 1.0;
    float data_eff_dn = 1.0;
    float simu_eff = 1.0;
    float simu_eff_up = 1.0;
    float simu_eff_dn = 1.0;
    if (miniso < 0.1) {
      data_eff = data_eff_0p10;
      simu_eff = simu_eff_0p10;
      data_eff_up = (data_eff_0p10+data_unc_0p10);
      data_eff_dn = (data_eff_0p10-data_unc_0p10);
      simu_eff_up = (simu_eff_0p10+simu_unc_0p10);
      simu_eff_dn = (simu_eff_0p10-simu_unc_0p10);
    }
    else if (miniso < 0.15) {
      data_eff = (data_eff_0p15-data_eff_0p10);
      data_eff_up = ((data_eff_0p15+data_unc_0p15)
                    -(data_eff_0p10+data_unc_0p10));
      data_eff_dn = ((data_eff_0p15-data_unc_0p15)
                    -(data_eff_0p10-data_unc_0p10));
      simu_eff = (simu_eff_0p15-simu_eff_0p10);
      simu_eff_up = ((simu_eff_0p15+simu_unc_0p15)
                    -(simu_eff_0p10+simu_unc_0p10));
      simu_eff_dn = ((simu_eff_0p15-simu_unc_0p15)
                    -(simu_eff_0p10-simu_unc_0p10));
    }
    else {
      data_eff = (1.0-data_eff_0p15);
      data_eff_up = (1.0-(data_eff_0p15+data_unc_0p15));
      data_eff_dn = (1.0-(data_eff_0p15-data_unc_0p15));
      simu_eff = (1.0-simu_eff_0p15);
      simu_eff_up = (1.0-(simu_eff_0p15+simu_unc_0p15));
      simu_eff_dn = (1.0-(simu_eff_0p15-simu_unc_0p15));
    }
    data_eff = max(data_eff, 0.0f);
    data_eff_up = max(data_eff_up, 0.0f);
    data_eff_dn = max(data_eff_dn, 0.0f);
    simu_eff = max(simu_eff, 0.0f);
    simu_eff_up = max(simu_eff_up, 0.0f);
    simu_eff_dn = max(simu_eff_dn, 0.0f);
    sf = data_eff/simu_eff;
    sf_up = data_eff_up/simu_eff_dn;
    sf_dn = data_eff_dn/simu_eff_up;
    if (isinf(sf) || isnan(sf)) sf = 1.0;
    if (isinf(sf_up) || isnan(sf_up)) sf_up = 1.0;
    if (isinf(sf_dn) || isnan(sf_dn)) sf_dn = 1.0;
    sf_tot *= sf;
    sf_tot_up *= sf_up;
    sf_tot_dn *= sf_dn;
  }
  pico.out_w_el() *= sf_tot;
  pico.out_sys_el()[0] *= sf_tot_up;
  pico.out_sys_el()[1] *= sf_tot_dn;
}

// Photon Total Scale Factors
void EventWeighter::PhotonSF(pico_tree &pico){
  double sf_tot = 1.0;
  double sf_tot_idup = 1.0;
  double sf_tot_iddn = 1.0;
  double sf_tot_evup = 1.0;
  double sf_tot_evdn = 1.0;
  //loop over reco photons since ~100% reco efficiency, only pt/eta cuts 
  //between NanoAOD and pico
  for (unsigned iph = 0; iph < pico.out_photon_pt().size(); iph++) {
    if (pico.out_photon_drmin().at(iph) < 0.3f) continue;
    if (pico.out_photon_pflavor().at(iph) != 1) continue;
    float pt = pico.out_photon_pt().at(iph);
    float eta = pico.out_photon_eta().at(iph);
    float phi = pico.out_photon_phi().at(iph);
    float r9 = pico.out_photon_r9().at(iph);
    if(r9<0) r9=0;//soft fix
    string category = "";
    if (pico.out_photon_isScEtaEB().at(iph))
      category += "EB";
    else
      category += "EE";
    if (r9 > 0.94f)
      category += "HighR9";
    else
      category += "LowR9";
    float ev_sf(1.0), ev_sfup(1.0), ev_sfdn(1.0);
    if (year_=="2016APV" || year_=="2016" 
        || year_=="2017" || year_=="2018") {
      ev_sf = map_photon_csev_->evaluate({key_, "sf", "MVA", category});
      ev_sfup = map_photon_csev_->evaluate({key_, "sfup", "MVA", category});
      ev_sfdn = map_photon_csev_->evaluate({key_, "sfdown", "MVA", category});
    }
    else if (year_=="2022" || year_=="2022EE" 
             || year_=="2023" || year_=="2023BPix") {
      ev_sf = map_photon_csev_->evaluate({key_, "sf", "MVA80", eta, r9});
      ev_sfup = map_photon_csev_->evaluate({key_, "sfup", "MVA80", eta, r9});
      ev_sfdn = map_photon_csev_->evaluate({key_, "sfdown", "MVA80", eta, r9});
    }
    string wpstring = "wp80";
    if (pt<20.0f)
      wpstring = "wp80Below20";
    float id_sf = 1.0;
    float id_sfup = 1.0;
    float id_sfdn = 1.0;
    if (!(year_ == "2023" || year_ == "2023BPix") && pt < 20.0f) {
      id_sf = (map_photon_id_low_pass_->evaluate({pt, eta}));
      float id_unc = map_photon_id_low_pass_unc_->evaluate({pt, eta});
      id_sfup = id_sf+id_unc;
      id_sfdn = id_sf-id_unc;
    }
    else if (year_=="2023"||year_=="2023BPix") {
      id_sf = map_photon_id_->evaluate({key_, "sf", wpstring, eta, pt, phi});
      id_sfup = map_photon_id_->evaluate({key_, "sfup", wpstring, eta, pt, phi});
      id_sfdn = map_photon_id_->evaluate({key_, "sfdown", wpstring, eta, pt, phi});
    }
    else {
      id_sf = map_photon_id_->evaluate({key_, "sf", wpstring, eta, pt});
      id_sfup = map_photon_id_->evaluate({key_, "sfup", wpstring, eta, pt});
      id_sfdn = map_photon_id_->evaluate({key_, "sfdown", wpstring, eta, pt});
    }
    float pass_sf = id_sf*ev_sf;
    float pass_sf_idup = id_sfup*ev_sf;
    float pass_sf_iddn = id_sfdn*ev_sf;
    float pass_sf_evup = id_sf*ev_sfup;
    float pass_sf_evdn = id_sf*ev_sfdn;
    float mc_eff = map_photon_mceff_->evaluate({pt, eta});
    float fail_sf = 1.0;
    float fail_sf_idup = 1.0;
    float fail_sf_iddn = 1.0;
    float fail_sf_evup = 1.0;
    float fail_sf_evdn = 1.0;
    if (mc_eff < 1.0) {
      fail_sf = (1.0-pass_sf*mc_eff)/(1.0-mc_eff);
      fail_sf_idup = (1.0-pass_sf_idup*mc_eff)/(1.0-mc_eff);
      fail_sf_iddn = (1.0-pass_sf_iddn*mc_eff)/(1.0-mc_eff);
      fail_sf_evup = (1.0-pass_sf_evup*mc_eff)/(1.0-mc_eff);
      fail_sf_evdn = (1.0-pass_sf_evdn*mc_eff)/(1.0-mc_eff);
    }
    if (pico.out_photon_sig().at(iph)) {
      sf_tot *= pass_sf;
      sf_tot_idup *= pass_sf_idup;
      sf_tot_iddn *= pass_sf_iddn;
      sf_tot_evup *= pass_sf_evup;
      sf_tot_evdn *= pass_sf_evdn;
    }
    else {
      sf_tot *= fail_sf;
      sf_tot_idup *= fail_sf_idup;
      sf_tot_iddn *= fail_sf_iddn;
      sf_tot_evup *= fail_sf_evup;
      sf_tot_evdn *= fail_sf_evdn;
    }
  }
  pico.out_w_photon() = sf_tot;
  pico.out_sys_photon().resize(2,1.); 
  pico.out_sys_photon()[0] = sf_tot_idup;
  pico.out_sys_photon()[1] = sf_tot_iddn;
  pico.out_sys_photon_csev().resize(2,1.); 
  pico.out_sys_photon_csev()[0] = sf_tot_evup;
  pico.out_sys_photon_csev()[1] = sf_tot_evdn;
}

// Muon Scale Factors
void EventWeighter::MuonSF(pico_tree &pico){
  float sf_tot = 1.0;
  float sf_tot_up = 1.0;
  float sf_tot_dn = 1.0;
  //SF logic: for each true electron, weight by prob(MC)/prob(data) which is
  //prob(reco and pass id) or prob(fail)=1-prob(reco and pass id)
  for (unsigned imc = 0; imc < pico.out_mc_id().size(); imc++) {
    if (abs(pico.out_mc_id().at(imc))==13 
        && ((pico.out_mc_statusflag().at(imc) & 0x1)!=0)
        && ((pico.out_mc_statusflag().at(imc) & 0x2000)!=0)) {
      //is prompt muon and last copy
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
      float sf = 1.0;
      float sf_up = 1.0;
      float sf_dn = 1.0;
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

// note: call after MuonSF
void EventWeighter::MuonMinisoSF(pico_tree &pico){
  float sf_tot = 1.0;
  float sf_tot_up = 1.0;
  float sf_tot_dn = 1.0;
  for (unsigned imu = 0; imu < pico.out_mu_pt().size(); imu++) {
    if (!pico.out_mu_sig().at(imu)) continue;
    float pt = pico.out_mu_pt().at(imu);
    float eta = pico.out_mu_eta().at(imu);
    float miniso = pico.out_mu_miniso().at(imu);
    float sf = 1.0;
    float sf_up = 1.0;
    float sf_dn = 1.0;
    float data_eff_0p10 = cs_mu_iso0p10_->at("effdata")->evaluate({pt,eta});
    float simu_eff_0p10 = cs_mu_iso0p10_->at("effmc")->evaluate({pt,eta});
    float data_eff_0p15 = cs_mu_iso0p15_->at("effdata")->evaluate({pt,eta});
    float simu_eff_0p15 = cs_mu_iso0p15_->at("effmc")->evaluate({pt,eta});
    float data_unc_0p10 = cs_mu_iso0p10_->at("systdata")->evaluate({pt,eta});
    float simu_unc_0p10 = cs_mu_iso0p10_->at("systmc")->evaluate({pt,eta});
    float data_unc_0p15 = cs_mu_iso0p15_->at("systdata")->evaluate({pt,eta});
    float simu_unc_0p15 = cs_mu_iso0p15_->at("systmc")->evaluate({pt,eta});
    float data_eff = 1.0;
    float data_eff_up = 1.0;
    float data_eff_dn = 1.0;
    float simu_eff = 1.0;
    float simu_eff_up = 1.0;
    float simu_eff_dn = 1.0;
    if (miniso < 0.1) {
      data_eff = data_eff_0p10;
      simu_eff = simu_eff_0p10;
      data_eff_up = (data_eff_0p10+data_unc_0p10);
      data_eff_dn = (data_eff_0p10-data_unc_0p10);
      simu_eff_up = (simu_eff_0p10+simu_unc_0p10);
      simu_eff_dn = (simu_eff_0p10-simu_unc_0p10);
    }
    else if (miniso < 0.15) {
      data_eff = (data_eff_0p15-data_eff_0p10);
      data_eff_up = ((data_eff_0p15+data_unc_0p15)
                    -(data_eff_0p10+data_unc_0p10));
      data_eff_dn = ((data_eff_0p15-data_unc_0p15)
                    -(data_eff_0p10-data_unc_0p10));
      simu_eff = (simu_eff_0p15-simu_eff_0p10);
      simu_eff_up = ((simu_eff_0p15+simu_unc_0p15)
                    -(simu_eff_0p10+simu_unc_0p10));
      simu_eff_dn = ((simu_eff_0p15-simu_unc_0p15)
                    -(simu_eff_0p10-simu_unc_0p10));
    }
    else {
      data_eff = (1.0-data_eff_0p15);
      data_eff_up = (1.0-(data_eff_0p15+data_unc_0p15));
      data_eff_dn = (1.0-(data_eff_0p15-data_unc_0p15));
      simu_eff = (1.0-simu_eff_0p15);
      simu_eff_up = (1.0-(simu_eff_0p15+simu_unc_0p15));
      simu_eff_dn = (1.0-(simu_eff_0p15-simu_unc_0p15));
    }
    data_eff = max(data_eff, 0.0f);
    data_eff_up = max(data_eff_up, 0.0f);
    data_eff_dn = max(data_eff_dn, 0.0f);
    simu_eff = max(simu_eff, 0.0f);
    simu_eff_up = max(simu_eff_up, 0.0f);
    simu_eff_dn = max(simu_eff_dn, 0.0f);
    sf = data_eff/simu_eff;
    sf_up = data_eff_up/simu_eff_dn;
    sf_dn = data_eff_dn/simu_eff_up;
    sf_tot *= sf;
    sf_tot_up *= sf_up;
    sf_tot_dn *= sf_dn;
  }
  pico.out_w_mu() *= sf_tot;
  pico.out_sys_mu()[0] *= sf_tot_up;
  pico.out_sys_mu()[1] *= sf_tot_dn;
}

// Pileup Scale Factors
void EventWeighter::PileupSF(pico_tree &pico){
  pico.out_w_pu() = min(map_pileup_->evaluate({float(pico.out_npu_tru()), 
      "nominal"}),10.0);
  pico.out_sys_pu().resize(2, 1.);
  pico.out_sys_pu()[0] = min(map_pileup_->evaluate({float(pico.out_npu_tru()), 
      "up"}),10.0);
  pico.out_sys_pu()[1] = min(map_pileup_->evaluate({float(pico.out_npu_tru()), 
      "down"}),10.0);
}

// b-tagging Scale Factors
void EventWeighter::bTaggingSF(pico_tree &pico){
  float sf_tot_nm = 1.0;
  float sf_tot_up_bc = 1.0;
  float sf_tot_dn_bc = 1.0;
  float sf_tot_up_udsg = 1.0;
  float sf_tot_dn_udsg = 1.0;
  float sf_tot_up_uncorr_bc = 1.0;
  float sf_tot_dn_uncorr_bc = 1.0;
  float sf_tot_up_uncorr_udsg = 1.0;
  float sf_tot_dn_uncorr_udsg = 1.0;
  float sf_tot_wpm = 1.0;
  for (unsigned ijet = 0; ijet < pico.out_jet_pt().size(); ijet++) {

    if(pico.out_jet_isgood().at(ijet) 
       && abs(pico.out_jet_eta().at(ijet)) < 2.4f){

      //get true flavor
      int jet_flavor = abs(pico.out_jet_hflavor().at(ijet));
      if (jet_flavor != 5 && jet_flavor != 4) jet_flavor = 0;
      correction::Correction::Ref *btag_map = &map_btag_;
      string mc_string = "Btag_b_WP";
      if (jet_flavor == 4) mc_string = "Btag_c_WP";
      if (jet_flavor == 0) {
        mc_string = "Btag_uds_WP";
        btag_map = &map_udsgtag_;
      }
      float pt = pico.out_jet_pt().at(ijet);
      float eta = pico.out_jet_eta().at(ijet);
      float abseta = fabs(pico.out_jet_eta().at(ijet));

      //calculate probability to be in exclusive categories (i.e. tight, 
      // medium-but-not-tight, loose-but-not-medium, not-loose) in data and MC.
      //This can be done by combining SFs (data/MC) with MC efficiencies
      //Then use prob(data)/prob(MC) to get final SFs
      float cat_data_eff(1.0), cat_data_eff_up(1.0), cat_data_eff_dn(1.0);
      float cat_data_eff_up_uncorr(1.0), cat_data_eff_dn_uncorr(1.0);
      float cat_mc_eff(1.0);
      //float cat_mc_eff_up(1.0), cat_mc_eff_dn(1.0);
      float cat_mc_eff_wpm(1.0), cat_data_eff_wpm(1.0);
      float t_sf = (*btag_map)->evaluate({"central", "T", 
          jet_flavor, abseta, pt});
      float t_mc_eff = cs_btag_mceff_->at((mc_string+"tight_MCeff").c_str())
          ->evaluate({"effmc", eta, pt});
      //float t_mc_syst = cs_btag_mceff_->at((mc_string+"tight_MCeff").c_str())
      //    ->evaluate({"systmc", eta, pt});
      float m_sf = (*btag_map)->evaluate({"central", "M", 
          jet_flavor, abseta, pt});
      float m_mc_eff = cs_btag_mceff_->at((mc_string+"medium_MCeff").c_str())
          ->evaluate({"effmc", eta, pt});
      //float m_mc_syst = cs_btag_mceff_->at((mc_string+"medium_MCeff").c_str())
      //    ->evaluate({"systmc", eta, pt});
      float l_sf = (*btag_map)->evaluate({"central", "L", 
          jet_flavor, abseta, pt});
      float l_mc_eff = cs_btag_mceff_->at((mc_string+"loose_MCeff").c_str())
          ->evaluate({"effmc", eta, pt});
      //float l_mc_syst = cs_btag_mceff_->at((mc_string+"loose_MCeff").c_str())
      //    ->evaluate({"systmc", eta, pt});
      float t_sf_up(t_sf), t_sf_dn(t_sf);
      float t_sf_up_uncorr(t_sf), t_sf_dn_uncorr(t_sf);
      float m_sf_up(m_sf), m_sf_dn(m_sf);
      float m_sf_up_uncorr(m_sf), m_sf_dn_uncorr(m_sf);
      float l_sf_up(l_sf), l_sf_dn(l_sf);
      float l_sf_up_uncorr(l_sf), l_sf_dn_uncorr(l_sf);
      t_sf_up = (*btag_map)->evaluate({"up_correlated", "T", 
        jet_flavor, abseta, pt});
      t_sf_dn = (*btag_map)->evaluate({"down_correlated", "T",
        jet_flavor, abseta, pt});
      t_sf_up_uncorr = (*btag_map)->evaluate({"up_uncorrelated", "T", 
        jet_flavor, abseta, pt});
      t_sf_dn_uncorr = (*btag_map)->evaluate({"down_uncorrelated", "T", 
        jet_flavor, abseta, pt});
      m_sf_up = (*btag_map)->evaluate({"up_correlated", "M", 
        jet_flavor, abseta, pt});
      m_sf_dn = (*btag_map)->evaluate({"down_correlated", "M", 
        jet_flavor, abseta, pt});
      m_sf_up_uncorr = (*btag_map)->evaluate({"up_uncorrelated", "M", 
        jet_flavor, abseta, pt});
      m_sf_dn_uncorr = (*btag_map)->evaluate({"down_uncorrelated", "M", 
        jet_flavor, abseta, pt});
      l_sf_up = (*btag_map)->evaluate({"up_correlated", "L", 
        jet_flavor, abseta, pt});
      l_sf_dn = (*btag_map)->evaluate({"down_correlated", "L", 
        jet_flavor, abseta, pt});
      l_sf_up_uncorr = (*btag_map)->evaluate({"up_uncorrelated", "L", 
        jet_flavor, abseta, pt});
      l_sf_dn_uncorr = (*btag_map)->evaluate({"down_uncorrelated", "L", 
        jet_flavor, abseta, pt});
      //currently, do not propoagate MC stats (negligible WRT SFs)
      if (pico.out_jet_deepflav().at(ijet) > btag_wp_tight_) { 
        cat_mc_eff = t_mc_eff;
        //cat_mc_eff_up = t_mc_eff+t_mc_syst;
        //cat_mc_eff_dn = t_mc_eff-t_mc_syst;
        cat_data_eff = t_mc_eff*t_sf;
        cat_data_eff_up = t_mc_eff*t_sf_up;
        cat_data_eff_dn = t_mc_eff*t_sf_dn;
        cat_data_eff = t_mc_eff*t_sf;
        cat_data_eff_up_uncorr = t_mc_eff*t_sf_up_uncorr;
        cat_data_eff_dn_uncorr = t_mc_eff*t_sf_dn_uncorr;
      }
      else if (pico.out_jet_deepflav().at(ijet) > btag_wp_medium_) {
        cat_mc_eff = m_mc_eff-t_mc_eff;
        //cat_mc_eff_up = m_mc_eff+m_mc_syst-t_mc_eff-t_mc_syst;
        //cat_mc_eff_dn = m_mc_eff-m_mc_syst-t_mc_eff+t_mc_syst;
        cat_data_eff = m_mc_eff*m_sf-t_mc_eff*t_sf;
        cat_data_eff_up = m_mc_eff*m_sf_up-t_mc_eff*t_sf_up;
        cat_data_eff_dn = m_mc_eff*m_sf_dn-t_mc_eff*t_sf_dn;
        cat_data_eff_up_uncorr = (m_mc_eff*m_sf_up_uncorr-t_mc_eff
                                  *t_sf_up_uncorr);
        cat_data_eff_dn_uncorr = (m_mc_eff*m_sf_dn_uncorr-t_mc_eff
                                  *t_sf_dn_uncorr);
      }
      else if (pico.out_jet_deepflav().at(ijet) > btag_wp_loose_) {
        cat_mc_eff = l_mc_eff-m_mc_eff;
        //cat_mc_eff_up = l_mc_eff+l_mc_syst-m_mc_eff-m_mc_syst;
        //cat_mc_eff_dn = l_mc_eff-l_mc_syst-m_mc_eff+m_mc_syst;
        cat_data_eff = l_mc_eff*l_sf-m_mc_eff*m_sf;
        cat_data_eff_up = l_mc_eff*l_sf_up-m_mc_eff*m_sf_up;
        cat_data_eff_dn = l_mc_eff*l_sf_dn-m_mc_eff*m_sf_dn;
        cat_data_eff_up_uncorr = (l_mc_eff*l_sf_up_uncorr
                                  -m_mc_eff*m_sf_up_uncorr);
        cat_data_eff_dn_uncorr = (l_mc_eff*l_sf_dn_uncorr
                                  -m_mc_eff*m_sf_dn_uncorr);
      }
      else {
        cat_mc_eff = 1.0-l_mc_eff;
        //cat_mc_eff_up = 1.0-(l_mc_eff+l_mc_syst);
        //cat_mc_eff_dn = 1.0-(l_mc_eff-l_mc_syst);
        cat_data_eff = 1.0-l_mc_eff*l_sf;
        cat_data_eff_up = 1.0-l_mc_eff*l_sf_up;
        cat_data_eff_dn = 1.0-l_mc_eff*l_sf_dn;
        cat_data_eff_up_uncorr = 1.0-l_mc_eff*l_sf_up_uncorr;
        cat_data_eff_dn_uncorr = 1.0-l_mc_eff*l_sf_dn_uncorr;
      }

      //currently we overwrite the systematic efficiencies (eff_*) with the
      //the single WP versions. These should be commented out to return to the
      //multi-WP versions
      if (pico.out_jet_deepflav().at(ijet) > btag_wp_medium_) {
        cat_mc_eff_wpm = m_mc_eff;
        //cat_mc_eff_up = m_mc_eff+m_mc_syst;
        //cat_mc_eff_dn = m_mc_eff-m_mc_syst;
        cat_data_eff_wpm = m_mc_eff*m_sf;
        cat_data_eff_up = m_mc_eff*m_sf_up;
        cat_data_eff_dn = m_mc_eff*m_sf_dn;
        cat_data_eff_up_uncorr = m_mc_eff*m_sf_up_uncorr;
        cat_data_eff_dn_uncorr = m_mc_eff*m_sf_dn_uncorr;
      }
      else {
        cat_mc_eff_wpm = 1.0-m_mc_eff;
        //cat_mc_eff_up = 1.0-(m_mc_eff+m_mc_syst);
        //cat_mc_eff_dn = 1.0-(m_mc_eff-m_mc_syst);
        cat_data_eff_wpm = 1.0-(m_mc_eff*m_sf);
        cat_data_eff_up = 1.0-(m_mc_eff*m_sf_up);
        cat_data_eff_dn = 1.0-(m_mc_eff*m_sf_dn);
        cat_data_eff_up_uncorr = 1.0-(m_mc_eff*m_sf_up_uncorr);
        cat_data_eff_dn_uncorr = 1.0-(m_mc_eff*m_sf_dn_uncorr);
      }

      //total SF is product of per-jet SFs
      float sf_nm = cat_data_eff/cat_mc_eff;
      float sf_nm_wpm = cat_data_eff_wpm/cat_mc_eff_wpm;
      float sf_up = cat_data_eff_up/cat_mc_eff_wpm;
      float sf_dn = cat_data_eff_dn/cat_mc_eff_wpm;
      float sf_up_uncorr = cat_data_eff_up_uncorr/cat_mc_eff_wpm;
      float sf_dn_uncorr = cat_data_eff_dn_uncorr/cat_mc_eff_wpm;
      if (isinf(sf_nm)||isnan(sf_nm)) sf_nm = 1.0;
      if (isinf(sf_nm_wpm)||isnan(sf_nm_wpm)) sf_nm_wpm = 1.0;
      if (isinf(sf_up)||isnan(sf_up)) sf_up = 1.0;
      if (isinf(sf_dn)||isnan(sf_dn)) sf_dn = 1.0;
      if (isinf(sf_up_uncorr)||isnan(sf_up_uncorr)) sf_up_uncorr = 1.0;
      if (isinf(sf_dn_uncorr)||isnan(sf_dn_uncorr)) sf_dn_uncorr = 1.0;
      sf_tot_nm *= sf_nm;
      sf_tot_wpm *= sf_nm_wpm;
      if (jet_flavor != 0) { //bottom and charm
        sf_tot_up_bc *= sf_up;
        sf_tot_dn_bc *= sf_dn;
        sf_tot_up_uncorr_bc *= sf_up_uncorr;
        sf_tot_dn_uncorr_bc *= sf_dn_uncorr;
        sf_tot_up_udsg *= sf_nm_wpm;
        sf_tot_dn_udsg *= sf_nm_wpm;
        sf_tot_up_uncorr_udsg *= sf_nm_wpm;
        sf_tot_dn_uncorr_udsg *= sf_nm_wpm;
      }
      else { //light flavor
        sf_tot_up_bc *= sf_nm_wpm;
        sf_tot_dn_bc *= sf_nm_wpm;
        sf_tot_up_uncorr_bc *= sf_nm_wpm;
        sf_tot_dn_uncorr_bc *= sf_nm_wpm;
        sf_tot_up_udsg *= sf_up;
        sf_tot_dn_udsg *= sf_dn;
        sf_tot_up_uncorr_udsg *= sf_up_uncorr;
        sf_tot_dn_uncorr_udsg *= sf_dn_uncorr;
      }
    } //jet is good
  } //loop over jets

  pico.out_w_bhig_df() = sf_tot_nm;
  pico.out_w_btag_df() = sf_tot_wpm;
  pico.out_sys_bchig().resize(2,1.); 
  pico.out_sys_udsghig().resize(2,1.); 
  pico.out_sys_bchig_uncorr().resize(2,1.); 
  pico.out_sys_udsghig_uncorr().resize(2,1.); 
  pico.out_sys_bchig()[0] = sf_tot_up_bc;
  pico.out_sys_bchig()[1] = sf_tot_dn_bc;
  pico.out_sys_udsghig()[0] = sf_tot_up_udsg;
  pico.out_sys_udsghig()[1] = sf_tot_dn_udsg;
  pico.out_sys_bchig_uncorr()[0] = sf_tot_up_uncorr_bc;
  pico.out_sys_bchig_uncorr()[1] = sf_tot_dn_uncorr_bc;
  pico.out_sys_udsghig_uncorr()[0] = sf_tot_up_uncorr_udsg;
  pico.out_sys_udsghig_uncorr()[1] = sf_tot_dn_uncorr_udsg;
}

//If we use mutiple WPs, this method is probably easiest to synchronize across
//groups, so I'll leave it here commented, but opt for the simpler approach
//as we are tentatively only using medium WP
//void EventWeighter::bTaggingSFshape(pico_tree &pico){
//  float sf_tot_nm = 1.0;
//  float sf_tot_up_lf = 1.0;
//  float sf_tot_dn_lf = 1.0;
//  float sf_tot_up_hf = 1.0;
//  float sf_tot_dn_hf = 1.0;
//  float sf_tot_up_lfstats1 = 1.0;
//  float sf_tot_dn_lfstats1 = 1.0;
//  float sf_tot_up_lfstats2 = 1.0;
//  float sf_tot_dn_lfstats2 = 1.0;
//  float sf_tot_up_lfstats1 = 1.0;
//  float sf_tot_dn_lfstats1 = 1.0;
//  float sf_tot_up_lfstats2 = 1.0;
//  float sf_tot_dn_lfstats2 = 1.0;
//  float sf_tot_up_cferr1 = 1.0;
//  float sf_tot_dn_cferr1 = 1.0;
//  float sf_tot_up_cferr2 = 1.0;
//  float sf_tot_dn_cferr2 = 1.0;
//  for (unsigned ijet = 0; ijet < pico.out_jet_pt().size(); ijet++) {
//
//    if(pico.out_jet_isgood().at(ijet) && abs(pico.out_jet_eta().at(ijet)) < 2.4){
//
//      int jet_flavor = abs(pico.out_jet_hflavor().at(ijet));
//      float abseta = fabs(pico.out_jet_eta().at(ijet));
//      float pt = pico.out_jet_pt().at(ijet);
//      float disc = pico.out_jet_deepflav().at(ijet);
//      if (jet_flavor != 5 && jet_flavor != 4) jet_flavor = 0;
//
//      float sf_central = map_btag_->evaluate({"central", jet_flavor, abseta, pt, disc});
//      sf_tot_nm *= sf_central;
//      if (jet_flavor != 4) {
//        sf_tot_up_lf *= map_btag_->evaluate({"up_lf", jet_flavor, abseta, pt, disc});
//        sf_tot_dn_lf *= map_btag_->evaluate({"dn_lf", jet_flavor, abseta, pt, disc});
//        sf_tot_up_hf *= map_btag_->evaluate({"up_hf", jet_flavor, abseta, pt, disc});
//        sf_tot_dn_hf *= map_btag_->evaluate({"dn_hf", jet_flavor, abseta, pt, disc});
//        sf_tot_up_lfstats1 *= map_btag_->evaluate({"up_lfstats1", jet_flavor, abseta, pt, disc});
//        sf_tot_dn_lfstats1 *= map_btag_->evaluate({"dn_lfstats1", jet_flavor, abseta, pt, disc});
//        sf_tot_up_hfstats1 *= map_btag_->evaluate({"up_hfstats1", jet_flavor, abseta, pt, disc});
//        sf_tot_dn_hfstats1 *= map_btag_->evaluate({"dn_hfstats1", jet_flavor, abseta, pt, disc});
//        sf_tot_up_lfstats2 *= map_btag_->evaluate({"up_lfstats2", jet_flavor, abseta, pt, disc});
//        sf_tot_dn_lfstats2 *= map_btag_->evaluate({"dn_lfstats2", jet_flavor, abseta, pt, disc});
//        sf_tot_up_hfstats2 *= map_btag_->evaluate({"up_hfstats2", jet_flavor, abseta, pt, disc});
//        sf_tot_dn_hfstats2 *= map_btag_->evaluate({"dn_hfstats2", jet_flavor, abseta, pt, disc});
//        sf_tot_up_cferr1 *= sf_central;
//        sf_tot_dn_cferr1 *= sf_central;
//        sf_tot_up_cferr2 *= sf_central;
//        sf_tot_dn_cferr2 *= sf_central;
//      }
//      else {
//        sf_tot_up_lf *= sf_central;
//        sf_tot_dn_lf *= sf_central;
//        sf_tot_up_hf *= sf_central;
//        sf_tot_dn_hf *= sf_central;
//        sf_tot_up_lfstats1 *= sf_central;
//        sf_tot_dn_lfstats1 *= sf_central;
//        sf_tot_up_hfstats1 *= sf_central;
//        sf_tot_dn_hfstats1 *= sf_central;
//        sf_tot_up_lfstats2 *= sf_central;
//        sf_tot_dn_lfstats2 *= sf_central;
//        sf_tot_up_hfstats2 *= sf_central;
//        sf_tot_dn_hfstats2 *= sf_central;
//        sf_up_cferr1 *= map_btag_->evaluate({"up_cferr1", jet_flavor, abseta, pt, disc});
//        sf_dn_cferr1 *= map_btag_->evaluate({"dn_cferr1", jet_flavor, abseta, pt, disc});
//        sf_up_cferr2 *= map_btag_->evaluate({"up_cferr2", jet_flavor, abseta, pt, disc});
//        sf_dn_cferr2 *= map_btag_->evaluate({"dn_cferr2", jet_flavor, abseta, pt, disc});
//      }
//
//      sf_tot_nm *= t_sf;
//      sf_tot
//    } //jet is good
//  } //loop over jets
//
//  pico.out_w_bhig_df() = sf_tot_nm;
//  pico.out_sys_bchig().resize(2,1.); 
//  pico.out_sys_bchig()[0] = sf_tot_up_bc;
//  pico.out_sys_bchig()[1] = sf_tot_dn_bc;
//  pico.out_sys_udsghig().resize(2,1.); 
//  pico.out_sys_udsghig()[0] = sf_tot_up_udsg;
//  pico.out_sys_udsghig()[1] = sf_tot_dn_udsg;
//}

// Photon shape SFs, call after photons have been produced
void EventWeighter::PhotonShapeSF(pico_tree &pico){
  //only apply to lead photon, and only if true photon
  if (pico.out_nphoton()==0) {
    pico.out_w_phshape() = 1.0;
    return;
  }
  if (pico.out_photon_pflavor().at(0) != 1) {
    pico.out_w_phshape() = 1.0;
    return;
  }
  vector<float> dnn_input = {pico.out_photon_pt().at(0),
      fabs(pico.out_photon_eta().at(0)), pico.out_photon_idmva().at(0),
      static_cast<float>(pico.out_photon_energyErr().at(0)/
      (pico.out_photon_pt().at(0)*TMath::CosH(pico.out_photon_eta().at(0))))};
  float dnn_output = ph_shape_weighter_->evaluate(dnn_input);
  float weight = (dnn_output/(1.0-dnn_output));
  if (fabs(weight) > 5.0)
    weight = 5.0*(weight/fabs(weight));
  if (isnan(weight) || isinf(weight)) weight = 1.0;
  pico.out_w_phshape() = weight;
}

// Gets DY fake photon weight
void EventWeighter::FakePhotonSF(pico_tree &pico) {
  if (!(pico.out_type() >= 6000 && pico.out_type() < 7000) 
        || pico.out_nphoton() < 1) {
    pico.out_w_fakephoton() = 1.0f;
    return;
  }
  if (pico.out_photon_pflavor()[0]==1) {
    pico.out_w_fakephoton() = 1.0f;
    return;
  }
  string run = "run3";
  int this_photon_isjet = 0;
  if (year_=="2016APV" || year_=="2016" || year_=="2017" || year_=="2018")
    run = "run2";
  if (pico.out_photon_hardprocess()[0])
    this_photon_isjet = 1;
  float ph_pt = pico.out_photon_pt()[0];
  float ph_abseta = fabs(pico.out_photon_eta()[0]);
  float sf_fp = map_fakephoton_->evaluate({run, this_photon_isjet, ph_pt, 
                                           ph_abseta});
  pico.out_w_fakephoton() = sf_fp;
}

// Z ISR/kinematics SF
void EventWeighter::ZISRSF(pico_tree &pico){
  //for DY and DYG samples
  if (!((pico.out_type() >= 6000 && pico.out_type() < 7000) ||
        (pico.out_type() >= 17000 && pico.out_type() < 18000))
        || pico.out_nllphoton() < 1) {
    pico.out_w_isr() = 1.0f;
    return;
  }

  vector<float> dnn_input = {pico.out_photon_idmva()[0],
                             pico.out_photon_mht_dphi()[0], 
                             static_cast<float>(pico.out_njet()), 
                             pico.out_ll_pt()[0], 
                             pico.out_photon_pt()[0], 
                             pico.out_llphoton_pt()[0], 
                             pico.out_mht(), 
                             pico.out_ht()}; 
  float dnn_output = zgbkg_isr_weighter_->evaluate(dnn_input);
  float weight = (dnn_output/(1.0-dnn_output));
  if (fabs(weight) > 5.0)
    weight = 5.0*(weight/fabs(weight));
  if (isnan(weight) || isinf(weight)) weight = 1.0;
  pico.out_w_isr() = weight;
}

// reweight ggF kinematics to NNLO
void EventWeighter::NNLOCorrection(pico_tree &pico){
  float sf = 1.0;
  if(pico.out_type()==200000){
    //generator level Higgs
    int hidx = -1;
    for (unsigned int i(0); i<pico.out_mc_pt().size(); i++){
      if (pico.out_mc_id().at(i)==25){
        hidx = i;
        break;
      }
    }
    if (hidx>=0){
      float h_pt = pico.out_mc_pt().at(hidx);
      sf = map_ggf_nnlo_->evaluate({h_pt});
    }
  }
  pico.out_w_nnlo() = sf;
}
