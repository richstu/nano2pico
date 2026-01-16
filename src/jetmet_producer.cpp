#include "jetmet_producer.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip> 
#include <iostream>
#include <vector>

#include "TLorentzVector.h"

#include "correction.h"
#include "hig_producer.hpp"
#include "met_producer.hpp"
#include "utilities.hpp"

using namespace std;

JetMetProducer::JetMetProducer(int year_, string year_string_, 
                               float nanoaod_version_, 
                               float min_jet_pt_, float max_jet_eta_, 
                               bool isData_, bool is_preUL_, 
                               bool verbose_) : 
    met_producer(MetProducer(year_, isData_, is_preUL_)) {
  year = year_;
  year_string = year_string_;
  isData = isData_;
  is_preUL = is_preUL_;
  verbose = verbose_;
  min_jet_pt = min_jet_pt_;
  max_jet_eta = max_jet_eta_;
  nanoaod_version = nanoaod_version_;
  rng_ = TRandom3(4357);
  if (year_string=="2016APV") {
    cs_jerc_ = correction::CorrectionSet::from_file("data/zgamma/2016preVFP_UL/jet_jerc_2016apv.json");
    map_jes_ = cs_jerc_->at("Summer19UL16APV_V7_MC_Total_AK4PFchs");
    map_jersf_ = cs_jerc_->at("Summer20UL16APV_JRV3_MC_ScaleFactor_AK4PFchs");
    map_jermc_ = cs_jerc_->at("Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs");
    map_jec_.push_back(cs_jerc_->compound().at("Summer19UL16APV_V7_MC_L1L2L3Res_AK4PFchs"));
    map_jec_l1_.push_back(cs_jerc_->at("Summer19UL16APV_V7_MC_L1FastJet_AK4PFchs"));
  }
  else if (year_string=="2016") {
    cs_jerc_ = correction::CorrectionSet::from_file("data/zgamma/2016postVFP_UL/jet_jerc_2016.json");
    map_jes_ = cs_jerc_->at("Summer19UL16_V7_MC_Total_AK4PFchs");
    map_jersf_ = cs_jerc_->at("Summer20UL16_JRV3_MC_ScaleFactor_AK4PFchs");
    map_jermc_ = cs_jerc_->at("Summer20UL16_JRV3_MC_PtResolution_AK4PFchs");
    map_jec_.push_back(cs_jerc_->compound().at("Summer19UL16_V7_MC_L1L2L3Res_AK4PFchs"));
    map_jec_l1_.push_back(cs_jerc_->at("Summer19UL16_V7_MC_L1FastJet_AK4PFchs"));
  }
  else if (year_string=="2017") {
    cs_jerc_ = correction::CorrectionSet::from_file("data/zgamma/2017_UL/jet_jerc_2017.json");
    map_jes_ = cs_jerc_->at("Summer19UL17_V5_MC_Total_AK4PFchs");
    map_jersf_ = cs_jerc_->at("Summer19UL17_JRV2_MC_ScaleFactor_AK4PFchs");
    map_jermc_ = cs_jerc_->at("Summer19UL17_JRV2_MC_PtResolution_AK4PFchs");
    map_jec_.push_back(cs_jerc_->compound().at("Summer19UL17_V5_MC_L1L2L3Res_AK4PFchs"));
    map_jec_l1_.push_back(cs_jerc_->at("Summer19UL17_V5_MC_L1FastJet_AK4PFchs"));
  }
  else if (year_string=="2018") {
    cs_jerc_ = correction::CorrectionSet::from_file("data/zgamma/2018_UL/jet_jerc_2018.json");
    map_jes_ = cs_jerc_->at("Summer19UL18_V5_MC_Total_AK4PFchs");
    map_jersf_ = cs_jerc_->at("Summer19UL18_JRV2_MC_ScaleFactor_AK4PFchs");
    map_jermc_ = cs_jerc_->at("Summer19UL18_JRV2_MC_PtResolution_AK4PFchs");
    map_jec_.push_back(cs_jerc_->compound().at("Summer19UL18_V5_MC_L1L2L3Res_AK4PFchs"));
    map_jec_l1_.push_back(cs_jerc_->at("Summer19UL18_V5_MC_L1FastJet_AK4PFchs"));
  }
  else if (year_string=="2022") {
    cs_jerc_ = correction::CorrectionSet::from_file("data/zgamma/2022/jet_jerc.json");
    if (isData) {
      map_jec_.push_back(cs_jerc_->compound().at("Summer22_22Sep2023_RunCD_V2_DATA_L1L2L3Res_AK4PFPuppi"));
      map_jec_l1_.push_back(cs_jerc_->at("Summer22_22Sep2023_RunCD_V2_DATA_L1FastJet_AK4PFPuppi"));
    }
    else {
      map_jes_ = cs_jerc_->at("Summer22_22Sep2023_V2_MC_Total_AK4PFPuppi");
      map_jersf_ = cs_jerc_->at("Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi");
      map_jermc_ = cs_jerc_->at("Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi");
      map_jec_.push_back(cs_jerc_->compound().at("Summer22_22Sep2023_V2_MC_L1L2L3Res_AK4PFPuppi"));
      map_jec_l1_.push_back(cs_jerc_->at("Summer22_22Sep2023_V2_MC_L1FastJet_AK4PFPuppi"));
    }

    in_file_jetveto_ = "data/zgamma/2022/jetvetomaps.json";
    cs_jetveto_ = correction::CorrectionSet::from_file(in_file_jetveto_);
    map_jetveto_ = cs_jetveto_->at("Summer22_23Sep2023_RunCD_V1");
  }
  else if (year_string=="2022EE") {
    cs_jerc_ = correction::CorrectionSet::from_file("data/zgamma/2022EE/jet_jerc.json");
    if (isData) {
      map_jec_.push_back(cs_jerc_->compound().at("Summer22EE_22Sep2023_RunE_V2_DATA_L1L2L3Res_AK4PFPuppi"));
      map_jec_l1_.push_back(cs_jerc_->at("Summer22EE_22Sep2023_RunE_V2_DATA_L1FastJet_AK4PFPuppi"));
      map_jec_.push_back(cs_jerc_->compound().at("Summer22EE_22Sep2023_RunF_V2_DATA_L1L2L3Res_AK4PFPuppi"));
      map_jec_l1_.push_back(cs_jerc_->at("Summer22EE_22Sep2023_RunF_V2_DATA_L1FastJet_AK4PFPuppi"));
      map_jec_.push_back(cs_jerc_->compound().at("Summer22EE_22Sep2023_RunG_V2_DATA_L1L2L3Res_AK4PFPuppi"));
      map_jec_l1_.push_back(cs_jerc_->at("Summer22EE_22Sep2023_RunG_V2_DATA_L1FastJet_AK4PFPuppi"));
      jec_run_start_.push_back(359022);
      jec_run_end_.push_back(360331);
      jec_run_start_.push_back(360332);
      jec_run_end_.push_back(362180);
      jec_run_start_.push_back(362350);
      jec_run_end_.push_back(362760);
    }
    else {
      map_jes_ = cs_jerc_->at("Summer22EE_22Sep2023_V2_MC_Total_AK4PFPuppi");
      map_jersf_ = cs_jerc_->at("Summer22EE_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi");
      map_jermc_ = cs_jerc_->at("Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi");
      map_jec_.push_back(cs_jerc_->compound().at("Summer22EE_22Sep2023_V2_MC_L1L2L3Res_AK4PFPuppi"));
      map_jec_l1_.push_back(cs_jerc_->at("Summer22EE_22Sep2023_V2_MC_L1FastJet_AK4PFPuppi"));
    }

    in_file_jetveto_ = "data/zgamma/2022EE/jetvetomaps.json";
    cs_jetveto_ = correction::CorrectionSet::from_file(in_file_jetveto_);
    map_jetveto_ = cs_jetveto_->at("Summer22EE_23Sep2023_RunEFG_V1");
  }
  else if (year_string=="2023") {
    cs_jerc_ = correction::CorrectionSet::from_file("data/zgamma/2023/jet_jerc.json");
    if (isData) {
      map_jec_.push_back(cs_jerc_->compound().at("Summer23Prompt23_V2_DATA_L1L2L3Res_AK4PFPuppi"));
      map_jec_l1_.push_back(cs_jerc_->at("Summer23Prompt23_V2_DATA_L1FastJet_AK4PFPuppi"));
    }
    else {
      map_jes_ = cs_jerc_->at("Summer23Prompt23_V2_MC_Total_AK4PFPuppi");
      map_jersf_ = cs_jerc_->at("Summer23Prompt23_RunCv1234_JRV1_MC_ScaleFactor_AK4PFPuppi");
      map_jermc_ = cs_jerc_->at("Summer23Prompt23_RunCv1234_JRV1_MC_PtResolution_AK4PFPuppi");
      map_jec_.push_back(cs_jerc_->compound().at("Summer23Prompt23_V2_MC_L1L2L3Res_AK4PFPuppi"));
      map_jec_l1_.push_back(cs_jerc_->at("Summer23Prompt23_V2_MC_L1FastJet_AK4PFPuppi"));
    }

    in_file_jetveto_ = "data/zgamma/2023/jetvetomaps.json";
    cs_jetveto_ = correction::CorrectionSet::from_file(in_file_jetveto_);
    map_jetveto_ = cs_jetveto_->at("Summer23Prompt23_RunC_V1");
  }
  else if (year_string=="2023BPix") {
    cs_jerc_ = correction::CorrectionSet::from_file("data/zgamma/2023BPix/jet_jerc.json");
    if (isData) {
      map_jec_.push_back(cs_jerc_->compound().at("Summer23BPixPrompt23_V3_DATA_L1L2L3Res_AK4PFPuppi"));
      map_jec_l1_.push_back(cs_jerc_->at("Summer23BPixPrompt23_V3_DATA_L1FastJet_AK4PFPuppi"));
    }
    else {
      map_jes_ = cs_jerc_->at("Summer23BPixPrompt23_V3_MC_Total_AK4PFPuppi");
      map_jersf_ = cs_jerc_->at("Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi");
      map_jermc_ = cs_jerc_->at("Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi");
      map_jec_.push_back(cs_jerc_->compound().at("Summer23BPixPrompt23_V3_MC_L1L2L3Res_AK4PFPuppi"));
      map_jec_l1_.push_back(cs_jerc_->at("Summer23BPixPrompt23_V3_MC_L1FastJet_AK4PFPuppi"));
    }

    in_file_jetveto_ = "data/zgamma/2023BPix/jetvetomaps.json";
    cs_jetveto_ = correction::CorrectionSet::from_file(in_file_jetveto_);
    map_jetveto_ = cs_jetveto_->at("Summer23BPixPrompt23_RunD_V1");
  }
  else if (year_string=="2024") {
    cs_jerc_ = correction::CorrectionSet::from_file("data/zgamma/2024/jet_jerc.json");
    if (isData) {
      map_jec_.push_back(cs_jerc_->compound().at("Summer24Prompt24_V2_DATA_L1L2L3Res_AK4PFPuppi"));
      map_jec_l1_.push_back(cs_jerc_->at("Summer24Prompt24_V2_DATA_L1FastJet_AK4PFPuppi"));
    }
    else {
      map_jes_ = cs_jerc_->at("Summer24Prompt24_V2_MC_Total_AK4PFPuppi");
      //jet_jerc json has these two branches from 2023BPix. . .
      map_jersf_ = cs_jerc_->at("Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi");
      map_jermc_ = cs_jerc_->at("Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi");
      map_jec_.push_back(cs_jerc_->compound().at("Summer24Prompt24_V2_MC_L1L2L3Res_AK4PFPuppi"));
      map_jec_l1_.push_back(cs_jerc_->at("Summer24Prompt24_V2_MC_L1FastJet_AK4PFPuppi"));
    }

    in_file_jetveto_ = "data/zgamma/2024/jetvetomaps.json";
    cs_jetveto_ = correction::CorrectionSet::from_file(in_file_jetveto_);
    map_jetveto_ = cs_jetveto_->at("Summer24Prompt24_RunBCDEFGHI_V1");

    in_file_jetid_ = "data/zgamma/2024/JetID_Run3_Rereco2022CDE_v2.json";
    cs_jetid_ = correction::CorrectionSet::from_file(in_file_jetid_);
    map_jetid_tight_ = cs_jetid_->at("AK4PUPPI_Tight");
    map_jetid_tightlepveto_ = cs_jetid_->at("AK4PUPPI_TightLeptonVeto");
  }
  else {
    cout << "WARNING: No dedicated JEC/JER, defaulting to 2023BPix NanoAODv12(!) treatment." << endl;

    cs_jerc_ = correction::CorrectionSet::from_file("data/zgamma/2023BPix/jet_jerc.json");
    if (isData) {
      map_jec_.push_back(cs_jerc_->compound().at("Summer23BPixPrompt23_V3_DATA_L1L2L3Res_AK4PFPuppi"));
      map_jec_l1_.push_back(cs_jerc_->at("Summer23BPixPrompt23_V3_DATA_L1FastJet_AK4PFPuppi"));
    }
    else {
      map_jes_ = cs_jerc_->at("Summer23BPixPrompt23_V3_MC_Total_AK4PFPuppi");
      map_jersf_ = cs_jerc_->at("Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi");
      map_jermc_ = cs_jerc_->at("Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi");
      map_jec_.push_back(cs_jerc_->compound().at("Summer23BPixPrompt23_V3_MC_L1L2L3Res_AK4PFPuppi"));
      map_jec_l1_.push_back(cs_jerc_->at("Summer23BPixPrompt23_V3_MC_L1FastJet_AK4PFPuppi"));
    }

    in_file_jetveto_ = "data/zgamma/2023BPix/jetvetomaps.json";
    cs_jetveto_ = correction::CorrectionSet::from_file(in_file_jetveto_);
    map_jetveto_ = cs_jetveto_->at("Summer23BPixPrompt23_RunD_V1");

  }
}

JetMetProducer::~JetMetProducer(){
}

float JetMetProducer::GetJEC(float jet_area, float jet_eta, float jet_phi, 
                             float jet_pt, float rho, unsigned int run, 
                             JECType jec_type) {
   if (year_string == "2023BPix") {
     if (jec_type == JECType::L1L2L3) {
       if (isData)
         return map_jec_[0]->evaluate({jet_area, jet_eta, jet_pt, rho, 
                                       jet_phi, static_cast<float>(run)});
       return map_jec_[0]->evaluate({jet_area, jet_eta, jet_pt, rho, jet_phi});
     }
     else {
       return map_jec_l1_[0]->evaluate({jet_area, jet_eta, jet_pt, rho});
     }
   }
   else if (year_string == "2023" && isData) {
     if (jec_type == JECType::L1L2L3) {
       return map_jec_[0]->evaluate({jet_area, jet_eta, jet_pt, rho, 
                                     static_cast<float>(run)});
     }
     else {
       return map_jec_l1_[0]->evaluate({jet_area, jet_eta, jet_pt, rho});
     }
   }
   else if (year_string == "2022EE" && isData) {
     bool found_era = false;
     unsigned int era_idx = 0;
     for (unsigned iera = 0; iera < jec_run_start_.size(); iera++) {
       if (run >= jec_run_start_[iera] && run <= jec_run_end_[iera]) {
         found_era = true;
         era_idx = iera;
         break;
       }
     }
     if (!found_era)
       throw runtime_error("2022EE run number out of bounds (JEC check).");
     if (jec_type == JECType::L1L2L3) {
       return map_jec_[era_idx]->evaluate({jet_area, jet_eta, jet_pt, rho});
     }
     else {
       return map_jec_l1_[era_idx]->evaluate({jet_area, jet_eta, jet_pt, rho});
     }
   }
   else if (year_string == "2024" || year_string == "2025"){
     if (jec_type == JECType::L1L2L3) {
       if (isData)
         return map_jec_[0]->evaluate({jet_area, jet_eta, jet_pt, rho,
                                       jet_phi, static_cast<float>(run)});
       return map_jec_[0]->evaluate({jet_area, jet_eta, jet_pt, rho, jet_phi});
     }
     else {
       return map_jec_l1_[0]->evaluate({jet_area, jet_eta, jet_pt, rho});
     }
   }
   else {
     if (jec_type == JECType::L1L2L3) {
       return map_jec_[0]->evaluate({jet_area, jet_eta, jet_pt, rho});
     }
     else {
       return map_jec_l1_[0]->evaluate({jet_area, jet_eta, jet_pt, rho});
     }
   }
}

// Note this function also writes out MET and its systematic uncertainties 
// from propagating JES/JER uncertainties
void JetMetProducer::PropagateJERC(nano_tree &nano, pico_tree &pico, 
                                   vector<float> &jet_nm_factor, 
                                   vector<float> &jer_up_factor,
                                   vector<float> &jer_dn_factor,
                                   vector<float> &jes_up_factor,
                                   vector<float> &jes_dn_factor) {

  if (year <= 2018 && isData) {
    //JECs already correct- no updating needed
    pico.out_met() = nano.MET_pt();
    pico.out_met_phi() = nano.MET_phi();
    jet_nm_factor.resize(nano.nJet(),1.0);
    WriteMetVariations(nano, pico);
    return;
  }

  //implementation originally based on the following:
  //https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/jetmetUncertainties.py
  float met_x, met_y;
  if (year <= 2018) {
    met_x = nano.MET_pt()*cos(nano.MET_phi());
    met_y = nano.MET_pt()*sin(nano.MET_phi());
  } else {
    // Start from Raw in run 3 to reapply JECs
    met_x = nano.RawPuppiMET_pt()*cos(nano.RawPuppiMET_phi());
    met_y = nano.RawPuppiMET_pt()*sin(nano.RawPuppiMET_phi());
  }
  float met_x_nom = met_x;
  float met_y_nom = met_y;
  float met_x_jesup = 0.0;
  float met_y_jesup = 0.0;
  float met_x_jesdn = 0.0;
  float met_y_jesdn = 0.0;
  float met_x_jerup = 0.0;
  float met_y_jerup = 0.0;
  float met_x_jerdn = 0.0;
  float met_y_jerdn = 0.0;

  // loop over regular jets and jets that didn't make it into slimmedjets 
  // (CorrT1METJets)
  for (int jet_type(0); jet_type<2; jet_type++) {

    vector<float> jet_type_pt, jet_type_eta, jet_type_phi, jet_type_area;
    vector<float> jet_type_rawfactor, jet_type_muonfactor;
    int jet_type_size(0);
    if (jet_type==0) {
      jet_type_area = nano.Jet_area();
      jet_type_pt = nano.Jet_pt();
      jet_type_eta = nano.Jet_eta();
      jet_type_phi = nano.Jet_phi();
      jet_type_rawfactor = nano.Jet_rawFactor();
      jet_type_muonfactor = nano.Jet_muonSubtrFactor();
      jet_type_size = nano.nJet();
    }
    else {
      jet_type_area = nano.CorrT1METJet_area();
      jet_type_pt = nano.CorrT1METJet_rawPt();
      jet_type_eta = nano.CorrT1METJet_eta();
      jet_type_phi = nano.CorrT1METJet_phi();
      jet_type_rawfactor.resize(nano.nCorrT1METJet(),0.0);
      jet_type_muonfactor = nano.CorrT1METJet_muonSubtrFactor();
      jet_type_size = nano.nCorrT1METJet();
    }
    float rho = 0.0f;
    if (year <= 2018)
      rho = nano.fixedGridRhoFastjetAll();
    else
      rho = nano.Rho_fixedGridRhoFastjetAll();
    
    for (int ijet(0); ijet<jet_type_size; ++ijet) {
      // Need new JECs for low pT jets (Run 2) and all jets (Run 3)
      float jec_default = 1.0/(1.0-jet_type_rawfactor[ijet]);
      float jec_l1 = 1.0;
      float jec = jec_default;
      float jet_raw_pt = jet_type_pt[ijet]/jec_default;
      if (jet_type==1 || year > 2018) {
        jec = GetJEC(jet_type_area[ijet],jet_type_eta[ijet],jet_type_phi[ijet],
                     jet_raw_pt,rho,nano.run(),JECType::L1L2L3);
        if (year > 2018)
          jec_l1 = GetJEC(jet_type_area[ijet],jet_type_eta[ijet],
                          jet_type_phi[ijet],jet_raw_pt,rho,nano.run(),
                          JECType::L1);
      }
      float jet_l1_pt = jet_raw_pt*jec_l1;
      float jet_l1l2l3_pt = jet_raw_pt*jec;
      float jet_raw_pt_nomu = jet_raw_pt*(1.0-jet_type_muonfactor[ijet]);
      float jet_l1l2l3_pt_nomu = jet_raw_pt_nomu*jec;
      // jets below 15 GeV considered in unclustered energy
      if (jet_type == 1 && jet_l1l2l3_pt_nomu < 15.0f) continue;

      float indiv_jer_nm(1.0), indiv_jer_up(1.0), indiv_jer_dn(1.0);
      // calculate JER (smearing) factors for MC only
      // https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#Smearing_procedures
      // JSONS found at
      // https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/JME
      if (!isData) {

        float sigmajer = map_jermc_->evaluate({jet_type_eta[ijet],
                                               jet_l1l2l3_pt,rho});
        float sjer_nom = 1.0f;
        float sjer_up = 1.0f;
        float sjer_dn = 1.0f;
        if (year <= 2018) {
          sjer_nom = map_jersf_->evaluate({jet_type_eta[ijet],"nom"});
          sjer_up = map_jersf_->evaluate({jet_type_eta[ijet],"up"});
          sjer_dn = map_jersf_->evaluate({jet_type_eta[ijet],"down"});
        }
        else {
          sjer_nom = map_jersf_->evaluate({jet_type_eta[ijet],jet_l1l2l3_pt,
                                           "nom"});
          sjer_up = map_jersf_->evaluate({jet_type_eta[ijet],jet_l1l2l3_pt,
                                          "up"});
          sjer_dn = map_jersf_->evaluate({jet_type_eta[ijet],jet_l1l2l3_pt,
                                          "down"});
        }

        bool found_genjet = false;
        float mindr = 999.0f;
        for (int igen(0); igen<nano.nGenJet(); ++igen) {
          float dr = dR(jet_type_eta[ijet], nano.GenJet_eta()[igen], 
                        jet_type_phi[ijet], nano.GenJet_phi()[igen]);
          float dpt = jet_l1l2l3_pt-nano.GenJet_pt()[igen];
          if (dr < 0.2f && fabs(dpt) < 3.0f*sigmajer*jet_l1l2l3_pt) {
            if (dr > mindr) continue;
            mindr = dr;
            found_genjet = true;
            indiv_jer_nm = (1.0+(sjer_nom-1.0)*dpt/jet_l1l2l3_pt);
            indiv_jer_up = (1.0+(sjer_up-1.0)*dpt/jet_l1l2l3_pt);
            indiv_jer_dn = (1.0+(sjer_dn-1.0)*dpt/jet_l1l2l3_pt);
          }
        }

        if (!found_genjet) {
          float rand = rng_.Gaus(0,sigmajer);
          indiv_jer_nm = 1.0+rand*sqrt(std::max(sjer_nom*sjer_nom-1.0,0.0));
          indiv_jer_up = 1.0+rand*sqrt(std::max(sjer_up*sjer_up-1.0,0.0));
          indiv_jer_dn = 1.0+rand*sqrt(std::max(sjer_dn*sjer_dn-1.0,0.0));
        }
        //turn off smearing for non-gen-jets with pT<50 and 2.5<|eta|<3 in 2017 and onward
        //this fixes an issue with PU jets in the horn region
        //roughly modified strategy 2 from VBF SUSY: 
        //https://indico.cern.ch/event/1046356/contributions/4397877/attachments/2259227/3834282/BrendaFabelaEnriquez_VBFSUSY_METstudies_JERCMeeting_June7_2021.pdf#page=7
        if (year>=2017) {
          if (!found_genjet && fabs(jet_type_eta[ijet])>2.5f && 
              fabs(jet_type_eta[ijet])<3.0f && jet_type_pt[ijet]<50.0f) {
            indiv_jer_nm = 1.0;
            indiv_jer_up = 1.0;
            indiv_jer_dn = 1.0;
          }
        } 
      }

      float jet_factor = jet_l1l2l3_pt*indiv_jer_nm/jet_type_pt[ijet];
      float jes_unc = 0.0;

      //Save values for regular jets 
      if (jet_type==0) {
        jet_nm_factor.push_back(jet_factor);
        if (!isData) {
          //Following NanoAOD-tools, JES uncertainties are evaluated post-JER
          jes_unc = map_jes_->evaluate({jet_type_eta[ijet],
                                        jet_l1l2l3_pt*indiv_jer_nm});
          jer_up_factor.push_back(indiv_jer_up/indiv_jer_nm);
          jer_dn_factor.push_back(indiv_jer_dn/indiv_jer_nm);
          jes_up_factor.push_back(1.0+jes_unc);
          jes_dn_factor.push_back(1.0-jes_unc);
        }
      }
      //propagate corrections to MET
      //propagate corrections for jets with pt>15 GeV after subtracting muons
      //(below this threshold, propagate from unclustered energy), and skip
      //jets with >90% EM energy
      //https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETRun2Corrections
      float emef = 0.0;
      if (jet_type==0)
        emef = nano.Jet_neEmEF()[ijet]+nano.Jet_chEmEF()[ijet];
      float jet_cosphi = cos(jet_type_phi[ijet]);
      float jet_sinphi = sin(jet_type_phi[ijet]);
      // Run 2: starting from T1 corrected MET  (i.e. L2L3-L1 applied)
      // Run 3: Apply L2L3-L1 corrections
      if (jet_l1l2l3_pt_nomu*indiv_jer_nm > 15.0f && fabs(jet_type_eta[ijet])<5.2f 
          && emef < 0.9f) {
        float jet_nom_pt = jet_type_pt[ijet]*jet_factor;
        if (year <= 2018) {
          met_x_nom -= jet_cosphi*jet_type_pt[ijet]*(jet_factor-1.0);
          met_y_nom -= jet_sinphi*jet_type_pt[ijet]*(jet_factor-1.0);
        }
        else {
          met_x_nom -= jet_cosphi*(jet_nom_pt-jet_l1_pt);
          met_y_nom -= jet_sinphi*(jet_nom_pt-jet_l1_pt);
        }
        if (!isData) {
          met_x_jerup -= (jet_cosphi*jet_nom_pt*(indiv_jer_up/indiv_jer_nm
                                                 -1.0));
          met_y_jerup -= (jet_sinphi*jet_nom_pt*(indiv_jer_up/indiv_jer_nm
                                                 -1.0));
          met_x_jerdn -= (jet_cosphi*jet_nom_pt*(indiv_jer_dn/indiv_jer_nm
                                                 -1.0));
          met_y_jerdn -= (jet_sinphi*jet_nom_pt*(indiv_jer_dn/indiv_jer_nm
                                                 -1.0));
          met_x_jesup -= (jet_cosphi*jet_nom_pt*(jes_unc));
          met_y_jesup -= (jet_sinphi*jet_nom_pt*(jes_unc));
          met_x_jesdn -= (jet_cosphi*jet_nom_pt*(-1.0*jes_unc));
          met_y_jesdn -= (jet_sinphi*jet_nom_pt*(-1.0*jes_unc));
        }
      }
    }
  }
  met_x_jerup += met_x_nom;
  met_y_jerup += met_y_nom;
  met_x_jerdn += met_x_nom;
  met_y_jerdn += met_y_nom;
  met_x_jesup += met_x_nom;
  met_y_jesup += met_y_nom;
  met_x_jesdn += met_x_nom;
  met_y_jesdn += met_y_nom;

  pico.out_met() = sqrt(met_x_nom*met_x_nom+met_y_nom*met_y_nom);
  pico.out_met_phi() = atan2(met_y_nom, met_x_nom);
  if (!isData) {
    pico.out_sys_met().resize(4,0.0);
    pico.out_sys_met_phi().resize(4,0.0);
    pico.out_sys_met()[0] = sqrt(met_x_jerup*met_x_jerup
                                 +met_y_jerup*met_y_jerup);
    pico.out_sys_met()[1] = sqrt(met_x_jerdn*met_x_jerdn
                                 +met_y_jerdn*met_y_jerdn);
    pico.out_sys_met()[2] = sqrt(met_x_jesup*met_x_jesup
                                 +met_y_jesup*met_y_jesup);
    pico.out_sys_met()[3] = sqrt(met_x_jesdn*met_x_jesdn
                                 +met_y_jesdn*met_y_jesdn);
    pico.out_sys_met_phi()[0] = atan2(met_y_jerup, met_x_jerup);
    pico.out_sys_met_phi()[1] = atan2(met_y_jerdn, met_x_jerdn);
    pico.out_sys_met_phi()[2] = atan2(met_y_jesup, met_x_jesup);
    pico.out_sys_met_phi()[3] = atan2(met_y_jesdn, met_x_jesdn);
  }

  WriteMetVariations(nano, pico);
}

//writes other MET variables and non-jet MET uncertainties
void JetMetProducer::WriteMetVariations(nano_tree &nano, pico_tree &pico) {
  pico.out_met_calo()    = nano.CaloMET_pt();
  pico.out_met_tk()      = nano.TkMET_pt();
  pico.out_met_tk_phi()  = nano.TkMET_phi();
  pico.out_met_tru()     = nano.GenMET_pt();
  pico.out_met_tru_phi() = nano.GenMET_phi();
  pico.out_ht_isr_me()   = nano.LHE_HTIncoming();
  if (isData) {
    return;
  }

  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETRun2Corrections
  float met_x = pico.out_met()*cos(pico.out_met_phi());
  float met_y = pico.out_met()*sin(pico.out_met_phi());
  float met_x_unclup, met_y_unclup, met_x_uncldn, met_y_uncldn;
  if (year <= 2018) {
    met_x_unclup = met_x + nano.MET_MetUnclustEnUpDeltaX();
    met_y_unclup = met_y + nano.MET_MetUnclustEnUpDeltaX();
    met_x_uncldn = met_x - nano.MET_MetUnclustEnUpDeltaX();
    met_y_uncldn = met_y - nano.MET_MetUnclustEnUpDeltaX();
  }
  else {
    met_x_unclup = met_x + (nano.PuppiMET_ptUnclusteredUp()
                            *cos(nano.PuppiMET_phiUnclusteredUp())
                            -nano.PuppiMET_pt()*cos(nano.PuppiMET_phi()));
    met_y_unclup = met_y + (nano.PuppiMET_ptUnclusteredUp()
                            *sin(nano.PuppiMET_phiUnclusteredUp())
                            -nano.PuppiMET_pt()*sin(nano.PuppiMET_phi()));
    met_x_uncldn = met_x + (nano.PuppiMET_ptUnclusteredDown()
                            *cos(nano.PuppiMET_phiUnclusteredDown())
                            -nano.PuppiMET_pt()*cos(nano.PuppiMET_phi()));
    met_y_uncldn = met_y + (nano.PuppiMET_ptUnclusteredDown()
                            *sin(nano.PuppiMET_phiUnclusteredDown())
                            -nano.PuppiMET_pt()*sin(nano.PuppiMET_phi()));
  }
  float met_x_leptonphotonup = met_x;
  float met_y_leptonphotonup = met_y;
  float met_x_leptonphotondn = met_x;
  float met_y_leptonphotondn = met_y;
  for (int iel = 0; iel < pico.out_nel(); iel++) {
    if (pico.out_el_sig()[iel]) {
      float unc_factor = 0.006; //EB
      if (fabs(pico.out_el_eta()[iel])>1.5051f) unc_factor = 0.015; //EE
      met_x_leptonphotonup += 
          unc_factor*cos(pico.out_el_phi()[iel])*pico.out_el_pt()[iel];
      met_y_leptonphotonup += 
          unc_factor*sin(pico.out_el_phi()[iel])*pico.out_el_pt()[iel];
      met_x_leptonphotondn -= 
          unc_factor*cos(pico.out_el_phi()[iel])*pico.out_el_pt()[iel];
      met_y_leptonphotondn -= 
          unc_factor*sin(pico.out_el_phi()[iel])*pico.out_el_pt()[iel];
    }
  }
  for (int iph = 0; iph < pico.out_nphoton(); iph++) {
    if (pico.out_photon_sig()[iph]) {
      float unc_factor = 0.006; //EB
      if (fabs(pico.out_photon_eta()[iph])>1.5051f) unc_factor = 0.015; //EE
      met_x_leptonphotonup += 
          unc_factor*cos(pico.out_photon_phi()[iph])*pico.out_photon_pt()[iph];
      met_y_leptonphotonup += 
          unc_factor*sin(pico.out_photon_phi()[iph])*pico.out_photon_pt()[iph];
      met_x_leptonphotondn -= 
          unc_factor*cos(pico.out_photon_phi()[iph])*pico.out_photon_pt()[iph];
      met_y_leptonphotondn -= 
          unc_factor*sin(pico.out_photon_phi()[iph])*pico.out_photon_pt()[iph];
    }
  }
  for (int imu = 0; imu < pico.out_nmu(); imu++) {
    if (pico.out_mu_sig()[imu]) {
      float unc_factor = 0.002;
      met_x_leptonphotonup += 
          unc_factor*cos(pico.out_mu_phi()[imu])*pico.out_mu_pt()[imu];
      met_y_leptonphotonup += 
          unc_factor*sin(pico.out_mu_phi()[imu])*pico.out_mu_pt()[imu];
      met_x_leptonphotondn -= 
          unc_factor*cos(pico.out_mu_phi()[imu])*pico.out_mu_pt()[imu];
      met_y_leptonphotondn -= 
          unc_factor*sin(pico.out_mu_phi()[imu])*pico.out_mu_pt()[imu];
    }
  }
  pico.out_sys_met().push_back(sqrt(met_x_unclup*met_x_unclup+met_y_unclup*met_y_unclup));
  pico.out_sys_met().push_back(sqrt(met_x_uncldn*met_x_uncldn+met_y_uncldn*met_y_uncldn));
  pico.out_sys_met().push_back(sqrt(met_x_leptonphotonup*met_x_leptonphotonup
                               +met_y_leptonphotonup*met_y_leptonphotonup));
  pico.out_sys_met().push_back(sqrt(met_x_leptonphotondn*met_x_leptonphotondn
                               +met_y_leptonphotondn*met_y_leptonphotondn));
  pico.out_sys_met_phi().push_back(atan2(met_y_unclup, met_x_unclup));
  pico.out_sys_met_phi().push_back(atan2(met_y_uncldn, met_x_uncldn));
  pico.out_sys_met_phi().push_back(atan2(met_y_leptonphotonup, met_x_leptonphotonup));
  pico.out_sys_met_phi().push_back(atan2(met_y_leptonphotondn, met_x_leptonphotondn));

}

vector<int> JetMetProducer::WriteJetMet(nano_tree &nano, pico_tree &pico, 
    vector<int> jet_islep_nano_idx, 
    vector<int> jet_isvlep_nano_idx,  
    vector<int> jet_isphoton_nano_idx,
    const vector<float> &btag_wpts, 
    const vector<float> &btag_df_wpts,
    const vector<float> &btag_upt_wpts, 
    bool isFastsim, 
    bool isSignal,
    vector<HiggsConstructionVariables> &sys_higvars){
  vector<int> sig_jet_nano_idx;
  pico.out_njet() = 0; pico.out_ht() = 0; pico.out_ht5() = 0; 
  pico.out_nbl() = 0; pico.out_nbm() = 0; pico.out_nbt() = 0; 
  pico.out_nbdfl() = 0; pico.out_nbdfm() = 0; pico.out_nbdft() = 0;
  pico.out_nbuptl() = 0; pico.out_nbuptm() = 0; pico.out_nbuptt() = 0; 
  pico.out_ngenjet() = 0;
  pico.out_ismapvetoevt() = false;
  pico.out_ishemvetoevt() = false;
  //add smearing to jets and calculate uncertainties
  vector<float> Jet_pt, Jet_mass;
  vector<float> jet_nm_factor;
  vector<float> jer_up_factor, jer_dn_factor, jes_up_factor, jes_dn_factor;
  float MET_pt, MET_phi;
  if ((nanoaod_version+0.01) < 9) {
    getJetWithJEC(nano, isFastsim, Jet_pt, Jet_mass);
    getMETWithJEC(nano, year, isFastsim, MET_pt, MET_phi, is_preUL);
    jet_nm_factor.resize(Jet_pt.size(),1.0);
    jer_up_factor.resize(Jet_pt.size(),1.0);
    jer_dn_factor.resize(Jet_pt.size(),1.0);
    jes_up_factor.resize(Jet_pt.size(),1.0);
    jes_dn_factor.resize(Jet_pt.size(),1.0);
    met_producer.WriteMet(nano, pico, isFastsim, isSignal, is_preUL);
  }
  else  {
    PropagateJERC(nano, pico, jet_nm_factor, jer_up_factor, 
                  jer_dn_factor, jes_up_factor, jes_dn_factor);
    MET_pt = pico.out_met();
    MET_phi = pico.out_met_phi();
    for(int ijet(0); ijet<nano.nJet(); ++ijet) {
      float jet_pt = nano.Jet_pt()[ijet];
      float jet_mass = nano.Jet_mass()[ijet];
      Jet_pt.push_back(jet_pt*jet_nm_factor[ijet]);
      Jet_mass.push_back(jet_mass*jet_nm_factor[ijet]);
    }
  }
  vector<int> Jet_jetId;
  if ((nanoaod_version+0.01) < 13) getJetId(nano, nanoaod_version, Jet_jetId);
  // Jet_jetId is left uninitialized for v15, filled later
  vector<int> Jet_hadronFlavour;
  if (!isData) getJet_hadronFlavour(nano, nanoaod_version, Jet_hadronFlavour);
  vector<int> Jet_partonFlavour;
  if (!isData) getJet_partonFlavour(nano, nanoaod_version, Jet_partonFlavour);
  vector<int> Jet_genJetIdx;
  if (!isData) getJet_genJetIdx(nano, nanoaod_version, Jet_genJetIdx);
  vector<int> GenJet_partonFlavour;
  if (!isData) getGenJet_partonFlavour(nano, nanoaod_version, GenJet_partonFlavour);
  // calculate MHT; needed when saving jet info
  TLorentzVector mht_vec;
  for(int ijet(0); ijet<nano.nJet(); ++ijet){
    if (Jet_pt[ijet] > min_jet_pt) {
      TLorentzVector ijet_v4;
      ijet_v4.SetPtEtaPhiM(Jet_pt[ijet], nano.Jet_eta()[ijet], nano.Jet_phi()[ijet], Jet_mass[ijet]);
      mht_vec -= ijet_v4;
    }
  }
  pico.out_mht() = mht_vec.Pt();
  pico.out_mht_phi() = mht_vec.Phi();
  vector<vector<float>> sys_jet_met_dphi;
  if ((nanoaod_version+0.01) < 9) {
    if (isSignal) {
      pico.out_sys_njet().resize(4,0);
      pico.out_sys_nbl().resize(4,0);
      pico.out_sys_nbm().resize(4,0);
      pico.out_sys_nbt().resize(4,0);
      pico.out_sys_ht().resize(4,0.0);
      sys_jet_met_dphi.resize(4,vector<float>({}));
      sys_higvars.resize(4, HiggsConstructionVariables());
      pico.out_sys_low_dphi_met().resize(4,false);
    }
  }
  else if (!isData) { 
    pico.out_sys_njet().resize(4,0);
    pico.out_sys_nbl().resize(4,0);
    pico.out_sys_nbm().resize(4,0);
    pico.out_sys_nbt().resize(4,0);
  }

  //calculate jet quality variables first to order pico list
  vector<bool> jet_pass_jetidFix;//jet ID tight
  vector<bool> jet_pass_PUjetid;
  vector<bool> jet_islep; 
  vector<bool> jet_isvlep; 
  vector<bool> jet_isphoton; 
  vector<bool> jet_inetahornveto;
  vector<bool> jet_inhemveto;
  vector<bool> jet_invetomap;
  vector<bool> jet_isgood_min; //!lep, !pho, eta cut, ID(fix), eta horn
  vector<bool> jet_isgood; //same as isgood_min + HEM/Run3 veto + pT cut
  for(int ijet(0); ijet<nano.nJet(); ++ijet){
    float jet_abseta = fabs(nano.Jet_eta()[ijet]);
    // jetid applied to only fullsim and data
    // Unclear if this is the case for PU jet ID, treating it this way for now
    //Check if this is correct for R2
    if (isFastsim) {
      jet_pass_jetidFix.push_back(true);
      jet_pass_PUjetid.push_back(true);
    }
    else {
      if (year >=2022 && (nanoaod_version+0.01) > 13){//Run3 NanoAODv15
         bool tightId = map_jetid_tight_->evaluate({
                        nano.Jet_eta()[ijet], nano.Jet_chHEF()[ijet], 
                        nano.Jet_neHEF()[ijet], nano.Jet_chEmEF()[ijet], 
                        nano.Jet_neEmEF()[ijet], nano.Jet_muEF()[ijet], 
                        static_cast<int>(nano.Jet_chMultiplicity()[ijet]),
                        static_cast<int>(nano.Jet_neMultiplicity()[ijet]), 
                        static_cast<int>(nano.Jet_chMultiplicity()[ijet]) + 
                        static_cast<int>(nano.Jet_neMultiplicity()[ijet])
                        });
         bool tightIdLepVeto = map_jetid_tightlepveto_->evaluate({
                        nano.Jet_eta()[ijet], nano.Jet_chHEF()[ijet], 
                        nano.Jet_neHEF()[ijet], nano.Jet_chEmEF()[ijet], 
                        nano.Jet_neEmEF()[ijet], nano.Jet_muEF()[ijet], 
                        static_cast<int>(nano.Jet_chMultiplicity()[ijet]),
                        static_cast<int>(nano.Jet_neMultiplicity()[ijet]),
                        static_cast<int>(nano.Jet_chMultiplicity()[ijet]) +
                        static_cast<int>(nano.Jet_neMultiplicity()[ijet])
                        });
         int jetIdBits = 0;
         if (tightId == 1) jetIdBits+=2;
         if (tightIdLepVeto) jetIdBits+=4;
         Jet_jetId.push_back(jetIdBits);
         jet_pass_jetidFix.push_back((Jet_jetId[ijet] >= 1));//Check this one again. . . 
      } else if (year < 2022) {//Run2 NanoAODv9
        jet_pass_jetidFix.push_back((Jet_jetId[ijet] >= 1));
      }  else {//Run 3 NanoAODv12
        if(jet_abseta<=2.7f){
          jet_pass_jetidFix.push_back(Jet_jetId[ijet] & (0b010));
        } else if (jet_abseta>2.7f && jet_abseta<=3.0f){
          jet_pass_jetidFix.push_back((Jet_jetId[ijet] & (0b010)) 
              && (nano.Jet_neHEF()[ijet] < 0.99f));
        } else {
          jet_pass_jetidFix.push_back((Jet_jetId[ijet] & (0b010)) 
              && (nano.Jet_neEmEF()[ijet] < 0.4f));
        }
      }
      if(year == 2016){//PU jet id is applied using JEC corrected jets, but not JES/JER
                       //https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetIDUL
        jet_pass_PUjetid.push_back( (nano.Jet_puId()[ijet] & (1 << 0)) || nano.Jet_pt()[ijet] > 50.f);
      } else if (year <= 2018){
        jet_pass_PUjetid.push_back( (nano.Jet_puId()[ijet] & (1 << 2)) || nano.Jet_pt()[ijet] > 50.f);
      } else {
        jet_pass_PUjetid.push_back(true);
      }
    }
    // check overlap with signal leptons (or photons)
    // N.B. photon collection is not filled for Higgsino analysis, so there is 
    // no overlap removal!
    jet_islep.push_back(find(jet_islep_nano_idx.begin(), 
        jet_islep_nano_idx.end(), ijet) != jet_islep_nano_idx.end());
    jet_isvlep.push_back(find(jet_isvlep_nano_idx.begin(), 
        jet_isvlep_nano_idx.end(), ijet) != jet_isvlep_nano_idx.end());
    jet_isphoton.push_back(find(jet_isphoton_nano_idx.begin(), 
        jet_isphoton_nano_idx.end(), ijet) != jet_isphoton_nano_idx.end()); 
    // eta horn veto
    jet_inetahornveto.push_back(year >=2017 && Jet_pt[ijet] < 50.0f 
        && jet_abseta > 2.5f && jet_abseta < 3.0f); 
    // For studying jetmaps and HEM vetos. Don't include eta veto yet in order 
    // to remove anomalous met events.
    bool isgood_min = !jet_islep.back() && !jet_isphoton.back() 
        && (jet_abseta < max_jet_eta) && jet_pass_jetidFix.back();
    // 2018 HEM veto. JetMET POG gives tighter selection than ours except 
    // 15 GeV cut. Apply our selection, with lower pT cut for veto events.
    bool isvetohem  = false;
    if (year==2018 && isgood_min && jet_pass_PUjetid.back() && nano.Jet_pt()[ijet]>15.0f 
        && nano.Jet_eta()[ijet]>-3.2f && nano.Jet_eta()[ijet]<-1.3f 
        && nano.Jet_phi()[ijet]>-1.57f && nano.Jet_phi()[ijet]<-0.87f){
      if(isData && nano.run()>=319077) isvetohem = true;
      else if(!isData && (nano.event()%10000>3564)) isvetohem = true;
      //same approximate portion of MC as to match the data percentage. 
      pico.out_ishemvetoevt()=isvetohem;
    }
    jet_inhemveto.push_back(isvetohem);
    // Run 3 jet veto maps. JetMET POG gives tighter selection than ours except 
    // 15 GeV cut. Apply our selection, with lower pT cut for veto events.
    float veto = 0.0f; 
    float phicorr;
    //a dumb addition because sometimes jet phi is slightly larger than pi
    if(nano.Jet_phi()[ijet]>3.1415926f) phicorr = 3.1415926f;
    else if (nano.Jet_phi()[ijet]<-3.1415926f) phicorr = -3.1415926f;
    else phicorr = nano.Jet_phi()[ijet];

    if (year>=2022 && nano.Jet_pt()[ijet]>15.0f && jet_abseta<5.191f) 
      veto = map_jetveto_->evaluate({"jetvetomap", nano.Jet_eta()[ijet],
                                     phicorr});
    if(veto!=0.0 && isgood_min) {
      pico.out_ismapvetoevt()=true;
      jet_invetomap.push_back(true);
    }
    else {
      jet_invetomap.push_back(false);
    }
    jet_isgood_min.push_back(isgood_min && !jet_inetahornveto.back());
    jet_isgood.push_back(jet_isgood_min.back() && !jet_invetomap.back() 
                         && !jet_inhemveto.back() 
                         && (Jet_pt[ijet] > min_jet_pt) && jet_pass_PUjetid.back());
  }
  //determine ordering based on isgood and pt
  vector<NanoOrderEntry> nano_entries;
  vector<int> ordered_nano_indices;
  for(int ijet(0); ijet<nano.nJet(); ++ijet){
    NanoOrderEntry nano_entry;
    nano_entry.nano_idx = ijet;
    nano_entry.pt = Jet_pt[ijet];
    nano_entry.is_sig = jet_isgood[ijet];
    nano_entries.push_back(nano_entry);
  }
  sort(nano_entries.begin(),nano_entries.end(), 
      [](NanoOrderEntry a, NanoOrderEntry b) {
        if (a.is_sig && !b.is_sig) return true;
        if (b.is_sig && !a.is_sig) return false;
        return (a.pt>b.pt);
      });
  for (NanoOrderEntry nano_entry : nano_entries)
    ordered_nano_indices.push_back(nano_entry.nano_idx);
  // saving jet info on all jets passing pt cut with any variation
  for(int ijet : ordered_nano_indices) {
    if (verbose) cout<<"Jet "<<ijet<<": pt = "<<setw(10)<<Jet_pt[ijet]
                                    <<" eta = "<<setw(10)<<nano.Jet_eta()[ijet]
                                    <<" phi = "<<setw(10)<<nano.Jet_phi()[ijet]
                                    <<" m = "<<setw(10)<<Jet_mass[ijet]
                                    <<endl;
    bool isgood_nopt = jet_isgood_min[ijet] && !jet_invetomap[ijet] 
                       && !jet_inhemveto[ijet];

    //sys_jetvar convention: [0] JER up, [1] JER down, [2] JEC up, [3] JEC down
    //for now, only save sys_ variables
    //old code for preUL Nanos; see GetJetUncertainties for newer Nano versions
    if (isSignal && (nanoaod_version+0.01)<9) {
      if (nano.Jet_pt_jerUp()[ijet] > min_jet_pt) {
        if (isgood_nopt) {
          pico.out_sys_njet()[0]++;
          pico.out_sys_ht()[0] += nano.Jet_pt_jerUp()[ijet];
          if (nano.Jet_btagDeepB()[ijet] > btag_wpts[0]) pico.out_sys_nbl()[0]++; 
          if (nano.Jet_btagDeepB()[ijet] > btag_wpts[1]) pico.out_sys_nbm()[0]++; 
          if (nano.Jet_btagDeepB()[ijet] > btag_wpts[2]) pico.out_sys_nbt()[0]++;
          sys_higvars[0].jet_deepcsv.push_back(nano.Jet_btagDeepB()[ijet]);
          TLorentzVector lv;
          lv.SetPtEtaPhiM(nano.Jet_pt_jerUp()[ijet], nano.Jet_eta()[ijet],
                          nano.Jet_phi()[ijet], nano.Jet_mass_jerUp()[ijet]);
          sys_higvars[0].jet_lv.push_back(lv);
        }
        if (fabs(nano.Jet_eta()[ijet]) < 2.4f)
          sys_jet_met_dphi.at(0).push_back(DeltaPhi(nano.Jet_phi()[ijet], pico.out_sys_met_phi()[0]));
      }
      if (nano.Jet_pt_jerDown()[ijet] > min_jet_pt) {
        if (isgood_nopt) {
          pico.out_sys_njet()[1]++;
          pico.out_sys_ht()[1] += nano.Jet_pt_jerDown()[ijet];
          if (nano.Jet_btagDeepB()[ijet] > btag_wpts[0]) pico.out_sys_nbl()[1]++; 
          if (nano.Jet_btagDeepB()[ijet] > btag_wpts[1]) pico.out_sys_nbm()[1]++; 
          if (nano.Jet_btagDeepB()[ijet] > btag_wpts[2]) pico.out_sys_nbt()[1]++;
          sys_higvars[1].jet_deepcsv.push_back(nano.Jet_btagDeepB()[ijet]);
          TLorentzVector lv;
          lv.SetPtEtaPhiM(nano.Jet_pt_jerDown()[ijet], nano.Jet_eta()[ijet],
                          nano.Jet_phi()[ijet], nano.Jet_mass_jerDown()[ijet]);
          sys_higvars[1].jet_lv.push_back(lv);
        }
        if (fabs(nano.Jet_eta()[ijet]) < 2.4)
          sys_jet_met_dphi.at(1).push_back(DeltaPhi(nano.Jet_phi()[ijet], pico.out_sys_met_phi()[1]));
      }
      if (nano.Jet_pt_jesTotalUp()[ijet] > min_jet_pt) {
        if (isgood_nopt) {
          pico.out_sys_njet()[2]++;
          pico.out_sys_ht()[2] += nano.Jet_pt_jesTotalUp()[ijet];
          if (nano.Jet_btagDeepB()[ijet] > btag_wpts[0]) pico.out_sys_nbl()[2]++; 
          if (nano.Jet_btagDeepB()[ijet] > btag_wpts[1]) pico.out_sys_nbm()[2]++; 
          if (nano.Jet_btagDeepB()[ijet] > btag_wpts[2]) pico.out_sys_nbt()[2]++;
          sys_higvars[2].jet_deepcsv.push_back(nano.Jet_btagDeepB()[ijet]);
          TLorentzVector lv;
          lv.SetPtEtaPhiM(nano.Jet_pt_jesTotalUp()[ijet], nano.Jet_eta()[ijet],
                          nano.Jet_phi()[ijet], nano.Jet_mass_jesTotalUp()[ijet]);
          sys_higvars[2].jet_lv.push_back(lv);
        }
        if (fabs(nano.Jet_eta()[ijet]) < 2.4f)
          sys_jet_met_dphi.at(2).push_back(DeltaPhi(nano.Jet_phi()[ijet], pico.out_sys_met_phi()[2]));
      }
      if (nano.Jet_pt_jesTotalDown()[ijet] > min_jet_pt) {
        if (isgood_nopt) {
          pico.out_sys_njet()[3]++;
          pico.out_sys_ht()[3] += nano.Jet_pt_jesTotalDown()[ijet];
          if (nano.Jet_btagDeepB()[ijet] > btag_wpts[0]) pico.out_sys_nbl()[3]++; 
          if (nano.Jet_btagDeepB()[ijet] > btag_wpts[1]) pico.out_sys_nbm()[3]++; 
          if (nano.Jet_btagDeepB()[ijet] > btag_wpts[2]) pico.out_sys_nbt()[3]++;
          sys_higvars[3].jet_deepcsv.push_back(nano.Jet_btagDeepB()[ijet]);
          TLorentzVector lv;
          lv.SetPtEtaPhiM(nano.Jet_pt_jesTotalDown()[ijet], nano.Jet_eta()[ijet],
                          nano.Jet_phi()[ijet], nano.Jet_mass_jesTotalDown()[ijet]);
          sys_higvars[3].jet_lv.push_back(lv);
        }
        if (fabs(nano.Jet_eta()[ijet]) < 2.4f)
          sys_jet_met_dphi.at(3).push_back(DeltaPhi(nano.Jet_phi()[ijet], pico.out_sys_met_phi()[3]));
      }
    }
    if (isData) {
      if (Jet_pt[ijet] <= min_jet_pt) continue;
    }
    else {
      if ((nanoaod_version+0.01)<9) {
        if (Jet_pt[ijet] <= min_jet_pt && 
            nano.Jet_pt_jerUp()[ijet]<=min_jet_pt &&
            nano.Jet_pt_jerDown()[ijet]<=min_jet_pt &&
            nano.Jet_pt_jesTotalUp()[ijet]<=min_jet_pt &&
            nano.Jet_pt_jesTotalDown()[ijet]<=min_jet_pt) continue;
      }
      else {
        if (Jet_pt[ijet] <= min_jet_pt && 
            (Jet_pt[ijet]*jer_up_factor[ijet]) <= min_jet_pt &&
            (Jet_pt[ijet]*jer_dn_factor[ijet]) <= min_jet_pt &&
            (Jet_pt[ijet]*jes_up_factor[ijet]) <= min_jet_pt &&
            (Jet_pt[ijet]*jes_dn_factor[ijet]) <= min_jet_pt) continue;
      }
    }
    if (!isData && (nanoaod_version+0.01)>9 && isgood_nopt) {
      if ((Jet_pt[ijet]*jer_up_factor[ijet]) > min_jet_pt) {
        pico.out_sys_njet()[0]++;
        if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[0]) 
          pico.out_sys_nbl()[0]++; 
        if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[1]) 
          pico.out_sys_nbm()[0]++; 
        if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[2]) 
          pico.out_sys_nbt()[0]++; 
      }
      if ((Jet_pt[ijet]*jer_dn_factor[ijet]) > min_jet_pt) {
        pico.out_sys_njet()[1]++;
        if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[0]) 
          pico.out_sys_nbl()[1]++; 
        if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[1]) 
          pico.out_sys_nbm()[1]++; 
        if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[2]) 
          pico.out_sys_nbt()[1]++; 
      }
      if ((Jet_pt[ijet]*jes_up_factor[ijet]) > min_jet_pt) {
        pico.out_sys_njet()[2]++;
        if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[0]) 
          pico.out_sys_nbl()[2]++; 
        if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[1]) 
          pico.out_sys_nbm()[2]++; 
        if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[2]) 
          pico.out_sys_nbt()[2]++; 
      }
      if ((Jet_pt[ijet]*jes_dn_factor[ijet]) > min_jet_pt) {
        pico.out_sys_njet()[3]++;
        if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[0]) 
          pico.out_sys_nbl()[3]++; 
        if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[1]) 
          pico.out_sys_nbm()[3]++; 
        if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[2]) 
          pico.out_sys_nbt()[3]++; 
      }
    }
    switch(year) {
      case 2016:
      case 2017:
      case 2018:
        pico.out_jet_pt().push_back(Jet_pt[ijet]);
        pico.out_jet_nanopt().push_back(nano.Jet_pt()[ijet]);
        pico.out_jet_eta().push_back(nano.Jet_eta()[ijet]);
        pico.out_jet_phi().push_back(nano.Jet_phi()[ijet]);
        pico.out_jet_m().push_back(Jet_mass[ijet]);
        pico.out_jet_breg_corr().push_back(nano.Jet_bRegCorr()[ijet]);
        pico.out_jet_breg_res().push_back(nano.Jet_bRegRes()[ijet]);
        pico.out_jet_deepcsv().push_back(nano.Jet_btagDeepB()[ijet]);
        pico.out_jet_deepflav().push_back(nano.Jet_btagDeepFlavB()[ijet]);
        pico.out_jet_ne_emef().push_back(nano.Jet_neEmEF()[ijet]);
        pico.out_jet_qgl().push_back(nano.Jet_qgl()[ijet]);
        pico.out_jet_islep().push_back(jet_islep[ijet]);
        pico.out_jet_isvlep().push_back(jet_isvlep[ijet]);
        pico.out_jet_isphoton().push_back(jet_isphoton[ijet]);
        pico.out_jet_isgood().push_back(jet_isgood[ijet]);
        pico.out_jet_isgood_min().push_back(jet_isgood_min[ijet]);
        pico.out_jet_isvetomap().push_back(jet_invetomap[ijet]);
        pico.out_jet_isvetohem().push_back(jet_inhemveto[ijet]);
        pico.out_jet_id().push_back(Jet_jetId[ijet]);
        pico.out_jet_mht_dphi().push_back(DeltaPhi(nano.Jet_phi()[ijet], mht_vec.Phi()));
        pico.out_jet_met_dphi().push_back(DeltaPhi(nano.Jet_phi()[ijet], MET_phi));
        pico.out_jet_puid().push_back(nano.Jet_puId()[ijet]);
        pico.out_jet_puid_disc().push_back(nano.Jet_puIdDisc()[ijet]);
        pico.out_jet_puid_pass().push_back(jet_pass_PUjetid[ijet]);
        if (!isData && isSignal) {
          pico.out_sys_jet_pt_jesup().push_back(
              Jet_pt[ijet]*jes_up_factor[ijet]);
          pico.out_sys_jet_pt_jesdn().push_back(
              Jet_pt[ijet]*jes_dn_factor[ijet]);
          pico.out_sys_jet_pt_jerup().push_back(
              Jet_pt[ijet]*jer_up_factor[ijet]);
          pico.out_sys_jet_pt_jerdn().push_back(
              Jet_pt[ijet]*jer_dn_factor[ijet]);
          pico.out_sys_jet_m_jesup().push_back(
              Jet_mass[ijet]*jes_up_factor[ijet]);
          pico.out_sys_jet_m_jesdn().push_back(
              Jet_mass[ijet]*jes_dn_factor[ijet]);
          pico.out_sys_jet_m_jerup().push_back(
              Jet_mass[ijet]*jer_up_factor[ijet]);
          pico.out_sys_jet_m_jerdn().push_back(
              Jet_mass[ijet]*jer_dn_factor[ijet]);
          pico.out_sys_jet_isgood_jesup().push_back(isgood_nopt
              && (Jet_pt[ijet]*jes_up_factor[ijet] > min_jet_pt));
          pico.out_sys_jet_isgood_jesdn().push_back(isgood_nopt
              && (Jet_pt[ijet]*jes_dn_factor[ijet] > min_jet_pt));
          pico.out_sys_jet_isgood_jerup().push_back(isgood_nopt
              && (Jet_pt[ijet]*jer_up_factor[ijet] > min_jet_pt));
          pico.out_sys_jet_isgood_jerdn().push_back(isgood_nopt
              && (Jet_pt[ijet]*jer_dn_factor[ijet] > min_jet_pt));
        }
        break;
      case 2022:
      case 2023:
      case 2024:
      case 2025:
        pico.out_jet_pt().push_back(Jet_pt[ijet]);
        pico.out_jet_nanopt().push_back(nano.Jet_pt()[ijet]);
        pico.out_jet_eta().push_back(nano.Jet_eta()[ijet]);
        pico.out_jet_phi().push_back(nano.Jet_phi()[ijet]);
        pico.out_jet_m().push_back(Jet_mass[ijet]);
        //pico.out_jet_breg_corr().push_back(nano.Jet_bRegCorr()[ijet]);
        //pico.out_jet_breg_res().push_back(nano.Jet_bRegRes()[ijet]);
        if (nanoaod_version < 11.89) pico.out_jet_deepcsv().push_back(nano.Jet_btagDeepB()[ijet]);
        if (nanoaod_version > 11.5) pico.out_jet_btagpnetb().push_back(nano.Jet_btagPNetB()[ijet]);
        if (nanoaod_version > 11.95 && nanoaod_version < 14.9) 
           pico.out_jet_btagak4b().push_back(nano.Jet_btagRobustParTAK4B()[ijet]);
        if(nanoaod_version+0.01 > 15)
           pico.out_jet_btaguptb().push_back(nano.Jet_btagUParTAK4B()[ijet]);
        pico.out_jet_deepflav().push_back(nano.Jet_btagDeepFlavB()[ijet]);
        pico.out_jet_ne_emef().push_back(nano.Jet_neEmEF()[ijet]);
        //pico.out_jet_qgl().push_back(nano.Jet_qgl()[ijet]);
        pico.out_jet_islep().push_back(jet_islep[ijet]);
        pico.out_jet_isvlep().push_back(jet_isvlep[ijet]);
        pico.out_jet_isphoton().push_back(jet_isphoton[ijet]);
        pico.out_jet_isgood().push_back(jet_isgood[ijet]);
        pico.out_jet_isgood_min().push_back(jet_isgood_min[ijet]);
        pico.out_jet_isvetomap().push_back(jet_invetomap[ijet]);
        pico.out_jet_isvetohem().push_back(jet_inhemveto[ijet]);
        pico.out_jet_id().push_back(Jet_jetId[ijet]);
        if (nanoaod_version >= 15) 
           pico.out_jet_puid_disc().push_back(nano.Jet_puIdDisc()[ijet]);
        pico.out_jet_mht_dphi().push_back(DeltaPhi(nano.Jet_phi()[ijet], mht_vec.Phi()));
        pico.out_jet_met_dphi().push_back(DeltaPhi(nano.Jet_phi()[ijet], MET_phi));
        pico.out_jet_puid_pass().push_back(true);
        //pico.out_jet_puid().push_back(nano.Jet_puId()[ijet]);
        //pico.out_jet_puid_disc().push_back(nano.Jet_puIdDisc()[ijet]);
        if (!isData) {
          pico.out_sys_jet_pt_jesup().push_back(
              Jet_pt[ijet]*jes_up_factor[ijet]);
          pico.out_sys_jet_pt_jesdn().push_back(
              Jet_pt[ijet]*jes_dn_factor[ijet]);
          pico.out_sys_jet_pt_jerup().push_back(
              Jet_pt[ijet]*jer_up_factor[ijet]);
          pico.out_sys_jet_pt_jerdn().push_back(
              Jet_pt[ijet]*jer_dn_factor[ijet]);
          pico.out_sys_jet_m_jesup().push_back(
              Jet_mass[ijet]*jes_up_factor[ijet]);
          pico.out_sys_jet_m_jesdn().push_back(
              Jet_mass[ijet]*jes_dn_factor[ijet]);
          pico.out_sys_jet_m_jerup().push_back(
              Jet_mass[ijet]*jer_up_factor[ijet]);
          pico.out_sys_jet_m_jerdn().push_back(
              Jet_mass[ijet]*jer_dn_factor[ijet]);
          pico.out_sys_jet_isgood_jesup().push_back(isgood_nopt
              && (Jet_pt[ijet]*jes_up_factor[ijet] > min_jet_pt));
          pico.out_sys_jet_isgood_jesdn().push_back(isgood_nopt
              && (Jet_pt[ijet]*jes_dn_factor[ijet] > min_jet_pt));
          pico.out_sys_jet_isgood_jerup().push_back(isgood_nopt
              && (Jet_pt[ijet]*jer_up_factor[ijet] > min_jet_pt));
          pico.out_sys_jet_isgood_jerdn().push_back(isgood_nopt
              && (Jet_pt[ijet]*jer_dn_factor[ijet] > min_jet_pt));
        }
        break;
      default:
        std::cout<<"Need code for new year in getZGammaJetBr in jetmet_producer.cpp"<<endl;
        exit(1);
    }
    if (!isData && Jet_hadronFlavour.size() > 0) {
      pico.out_jet_hflavor().push_back(Jet_hadronFlavour[ijet]);
    pico.out_jet_pflavor().push_back(Jet_partonFlavour[ijet]);
      pico.out_jet_genjet_idx().push_back(Jet_genJetIdx[ijet]);
    }
    // will be overwritten with the overlapping fat jet index, if such exists, in WriteFatJets
    pico.out_jet_fjet_idx().push_back(-999);
    //the jets for the higgs pair with smallest dm will be set to true in hig_producer
    pico.out_jet_h1d().push_back(false);
    pico.out_jet_h2d().push_back(false);

    if (!jet_islep[ijet] && !jet_isphoton[ijet]) 
      pico.out_ht5() += Jet_pt[ijet];
    if (jet_isgood[ijet]) {
      sig_jet_nano_idx.push_back(ijet);
      pico.out_njet()++;
      pico.out_ht() += Jet_pt[ijet];
      if ((nanoaod_version+0.01) < 12){
        if (nano.Jet_btagDeepB()[ijet] > btag_wpts[0]) pico.out_nbl()++; 
        if (nano.Jet_btagDeepB()[ijet] > btag_wpts[1]) pico.out_nbm()++; 
        if (nano.Jet_btagDeepB()[ijet] > btag_wpts[2]) pico.out_nbt()++;
      } else {
        if (nano.Jet_btagPNetB()[ijet] > btag_wpts[0]) pico.out_nbl()++;
        if (nano.Jet_btagPNetB()[ijet] > btag_wpts[1]) pico.out_nbm()++;
        if (nano.Jet_btagPNetB()[ijet] > btag_wpts[2]) pico.out_nbt()++;
      }
      if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[0]) pico.out_nbdfl()++; 
      if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[1]) pico.out_nbdfm()++; 
      if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[2]) pico.out_nbdft()++;
      if ((nanoaod_version+0.01)>15){
        if (nano.Jet_btagUParTAK4B()[ijet] > btag_upt_wpts[0]) pico.out_nbuptl()++;
        if (nano.Jet_btagUParTAK4B()[ijet] > btag_upt_wpts[1]) pico.out_nbuptm()++;
        if (nano.Jet_btagUParTAK4B()[ijet] > btag_upt_wpts[2]) pico.out_nbuptt()++;
      } 
    }
  } // end jet loop
  pico.out_low_dphi_mht_e5() = false;
  pico.out_low_dphi_met_e5() = false;
  for (unsigned ijet(0); ijet<pico.out_jet_mht_dphi().size(); ijet++){
    float cut_ = ijet<=1 ? 0.5 : 0.3;
    if (pico.out_jet_mht_dphi()[ijet]<=cut_) pico.out_low_dphi_mht_e5() = true;
    if (pico.out_jet_met_dphi()[ijet]<=cut_) pico.out_low_dphi_met_e5() = true;
    if (ijet==3) break;
  }

  pico.out_low_dphi_mht() = false;
  pico.out_low_dphi_met() = false;
  for (unsigned ijet(0), ijet_e24(0); ijet<pico.out_jet_mht_dphi().size(); ijet++){
    if (fabs(pico.out_jet_eta()[ijet]) > max_jet_eta) continue;
    float cut_ = ijet_e24<=1 ? 0.5 : 0.3;
    if (pico.out_jet_mht_dphi()[ijet]<=cut_) pico.out_low_dphi_mht() = true;
    if (pico.out_jet_met_dphi()[ijet]<=cut_) pico.out_low_dphi_met() = true;
    if (ijet_e24==3) break;
    ijet_e24++;
  }

  if (isSignal && (nanoaod_version+0.01)<9) {
    for (unsigned ijec(0); ijec < 4; ijec++) {
      for (unsigned ijet(0); ijet < sys_jet_met_dphi[ijec].size(); ijet++) {
        float cut_ = ijet<=1 ? 0.5 : 0.3;
        if (sys_jet_met_dphi[ijec][ijet] <= cut_) pico.out_sys_low_dphi_met()[ijec] = true;
        if (ijet==3) break;
      }
    }
  }
  // Saving GenJet information
  if (!isData) {
    pico.out_ngenjet() = nano.nGenJet();
    for(int ijet(0); ijet<nano.nGenJet(); ++ijet){
      pico.out_genjet_pt().push_back(nano.GenJet_pt()[ijet]);
      pico.out_genjet_eta().push_back(nano.GenJet_eta()[ijet]);
      pico.out_genjet_phi().push_back(nano.GenJet_phi()[ijet]);
      pico.out_genjet_m().push_back(nano.GenJet_mass()[ijet]);
      pico.out_genjet_pflavor().push_back(GenJet_partonFlavour[ijet]);
      pico.out_genjet_hflavor().push_back(int(nano.GenJet_hadronFlavour()[ijet]));
    }
  }
  if (verbose) cout<<"Done with AK4 jets"<<endl;
  return sig_jet_nano_idx;
}

void JetMetProducer::WriteFatJets(nano_tree &nano, pico_tree &pico){
  pico.out_nfjet() = 0; 
  vector<float> FatJet_btagDDBvL;
  vector<float> FatJet_particleNetWithMass_WvsQCD;
  vector<float> FatJet_particleNetWithMass_ZvsQCD;
  vector<float> FatJet_particleNetWithMass_TvsQCD;
  vector<float> FatJet_particleNet_mass;
  vector<int> FatJet_subJetIdx1;
  vector<int> FatJet_subJetIdx2;
  if(nanoaod_version < 13.1) getFatJet_btagDDBvL(nano, nanoaod_version, FatJet_btagDDBvL);
  getFatJet_subJetIdx1(nano, nanoaod_version, FatJet_subJetIdx1);
  getFatJet_subJetIdx2(nano, nanoaod_version, FatJet_subJetIdx2);
  if (nanoaod_version+0.01 > 9) {
    getFatJet_particleNetWithMass_WvsQCD(nano, nanoaod_version, FatJet_particleNetWithMass_WvsQCD);
    getFatJet_particleNetWithMass_ZvsQCD(nano, nanoaod_version, FatJet_particleNetWithMass_ZvsQCD);
    getFatJet_particleNetWithMass_TvsQCD(nano, nanoaod_version, FatJet_particleNetWithMass_TvsQCD);
    getFatJet_particleNet_mass(nano, nanoaod_version, FatJet_particleNet_mass);
  }

  for(int ifjet(0); ifjet<nano.nFatJet(); ++ifjet){
    if (verbose) cout<<"FatJet "<<ifjet<<": pt = "<<setw(10)<<nano.FatJet_pt()[ifjet]
                                       <<" eta = "<<setw(10)<<nano.FatJet_eta()[ifjet]
                                       <<" phi = "<<setw(10)<<nano.FatJet_phi()[ifjet]
                                       <<" m = "<<setw(10)<<nano.FatJet_mass()[ifjet]
                                       <<endl;

    pico.out_fjet_pt().push_back(nano.FatJet_pt()[ifjet]);
    pico.out_fjet_eta().push_back(nano.FatJet_eta()[ifjet]);
    pico.out_fjet_phi().push_back(nano.FatJet_phi()[ifjet]);
    pico.out_fjet_m().push_back(nano.FatJet_mass()[ifjet]);
    pico.out_fjet_msoftdrop().push_back(nano.FatJet_msoftdrop()[ifjet]);
    // Mass-decorrelated Deep Double B, H->bb vs QCD discriminator, endorsed by BTV
    if(nanoaod_version < 13.1) pico.out_fjet_deep_md_hbb_btv().push_back(FatJet_btagDDBvL[ifjet]);
    if(nanoaod_version < 13.1) pico.out_fjet_mva_hbb_btv().push_back(nano.FatJet_btagHbb()[ifjet]);
    if (nanoaod_version+0.01 < 11.9) {
      // Mass-decorrelated DeepAK8, H->bb vs QCD discriminator, endorsed by JME
      pico.out_fjet_deep_md_hbb_jme().push_back(nano.FatJet_deepTagMD_HbbvsQCD()[ifjet]);
      pico.out_fjet_deep_md_tvsqcd().push_back(nano.FatJet_deepTagMD_TvsQCD()[ifjet]);
      pico.out_fjet_deep_tvsqcd().push_back(nano.FatJet_deepTag_TvsQCD()[ifjet]);
    }
    if (nanoaod_version+0.01 > 9) {
      pico.out_fjet_pnet_wtag().push_back(FatJet_particleNetWithMass_WvsQCD[ifjet]);
      pico.out_fjet_pnet_ztag().push_back(FatJet_particleNetWithMass_ZvsQCD[ifjet]);
      pico.out_fjet_pnet_ttag().push_back(FatJet_particleNetWithMass_TvsQCD[ifjet]);
      pico.out_fjet_pnet_m().push_back(FatJet_particleNet_mass[ifjet]);
    }
    pico.out_fjet_subjet_idx1().push_back(FatJet_subJetIdx1[ifjet]);
    pico.out_fjet_subjet_idx2().push_back(FatJet_subJetIdx2[ifjet]);
    for(unsigned ijet(0); ijet<pico.out_jet_pt().size(); ++ijet){
      if (dR(pico.out_jet_eta()[ijet], nano.FatJet_eta()[ifjet], pico.out_jet_phi()[ijet], nano.FatJet_phi()[ifjet])<0.8f)
        pico.out_jet_fjet_idx()[ijet] = ifjet;
    }
    pico.out_nfjet()++;
  }
  if (verbose) cout<<"Done with fat jets"<<endl;
  return;
}

void JetMetProducer::WriteSubJets(nano_tree &nano, pico_tree &pico){
  pico.out_nsubfjet() = 0; 
  set<int> matched_ak4_jets;
  for(int isubj(0); isubj<nano.nSubJet(); ++isubj){
    pico.out_subfjet_pt().push_back(nano.SubJet_pt()[isubj]);
    pico.out_subfjet_eta().push_back(nano.SubJet_eta()[isubj]);
    pico.out_subfjet_phi().push_back(nano.SubJet_phi()[isubj]);
    pico.out_subfjet_m().push_back(nano.SubJet_mass()[isubj]);
    if (nanoaod_version < 13.1) pico.out_subfjet_deepcsv().push_back(nano.SubJet_btagDeepB()[isubj]);
    pico.out_subfjet_raw_factor().push_back(nano.SubJet_rawFactor()[isubj]);
    pico.out_subfjet_tau1().push_back(nano.SubJet_tau1()[isubj]);
    pico.out_subfjet_tau2().push_back(nano.SubJet_tau2()[isubj]);
    pico.out_subfjet_tau3().push_back(nano.SubJet_tau3()[isubj]);
    pico.out_subfjet_tau4().push_back(nano.SubJet_tau4()[isubj]);
    //match to highest pT ak4 jet that has not already been matched to a previous subjet
    float mindr(999.); int closest_jet(-1);
    for(unsigned ijet(0); ijet<pico.out_jet_pt().size(); ++ijet){
      float idr = dR(pico.out_jet_eta()[ijet], nano.SubJet_eta()[isubj], pico.out_jet_phi()[ijet], nano.SubJet_phi()[isubj]);
      if (idr<mindr){
        mindr = idr;
        if (mindr<0.4f && matched_ak4_jets.find(ijet)==matched_ak4_jets.end()) {
          closest_jet = ijet;
          matched_ak4_jets.insert(ijet);
          break; 
        }
      }
    }
    pico.out_subfjet_jet_idx().push_back(closest_jet);

    if (verbose) cout<<"SubJet "<<isubj<<": pt = "<<setw(10)<<nano.SubJet_pt()[isubj]
                                       <<" eta = "<<setw(10)<<nano.SubJet_eta()[isubj]
                                       <<" phi = "<<setw(10)<<nano.SubJet_phi()[isubj]
                                       <<" m = "<<setw(10)<<nano.SubJet_mass()[isubj]
                                       <<" ijet = "<<closest_jet<<endl;
    pico.out_nsubfjet()++;
  }
  if (verbose) cout<<"Done with subjets"<<endl;
  return;
}

void JetMetProducer::WriteJetSystemPt(nano_tree &nano, pico_tree &pico, 
                                      vector<int> &sig_jet_nano_idx, const float &btag_wpt, bool isFastsim){

  vector<float> Jet_pt, Jet_mass;
  getJetWithJEC(nano, isFastsim, Jet_pt, Jet_mass);

  TLorentzVector jetsys_v4, jetsys_nob_v4;
  int njet_nob(0);
  for (auto &idx: sig_jet_nano_idx) {
    TLorentzVector ijet_v4;
    ijet_v4.SetPtEtaPhiM(Jet_pt[idx], nano.Jet_eta()[idx], nano.Jet_phi()[idx], Jet_mass[idx]);
    jetsys_v4 += ijet_v4;

    if (nanoaod_version+0.1 < 12) {
      if (nano.Jet_btagDeepB()[idx] <= btag_wpt){
        njet_nob++;
        jetsys_nob_v4 += ijet_v4;
      }
    } else {
      if (nano.Jet_btagPNetB()[idx] <= btag_wpt){
        njet_nob++;
        jetsys_nob_v4 += ijet_v4;     
      }
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


