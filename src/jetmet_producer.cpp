#include "jetmet_producer.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip> 
#include <iostream>
#include <vector>

#include "TLorentzVector.h"

#include "correction.hpp"
#include "hig_producer.hpp"
#include "met_producer.hpp"
#include "utilities.hpp"

using namespace std;

JetMetProducer::JetMetProducer(int year_, float nanoaod_version_, 
                               float min_jet_pt_, float max_jet_eta_, 
                               bool isData_, bool preVFP, bool is_preUL_, 
                               bool verbose_) : 
    met_producer(MetProducer(year_, isData_, is_preUL_)) {
  year = year_;
  isData = isData_;
  is_preUL = is_preUL_;
  verbose = verbose_;
  min_jet_pt = min_jet_pt_;
  max_jet_eta = max_jet_eta_;
  nanoaod_version = nanoaod_version_;
  rng_ = TRandom3(0);
  if (year==2016 && preVFP) {
    cs_jerc_ = correction::CorrectionSet::from_file("data/zgamma/2016preVFP_UL/jet_jerc_2016apv.json");
    //despite strange name, this map does have JES variations and should be evaluated w.r.t. corrected jet pt
    map_jes_ = cs_jerc_->at("Summer19UL16APV_V7_MC_Total_AK4PFchs");
    map_jersf_ = cs_jerc_->at("Summer20UL16APV_JRV3_MC_ScaleFactor_AK4PFchs");
    map_jermc_ = cs_jerc_->at("Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs");
    map_jec_ = cs_jerc_->compound().at("Summer19UL16APV_V7_MC_L1L2L3Res_AK4PFchs");
  }
  else if (year==2016) {
    cs_jerc_ = correction::CorrectionSet::from_file("data/zgamma/2016postVFP_UL/jet_jerc_2016.json");
    map_jes_ = cs_jerc_->at("Summer19UL16_V7_MC_Total_AK4PFchs");
    map_jersf_ = cs_jerc_->at("Summer20UL16_JRV3_MC_ScaleFactor_AK4PFchs");
    map_jermc_ = cs_jerc_->at("Summer20UL16_JRV3_MC_PtResolution_AK4PFchs");
    map_jec_ = cs_jerc_->compound().at("Summer19UL16_V7_MC_L1L2L3Res_AK4PFchs");
  }
  else if (year==2017) {
    cs_jerc_ = correction::CorrectionSet::from_file("data/zgamma/2017_UL/jet_jerc_2017.json");
    map_jes_ = cs_jerc_->at("Summer19UL17_V5_MC_Total_AK4PFchs");
    map_jersf_ = cs_jerc_->at("Summer19UL17_JRV2_MC_ScaleFactor_AK4PFchs");
    map_jermc_ = cs_jerc_->at("Summer19UL17_JRV2_MC_PtResolution_AK4PFchs");
    map_jec_ = cs_jerc_->compound().at("Summer19UL17_V5_MC_L1L2L3Res_AK4PFchs");
  }
  else if (year==2018) {
    cs_jerc_ = correction::CorrectionSet::from_file("data/zgamma/2018_UL/jet_jerc_2018.json");
    map_jes_ = cs_jerc_->at("Summer19UL18_V5_MC_Total_AK4PFchs");
    map_jersf_ = cs_jerc_->at("Summer19UL18_JRV2_MC_ScaleFactor_AK4PFchs");
    map_jermc_ = cs_jerc_->at("Summer19UL18_JRV2_MC_PtResolution_AK4PFchs");
    map_jec_ = cs_jerc_->compound().at("Summer19UL18_V5_MC_L1L2L3Res_AK4PFchs");
  }
  else {
    std::cout << "WARNING: No dedicated JES/JER uncertainties, defaulting to 2018" << std::endl;
    cs_jerc_ = correction::CorrectionSet::from_file("data/zgamma/2018_UL/jet_jerc_2018.json");
    map_jes_ = cs_jerc_->at("Summer19UL18_V5_MC_Total_AK4PFchs");
    map_jersf_ = cs_jerc_->at("Summer19UL18_JRV2_MC_ScaleFactor_AK4PFchs");
    map_jermc_ = cs_jerc_->at("Summer19UL18_JRV2_MC_PtResolution_AK4PFchs");
    map_jec_ = cs_jerc_->compound().at("Summer19UL18_V5_MC_L1L2L3Res_AK4PFchs");
  }
  in_file_jetveto_ = "data/zgamma/2022/jetvetomaps.json";
  cs_jetveto_ = correction::CorrectionSet::from_file(in_file_jetveto_);
}

JetMetProducer::~JetMetProducer(){
}

//Note this function also writes out MET and its systematic uncertainties from propagating JES/JER uncertainties
void JetMetProducer::GetJetUncertainties(nano_tree &nano, pico_tree &pico, 
                                      vector<float> &jer_nm_factor, 
                                      vector<float> &jer_up_factor,
                                      vector<float> &jer_dn_factor,
                                      vector<float> &jes_up_factor,
                                      vector<float> &jes_dn_factor) {

  //do not call this function for data
  //implementation basically follows the following without the 2017 EE fix
  //https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/jetmetUncertainties.py

  float met_x = nano.MET_pt()*cos(nano.MET_phi());
  float met_y = nano.MET_pt()*sin(nano.MET_phi());
  float met_x_nom = met_x;
  float met_y_nom = met_y;
  float met_x_jesup = met_x;
  float met_y_jesup = met_y;
  float met_x_jesdn = met_x;
  float met_y_jesdn = met_y;
  float met_x_jerup = met_x;
  float met_y_jerup = met_y;
  float met_x_jerdn = met_x;
  float met_y_jerdn = met_y;

  //loop over regular jets and jets that didn't make it into slimmedjets (CorrT1METJets)
  for (int jet_type(0); jet_type<2; jet_type++) {

    vector<float> jet_type_pt, jet_type_eta, jet_type_phi;
    vector<float> jet_type_rawfactor, jet_type_muonfactor;
    int jet_type_size(0);
    if (jet_type==0) {
      jet_type_pt = nano.Jet_pt();
      jet_type_eta = nano.Jet_eta();
      jet_type_phi = nano.Jet_phi();
      jet_type_rawfactor = nano.Jet_rawFactor();
      jet_type_muonfactor = nano.Jet_muonSubtrFactor();
      jet_type_size = nano.nJet();
    }
    else {
      jet_type_pt = nano.CorrT1METJet_rawPt();
      jet_type_eta = nano.CorrT1METJet_eta();
      jet_type_phi = nano.CorrT1METJet_phi();
      jet_type_rawfactor.resize(nano.nCorrT1METJet(),0.0);
      jet_type_muonfactor = nano.CorrT1METJet_muonSubtrFactor();
      jet_type_size = nano.nCorrT1METJet();
    }

    for (int ijet(0); ijet<jet_type_size; ++ijet) {

      //Get JECs and correct nominal pt for CorrT1METJets
      float jec = 1.0/(1.0-jet_type_rawfactor[ijet]);
      if (jet_type==1) {
        jec = map_jec_->evaluate({
              nano.CorrT1METJet_area()[ijet],jet_type_eta[ijet],
              jet_type_pt[ijet],nano.fixedGridRhoFastjetAll()});
        jet_type_pt[ijet] = jet_type_pt[ijet]*jec;
      }
      float jet_raw_pt = jet_type_pt[ijet]/jec;
      float jet_raw_pt_nomu = jet_raw_pt*(1.0-jet_type_muonfactor[ijet]);
      float jet_l1l2l3_pt_nomu = jet_raw_pt_nomu*jec;
      if (jet_type == 1 && jet_l1l2l3_pt_nomu < 15) continue;

      //calculate JER (smearing) factors
      //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#Smearing_procedures
      float sigmajer = map_jermc_->evaluate({jet_type_eta[ijet],jet_type_pt[ijet],
                                             nano.fixedGridRhoFastjetAll()});
      float sjer_nom = map_jersf_->evaluate({jet_type_eta[ijet],"nom"});
      float sjer_up = map_jersf_->evaluate({jet_type_eta[ijet],"up"});
      float sjer_dn = map_jersf_->evaluate({jet_type_eta[ijet],"down"});
      float indiv_jer_nm(1.0), indiv_jer_up(1.0), indiv_jer_dn(1.0);

      bool found_genjet = false;
      for (int igen(0); igen<nano.nGenJet(); ++igen) {
        float dr = dR(jet_type_eta[ijet], nano.GenJet_eta()[igen], jet_type_phi[ijet], nano.GenJet_phi()[ijet]);
        float dpt = jet_type_pt[ijet]-nano.GenJet_pt()[igen];
        if (dr < 0.2 && fabs(dpt) < 3.0*sigmajer*jet_type_pt[ijet]) {
          found_genjet = true;
          indiv_jer_nm = (1.0+(sjer_nom-1.0)*dpt/jet_type_pt[ijet]);
          indiv_jer_up = (1.0+(sjer_up-1.0)*dpt/jet_type_pt[ijet]);
          indiv_jer_dn = (1.0+(sjer_dn-1.0)*dpt/jet_type_pt[ijet]);
          break;
        }
      }

      if (!found_genjet) {
        indiv_jer_nm = 1.0+rng_.Gaus(0,sigmajer)*sqrt(std::max(sjer_nom*sjer_nom-1.0,0.0));
        indiv_jer_up = 1.0+rng_.Gaus(0,sigmajer)*sqrt(std::max(sjer_up*sjer_up-1.0,0.0));
        indiv_jer_dn = 1.0+rng_.Gaus(0,sigmajer)*sqrt(std::max(sjer_dn*sjer_dn-1.0,0.0));
      }

      //Following NanoAOD-tools, JES uncertainties are evaluated post-smearing
      float jes_unc = map_jes_->evaluate({jet_type_eta[ijet],jet_type_pt[ijet]*indiv_jer_nm});

      //Save values for regular jets
      if (jet_type==0) {
        jer_nm_factor.push_back(indiv_jer_nm);
        jer_up_factor.push_back(indiv_jer_up);
        jer_dn_factor.push_back(indiv_jer_dn);
        jes_up_factor.push_back(1.0+jes_unc);
        jes_dn_factor.push_back(1.0-jes_unc);
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
      //starting from T1 corrected MET in contrast with NanoAOD-tools
      //i.e. L2L3-L1 already done, just need to worry about variations
      if (jet_l1l2l3_pt_nomu > 15 && fabs(jet_type_eta[ijet])<5.2 && emef < 0.9) {
        met_x_nom -= jet_cosphi*(jet_type_pt[ijet]*indiv_jer_nm-jet_type_pt[ijet]);
        met_y_nom -= jet_sinphi*(jet_type_pt[ijet]*indiv_jer_nm-jet_type_pt[ijet]);
        met_x_jerup -= jet_cosphi*(jet_type_pt[ijet]*indiv_jer_up-jet_type_pt[ijet]);
        met_y_jerup -= jet_sinphi*(jet_type_pt[ijet]*indiv_jer_up-jet_type_pt[ijet]);
        met_x_jerdn -= jet_cosphi*(jet_type_pt[ijet]*indiv_jer_dn-jet_type_pt[ijet]);
        met_y_jerdn -= jet_sinphi*(jet_type_pt[ijet]*indiv_jer_dn-jet_type_pt[ijet]);
        met_x_jesup -= jet_cosphi*(jet_type_pt[ijet]*indiv_jer_nm*(1.0+jes_unc)-jet_type_pt[ijet]);
        met_y_jesup -= jet_sinphi*(jet_type_pt[ijet]*indiv_jer_nm*(1.0+jes_unc)-jet_type_pt[ijet]);
        met_x_jesdn -= jet_cosphi*(jet_type_pt[ijet]*indiv_jer_nm*(1.0-jes_unc)-jet_type_pt[ijet]);
        met_y_jesdn -= jet_sinphi*(jet_type_pt[ijet]*indiv_jer_nm*(1.0-jes_unc)-jet_type_pt[ijet]);
      }
    }
  }

  pico.out_sys_met().resize(4,0.0);
  pico.out_sys_met_phi().resize(4,0.0);
  pico.out_met() = sqrt(met_x_nom*met_x_nom+met_y_nom*met_y_nom);
  pico.out_sys_met()[0] = sqrt(met_x_jerup*met_x_jerup+met_y_jerup*met_y_jerup);
  pico.out_sys_met()[1] = sqrt(met_x_jerdn*met_x_jerdn+met_y_jerdn*met_y_jerdn);
  pico.out_sys_met()[2] = sqrt(met_x_jesup*met_x_jesup+met_y_jesup*met_y_jesup);
  pico.out_sys_met()[3] = sqrt(met_x_jesdn*met_x_jesdn+met_y_jesdn*met_y_jesdn);
  pico.out_met_phi() = atan2(met_y_nom, met_x_nom);
  pico.out_sys_met_phi()[0] = atan2(met_y_jerup, met_x_jerup);
  pico.out_sys_met_phi()[1] = atan2(met_y_jerdn, met_x_jerdn);
  pico.out_sys_met_phi()[2] = atan2(met_y_jesup, met_x_jesup);
  pico.out_sys_met_phi()[3] = atan2(met_y_jesdn, met_x_jesdn);
}

//replaces MET_producer for UL and beyond
//writes other MET variables and non-jet MET uncertainties
void JetMetProducer::WriteMet(nano_tree &nano, pico_tree &pico) {
  pico.out_met_calo()    = nano.CaloMET_pt();
  pico.out_met_tru()     = nano.GenMET_pt();
  pico.out_met_tru_phi() = nano.GenMET_phi();
  pico.out_ht_isr_me()   = nano.LHE_HTIncoming();
  if (isData) {
    pico.out_met() = nano.MET_pt();
    pico.out_met_phi() = nano.MET_phi();
    return;
  }

  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETRun2Corrections
  float met_x = pico.out_met()*cos(pico.out_met_phi());
  float met_y = pico.out_met()*sin(pico.out_met_phi());
  float met_x_unclup = met_x + nano.MET_MetUnclustEnUpDeltaX();
  float met_y_unclup = met_y + nano.MET_MetUnclustEnUpDeltaX();
  float met_x_uncldn = met_x - nano.MET_MetUnclustEnUpDeltaX();
  float met_y_uncldn = met_y - nano.MET_MetUnclustEnUpDeltaX();
  float met_x_leptonphotonup = met_x;
  float met_y_leptonphotonup = met_y;
  float met_x_leptonphotondn = met_x;
  float met_y_leptonphotondn = met_y;
  for (int iel = 0; iel < pico.out_nel(); iel++) {
    if (pico.out_el_sig()[iel]) {
      float unc_factor = 0.006; //EB
      if (fabs(pico.out_el_eta()[iel])>1.5051) unc_factor = 0.015; //EE
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
      if (fabs(pico.out_photon_eta()[iph])>1.5051) unc_factor = 0.015; //EE
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
                                   bool isFastsim, 
                                   bool isSignal,
                                   bool is2022preEE,
                                   vector<HiggsConstructionVariables> &sys_higvars){
  vector<int> sig_jet_nano_idx;
  pico.out_njet() = 0; pico.out_ht() = 0; pico.out_ht5() = 0; 
  pico.out_nbl() = 0; pico.out_nbm() = 0; pico.out_nbt() = 0; 
  pico.out_nbdfl() = 0; pico.out_nbdfm() = 0; pico.out_nbdft() = 0; 
  pico.out_ngenjet() = 0;

  //add smearing to jets and calculate uncertainties
  vector<float> Jet_pt, Jet_mass;
  vector<float> jer_nm_factor, jer_up_factor, jer_dn_factor, jes_up_factor, jes_dn_factor;
  float MET_pt, MET_phi;
  if (isData) {
    WriteMet(nano, pico);
    Jet_pt = nano.Jet_pt();
    Jet_mass = nano.Jet_mass();
    jer_nm_factor.resize(Jet_pt.size(),1.0);
    jer_up_factor.resize(Jet_pt.size(),1.0);
    jer_dn_factor.resize(Jet_pt.size(),1.0);
    jes_up_factor.resize(Jet_pt.size(),1.0);
    jes_dn_factor.resize(Jet_pt.size(),1.0);
  }
  else if (is_preUL && year <= 2018) {
    getJetWithJEC(nano, isFastsim, Jet_pt, Jet_mass);
    getMETWithJEC(nano, year, isFastsim, MET_pt, MET_phi, is_preUL);
    jer_nm_factor.resize(Jet_pt.size(),1.0);
    jer_up_factor.resize(Jet_pt.size(),1.0);
    jer_dn_factor.resize(Jet_pt.size(),1.0);
    jes_up_factor.resize(Jet_pt.size(),1.0);
    jes_dn_factor.resize(Jet_pt.size(),1.0);
    met_producer.WriteMet(nano, pico, isFastsim, isSignal, is_preUL);
  }
  else  {
    GetJetUncertainties(nano, pico, jer_nm_factor, jer_up_factor, jer_dn_factor, 
                        jes_up_factor, jes_dn_factor);
    MET_pt = pico.out_met();
    MET_phi = pico.out_met_phi();
    WriteMet(nano, pico);
    for(int ijet(0); ijet<nano.nJet(); ++ijet) {
      Jet_pt.push_back(nano.Jet_pt()[ijet]*jer_nm_factor[ijet]);
      Jet_mass.push_back(nano.Jet_mass()[ijet]*jer_nm_factor[ijet]);
    }
  }
  vector<int> Jet_jetId;
  getJetId(nano, year, Jet_jetId);
  
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

  if (is_preUL && isSignal) {
    pico.out_sys_njet().resize(4,0);
    pico.out_sys_nbl().resize(4,0);
    pico.out_sys_nbm().resize(4,0);
    pico.out_sys_nbt().resize(4,0);
    pico.out_sys_ht().resize(4,0.0);
    sys_jet_met_dphi.resize(4,vector<float>({}));
    sys_higvars.resize(4, HiggsConstructionVariables());
    pico.out_sys_low_dphi_met().resize(4,false);
  }
  else if (!is_preUL) {
    pico.out_sys_njet().resize(4,0);
    //TODO add variations for deepflavor b-tagging
  }

  // saving jet info on all jets passing pt cut, including endcap
  for(int ijet(0); ijet<nano.nJet(); ++ijet){
    if (verbose) cout<<"Jet "<<ijet<<": pt = "<<setw(10)<<Jet_pt[ijet]
                                    <<" eta = "<<setw(10)<<nano.Jet_eta()[ijet]
                                    <<" phi = "<<setw(10)<<nano.Jet_phi()[ijet]
                                    <<" m = "<<setw(10)<<Jet_mass[ijet]
                                    <<endl;


    // check overlap with signal leptons (or photons)
    bool islep = find(jet_islep_nano_idx.begin(), jet_islep_nano_idx.end(), ijet) != jet_islep_nano_idx.end();
    // check overlap with veto leptons
    bool isvlep = find(jet_isvlep_nano_idx.begin(), jet_isvlep_nano_idx.end(), ijet) != jet_isvlep_nano_idx.end();
    // N.B. photon collection is not filled for Higgsino analysis, so there is no overlap removal!
    bool isphoton = find(jet_isphoton_nano_idx.begin(), jet_isphoton_nano_idx.end(), ijet) != jet_isphoton_nano_idx.end();
    // jetid applied to only full sim and data
    bool pass_jetid = true;
    if (!isFastsim) if (Jet_jetId[ijet] <1) pass_jetid = false;

    bool isvetojet = false;

    float veto = 0; 
    float vetoEE = 0;
    double phicorr;
    if(nano.Jet_phi().at(ijet)>3.1415926){ //a dumb addition because sometimes jet phi is slightly larger than pi
      phicorr = 3.1415926;
    } else if (nano.Jet_phi().at(ijet)<-3.1415926){
      phicorr = -3.1415926;
    } else {
      phicorr = nano.Jet_phi().at(ijet);
    }

    if (fabs(nano.Jet_eta()[ijet])<5.191) {
      if (year==2022 && is2022preEE==true){
        map_jetveto_ = cs_jetveto_->at("Winter22Run3_RunCD_V1");
        veto = map_jetveto_->evaluate({"jetvetomap", nano.Jet_eta().at(ijet),phicorr});
      } else if (year==2022 && is2022preEE==false){
        map_jetveto_ = cs_jetveto_->at("Winter22Run3_RunE_V1");
        veto = map_jetveto_->evaluate({"jetvetomap", nano.Jet_eta().at(ijet),phicorr});
        vetoEE = map_jetveto_->evaluate({"jetvetomap_eep", nano.Jet_eta().at(ijet),phicorr});
      }
    }
    if(veto != 0.0 || vetoEE != 0.0) {
      isvetojet = true;
    }

    //don't include pt in isgood until later in order to save systematics
    bool isgood = !islep && !isphoton && (fabs(nano.Jet_eta()[ijet]) <= max_jet_eta) && pass_jetid && !isvetojet;

    //sys_jetvar convention: [0] JER up, [1] JER down, [2] JEC up, [3] JEC down
    //for now, only save sys_ variables
    if (isSignal && is_preUL) {
      if (nano.Jet_pt_jerUp()[ijet] > min_jet_pt) {
        if (isgood) {
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
        if (fabs(nano.Jet_eta()[ijet]) < 2.4)
          sys_jet_met_dphi.at(0).push_back(DeltaPhi(nano.Jet_phi()[ijet], pico.out_sys_met_phi()[0]));
      }
      if (nano.Jet_pt_jerDown()[ijet] > min_jet_pt) {
        if (isgood) {
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
        if (isgood) {
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
        if (fabs(nano.Jet_eta()[ijet]) < 2.4)
          sys_jet_met_dphi.at(2).push_back(DeltaPhi(nano.Jet_phi()[ijet], pico.out_sys_met_phi()[2]));
      }
      if (nano.Jet_pt_jesTotalDown()[ijet] > min_jet_pt) {
        if (isgood) {
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
        if (fabs(nano.Jet_eta()[ijet]) < 2.4)
          sys_jet_met_dphi.at(3).push_back(DeltaPhi(nano.Jet_phi()[ijet], pico.out_sys_met_phi()[3]));
      }
    }

    if (isvetojet) continue;
    if (Jet_pt[ijet] <= min_jet_pt && 
        (Jet_pt[ijet]*jer_up_factor[ijet]/jer_nm_factor[ijet]) <= min_jet_pt &&
        (Jet_pt[ijet]*jer_dn_factor[ijet]/jer_nm_factor[ijet]) <= min_jet_pt &&
        (Jet_pt[ijet]*jes_up_factor[ijet]) <= min_jet_pt) continue;

    if (!is_preUL && isgood) {
      if ((Jet_pt[ijet]*jer_up_factor[ijet]/jer_nm_factor[ijet]) > min_jet_pt)
        pico.out_sys_njet()[0]++;
      if ((Jet_pt[ijet]*jer_dn_factor[ijet]/jer_nm_factor[ijet]) > min_jet_pt)
        pico.out_sys_njet()[1]++;
      if ((Jet_pt[ijet]*jes_up_factor[ijet]) > min_jet_pt)
        pico.out_sys_njet()[2]++;
      if ((Jet_pt[ijet]*jes_dn_factor[ijet]) > min_jet_pt)
        pico.out_sys_njet()[3]++;
    }

    //after this point isgood means with nominal (smeared) pt
    isgood = isgood && (Jet_pt[ijet] >= min_jet_pt);

    switch(year) {
      case 2016:
      case 2017:
      case 2018:
        pico.out_jet_pt().push_back(Jet_pt[ijet]);
        pico.out_jet_eta().push_back(nano.Jet_eta()[ijet]);
        pico.out_jet_phi().push_back(nano.Jet_phi()[ijet]);
        pico.out_jet_m().push_back(Jet_mass[ijet]);
        pico.out_jet_breg_corr().push_back(nano.Jet_bRegCorr()[ijet]);
        pico.out_jet_breg_res().push_back(nano.Jet_bRegRes()[ijet]);
        pico.out_jet_deepcsv().push_back(nano.Jet_btagDeepB()[ijet]);
        pico.out_jet_deepflav().push_back(nano.Jet_btagDeepFlavB()[ijet]);
        pico.out_jet_ne_emef().push_back(nano.Jet_neEmEF()[ijet]);
        pico.out_jet_qgl().push_back(nano.Jet_qgl()[ijet]);
        pico.out_jet_islep().push_back(islep);
        pico.out_jet_isvlep().push_back(isvlep);
        pico.out_jet_isphoton().push_back(isphoton);
        pico.out_jet_isgood().push_back(isgood);
        pico.out_jet_id().push_back(Jet_jetId[ijet]);
        pico.out_jet_mht_dphi().push_back(DeltaPhi(nano.Jet_phi()[ijet], mht_vec.Phi()));
        pico.out_jet_met_dphi().push_back(DeltaPhi(nano.Jet_phi()[ijet], MET_phi));
        pico.out_jet_puid().push_back(nano.Jet_puId()[ijet]);
        pico.out_jet_puid_disc().push_back(nano.Jet_puIdDisc()[ijet]);
        pico.out_sys_jet_pt_jesup().push_back(Jet_pt[ijet]*jes_up_factor[ijet]);
        pico.out_sys_jet_pt_jesdn().push_back(Jet_pt[ijet]*jes_dn_factor[ijet]);
        pico.out_sys_jet_pt_jerup().push_back(Jet_pt[ijet]*jer_up_factor[ijet]/jer_nm_factor[ijet]);
        pico.out_sys_jet_pt_jerdn().push_back(Jet_pt[ijet]*jer_dn_factor[ijet]/jer_nm_factor[ijet]);
        pico.out_sys_jet_m_jesup().push_back(Jet_mass[ijet]*jes_up_factor[ijet]);
        pico.out_sys_jet_m_jesdn().push_back(Jet_mass[ijet]*jes_dn_factor[ijet]);
        pico.out_sys_jet_m_jerup().push_back(Jet_mass[ijet]*jer_up_factor[ijet]/jer_nm_factor[ijet]);
        pico.out_sys_jet_m_jerdn().push_back(Jet_mass[ijet]*jer_dn_factor[ijet]/jer_nm_factor[ijet]);
        break;
      case 2022:
      case 2023:
        pico.out_jet_pt().push_back(Jet_pt[ijet]);
        pico.out_jet_eta().push_back(nano.Jet_eta()[ijet]);
        pico.out_jet_phi().push_back(nano.Jet_phi()[ijet]);
        pico.out_jet_m().push_back(Jet_mass[ijet]);
        //pico.out_jet_breg_corr().push_back(nano.Jet_bRegCorr()[ijet]);
        //pico.out_jet_breg_res().push_back(nano.Jet_bRegRes()[ijet]);
        pico.out_jet_deepcsv().push_back(nano.Jet_btagDeepB()[ijet]);
        pico.out_jet_deepflav().push_back(nano.Jet_btagDeepFlavB()[ijet]);
        pico.out_jet_ne_emef().push_back(nano.Jet_neEmEF()[ijet]);
        //pico.out_jet_qgl().push_back(nano.Jet_qgl()[ijet]);
        pico.out_jet_islep().push_back(islep);
        pico.out_jet_isvlep().push_back(isvlep);
        pico.out_jet_isphoton().push_back(isphoton);
        pico.out_jet_isgood().push_back(isgood);
        pico.out_jet_id().push_back(Jet_jetId[ijet]);
        pico.out_jet_mht_dphi().push_back(DeltaPhi(nano.Jet_phi()[ijet], mht_vec.Phi()));
        pico.out_jet_met_dphi().push_back(DeltaPhi(nano.Jet_phi()[ijet], MET_phi));
        //pico.out_jet_puid().push_back(nano.Jet_puId()[ijet]);
        //pico.out_jet_puid_disc().push_back(nano.Jet_puIdDisc()[ijet]);
        pico.out_sys_jet_pt_jesup().push_back(Jet_pt[ijet]*jes_up_factor[ijet]);
        pico.out_sys_jet_pt_jesdn().push_back(Jet_pt[ijet]*jes_dn_factor[ijet]);
        pico.out_sys_jet_pt_jerup().push_back(Jet_pt[ijet]*jer_up_factor[ijet]/jer_nm_factor[ijet]);
        pico.out_sys_jet_pt_jerdn().push_back(Jet_pt[ijet]*jer_dn_factor[ijet]/jer_nm_factor[ijet]);
        pico.out_sys_jet_m_jesup().push_back(Jet_mass[ijet]*jes_up_factor[ijet]);
        pico.out_sys_jet_m_jesdn().push_back(Jet_mass[ijet]*jes_dn_factor[ijet]);
        pico.out_sys_jet_m_jerup().push_back(Jet_mass[ijet]*jer_up_factor[ijet]/jer_nm_factor[ijet]);
        pico.out_sys_jet_m_jerdn().push_back(Jet_mass[ijet]*jer_dn_factor[ijet]/jer_nm_factor[ijet]);
        break;
      default:
        std::cout<<"Need code for new year in getZGammaJetBr in jetmet_producer.cpp"<<endl;
        exit(1);
    }

    if (!isData) {
      pico.out_jet_hflavor().push_back(nano.Jet_hadronFlavour()[ijet]);
      pico.out_jet_pflavor().push_back(nano.Jet_partonFlavour()[ijet]);
      pico.out_jet_genjet_idx().push_back(nano.Jet_genJetIdx()[ijet]);
    }
    
    // will be overwritten with the overlapping fat jet index, if such exists, in WriteFatJets
    pico.out_jet_fjet_idx().push_back(-999);

    //the jets for the higgs pair with smallest dm will be set to true in hig_producer
    pico.out_jet_h1d().push_back(false);
    pico.out_jet_h2d().push_back(false);

    if (!islep && !isphoton) pico.out_ht5() += Jet_pt[ijet];

    if (isgood) {
      sig_jet_nano_idx.push_back(ijet);
      pico.out_njet()++;
      pico.out_ht() += Jet_pt[ijet];
      if (nano.Jet_btagDeepB()[ijet] > btag_wpts[0]) pico.out_nbl()++; 
      if (nano.Jet_btagDeepB()[ijet] > btag_wpts[1]) pico.out_nbm()++; 
      if (nano.Jet_btagDeepB()[ijet] > btag_wpts[2]) pico.out_nbt()++;
      if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[0]) pico.out_nbdfl()++; 
      if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[1]) pico.out_nbdfm()++; 
      if (nano.Jet_btagDeepFlavB()[ijet] > btag_df_wpts[2]) pico.out_nbdft()++; 
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

  if (isSignal) {
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
      pico.out_genjet_pflavor().push_back(nano.GenJet_partonFlavour()[ijet]);
      pico.out_genjet_hflavor().push_back(int(nano.GenJet_hadronFlavour()[ijet]));
    }
  }

  if (verbose) cout<<"Done with AK4 jets"<<endl;
  return sig_jet_nano_idx;
}

void JetMetProducer::WriteFatJets(nano_tree &nano, pico_tree &pico){
  pico.out_nfjet() = 0; 
  vector<int> FatJet_btagDDBvL;
  getFatJet_btagDDBvL(nano, nanoaod_version, FatJet_btagDDBvL);

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
    pico.out_fjet_deep_md_hbb_btv().push_back(FatJet_btagDDBvL[ifjet]);
    pico.out_fjet_mva_hbb_btv().push_back(nano.FatJet_btagHbb()[ifjet]);
    // Mass-decorrelated DeepAk8, H->bb vs QCD discriminator, endorsed by JME
    pico.out_fjet_deep_md_hbb_jme().push_back(nano.FatJet_deepTagMD_HbbvsQCD()[ifjet]);
    pico.out_fjet_deep_md_tvsqcd().push_back(nano.FatJet_deepTagMD_TvsQCD()[ifjet]);
    pico.out_fjet_deep_tvsqcd().push_back(nano.FatJet_deepTag_TvsQCD()[ifjet]);

    pico.out_fjet_subjet_idx1().push_back(nano.FatJet_subJetIdx1()[ifjet]);
    pico.out_fjet_subjet_idx2().push_back(nano.FatJet_subJetIdx2()[ifjet]);

    for(unsigned ijet(0); ijet<pico.out_jet_pt().size(); ++ijet){
      if (dR(pico.out_jet_eta()[ijet], nano.FatJet_eta()[ifjet], pico.out_jet_phi()[ijet], nano.FatJet_phi()[ifjet])<0.8)
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

    pico.out_subfjet_deepcsv().push_back(nano.SubJet_btagDeepB()[isubj]);
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
        if (mindr<0.4 && matched_ak4_jets.find(ijet)==matched_ak4_jets.end()) {
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
                                   vector<int> &sig_jet_nano_idx, const float &btag_wpt, bool isFastsim) {

  vector<float> Jet_pt, Jet_mass;
  getJetWithJEC(nano, isFastsim, Jet_pt, Jet_mass);

  TLorentzVector jetsys_v4, jetsys_nob_v4;
  int njet_nob(0);
  for (auto &idx: sig_jet_nano_idx) {
    TLorentzVector ijet_v4;
    ijet_v4.SetPtEtaPhiM(Jet_pt[idx], nano.Jet_eta()[idx], nano.Jet_phi()[idx], Jet_mass[idx]);
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


