#include "photon_producer.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TRandom3.h"

#include "correction.hpp"
#include "utilities.hpp"

using namespace std;

PhotonProducer::PhotonProducer(string year_, bool isData_, 
                               float nanoaod_version_){
  year = year_;
  isData = isData_;
  nanoaod_version = nanoaod_version_;
  if (year=="2016APV") {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2016preVFP_UL/EGM_ScaleUnc.json");
    str_scale_syst_ = "2016preVFP";
    map_scale_syst_ = cs_scale_syst_->at("UL-EGM_ScaleUnc");
  }
  else if (year=="2016") {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2016postVFP_UL/EGM_ScaleUnc.json");
    str_scale_syst_ = "2016postVFP";
    map_scale_syst_ = cs_scale_syst_->at("UL-EGM_ScaleUnc");
  }
  else if (year=="2017") {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2017_UL/EGM_ScaleUnc.json");
    str_scale_syst_ = "2017";
    map_scale_syst_ = cs_scale_syst_->at("UL-EGM_ScaleUnc");
  }
  else if (year=="2018") {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2018_UL/EGM_ScaleUnc.json");
    str_scale_syst_ = "2018";
    map_scale_syst_ = cs_scale_syst_->at("UL-EGM_ScaleUnc");
  }
  else if (year=="2022") {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2022/photonSS_EtDependent.json");
    map_scale_ = cs_scale_syst_->compound().at(
        "EGMScale_Compound_Pho_2022preEE");
    map_smearing_ = cs_scale_syst_->at(
        "EGMSmearAndSyst_PhoPTsplit_2022preEE");
  }
  else if (year=="2022EE") {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2022EE/photonSS_EtDependent.json");
    map_scale_ = cs_scale_syst_->compound().at(
        "EGMScale_Compound_Pho_2022postEE");
    map_smearing_ = cs_scale_syst_->at(
        "EGMSmearAndSyst_PhoPTsplit_2022postEE");
  }
  else if (year=="2023") {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2023/photonSS_EtDependent.json");
    map_scale_ = cs_scale_syst_->compound().at(
        "EGMScale_Compound_Pho_2023preBPIX");
    map_smearing_ = cs_scale_syst_->at(
        "EGMSmearAndSyst_PhoPTsplit_2023preBPIX");
  }
  else if (year=="2023BPix") {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2023BPix/photonSS_EtDependent.json");
    map_scale_ = cs_scale_syst_->compound().at(
        "EGMScale_Compound_Pho_2023postBPIX");
    map_smearing_ = cs_scale_syst_->at(
        "EGMSmearAndSyst_PhoPTsplit_2023postBPIX");
  }
  else {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2023BPix/photonSS_EtDependent.json");
    map_scale_ = cs_scale_syst_->compound().at(
        "EGMScale_Compound_Pho_2023postBPIX");
    map_smearing_ = cs_scale_syst_->at(
        "EGMSmearAndSyst_PhoPTsplit_2023postBPIX");
    std::cout << "WARNING: No dedicated EGM scale/smearing JSONs, defaulting to 2023BPix" << std::endl;
  }
}

PhotonProducer::~PhotonProducer(){
}

bool PhotonProducer::IsSignal(nano_tree &nano, pico_tree &pico, int nano_idx, 
                              float scaleres_corr, float minLepDR, 
                              vector<int> &photon_el_pico_idx) {
  float pt = nano.Photon_pt()[nano_idx]*scaleres_corr;
  if (pt < SignalPhotonPtCut) return false;
  if (!(nano.Photon_isScEtaEB()[nano_idx] 
        || nano.Photon_isScEtaEE()[nano_idx])) return false;
  if (!nano.Photon_mvaID_WP80()[nano_idx]) return false;
  if (!nano.Photon_electronVeto()[nano_idx]) return false;
  if (!(minLepDR > 0.3f)) return false;
  if (!(photon_el_pico_idx[nano_idx]==-1 || !(pico.out_el_sig()[photon_el_pico_idx[nano_idx]]))) return false;
  return true;
}

vector<int> PhotonProducer::WritePhotons(nano_tree &nano, pico_tree &pico, vector<int> &jet_isphoton_nano_idx, vector<int> &sig_el_nano_idx, vector<int> &sig_mu_nano_idx, vector<int> &photon_el_pico_idx){
  pico.out_nphoton() = 0;
  // pico.out_nfsrphoton() = 0;
  vector<int> sig_photon_nano_idx;
  int nphotons(0);

  vector<int> FsrPhoton_muonIdx;
  getFsrPhoton_muonIdx(nano, nanoaod_version, FsrPhoton_muonIdx);
  vector<int> Photon_jetIdx;
  getPhoton_jetIdx(nano, nanoaod_version, Photon_jetIdx);
  vector<int> Photon_cutBased;
  getPhoton_cutBased(nano, nanoaod_version, Photon_cutBased);

  int hig019014_photon_idx = -1;
  //Set a default value 
  pico.out_photon_idx_hig019014() = -1;

  vector<int> GenPart_statusFlags;
  if (!isData)
    getGenPart_statusFlags(nano, nanoaod_version, GenPart_statusFlags);

  //calculate scale/resolution corrections and drmin/drmax
  vector<float> scaleres_corr;
  vector<float> scale_syst_up;
  vector<float> scale_syst_dn;
  vector<float> smear_syst_up;
  vector<float> smear_syst_dn;
  vector<float> photon_drmin;
  vector<float> photon_drmax;
  for(int iph(0); iph<nano.nPhoton(); ++iph){
    //deal with scale/smearing (systematics only for NanoAODv9 [run 2], full
    //correction for NanoAODv10+ [run3])
    float pt = nano.Photon_pt()[iph];
    float eta = nano.Photon_eta()[iph];
    float phi = nano.Photon_phi()[iph];
    if (year=="2016APV"||year=="2016"||year=="2017"||year=="2018") {
      scaleres_corr.push_back(1.0f);
      if (!isData) {
        scale_syst_up.push_back(map_scale_syst_->evaluate({str_scale_syst_,
            "scaleup",eta,nano.Photon_seedGain()[iph]}));
        scale_syst_dn.push_back(map_scale_syst_->evaluate({str_scale_syst_,
            "scaledown",eta,nano.Photon_seedGain()[iph]}));
        smear_syst_up.push_back(1.0f+nano.Photon_dEsigmaUp()[iph]);
        smear_syst_dn.push_back(1.0f+nano.Photon_dEsigmaDown()[iph]);
      }
    }
    else if ((year=="2022"||year=="2022EE"||year=="2023"||year=="2023BPix") 
             && pt>20) {
      float run = static_cast<float>(nano.run());
      float r9 = fmin(fmax(nano.Photon_r9()[iph],0.0),1.0);
      float seedGain = static_cast<float>(nano.Photon_seedGain()[iph]);
      if (isData) {
        //scale corrections applied to data
        scaleres_corr.push_back(map_scale_->evaluate({"scale",run,eta,r9,
            fabs(eta),pt,seedGain}));
      }
      else {
        //smearing corrections applied to MC, syst.s also calculated
        float rho = map_smearing_->evaluate({"smear",pt,r9,fabs(eta)});
        float err_rho = map_smearing_->evaluate({"esmear",pt,r9,fabs(eta)});
        float scale_unc = map_smearing_->evaluate({"escale",pt,r9,fabs(eta)});
        float rand = rng_.Gaus();
        scaleres_corr.push_back(1.0f+rand*rho);
        smear_syst_up.push_back(1.0f+rand*(rho+err_rho));
        smear_syst_dn.push_back(1.0f+rand*(rho-err_rho));
        scale_syst_up.push_back(1.0f+scale_unc);
        scale_syst_dn.push_back(1.0f-scale_unc);
      }
    }
    else {
      scaleres_corr.push_back(1.0f);
      smear_syst_up.push_back(1.0f);
      smear_syst_dn.push_back(1.0f);
      scale_syst_up.push_back(1.0f);
      scale_syst_dn.push_back(1.0f);
    }

    // Find min(dR) between photon and signal lepton
    double minLepDR(999.);
    double maxLepDR(0.);
    for(size_t iel(0); iel<sig_el_nano_idx.size(); iel++) {
      double tempDR = dR(eta, nano.Electron_eta()[sig_el_nano_idx.at(iel)],
                         phi, nano.Electron_phi()[sig_el_nano_idx.at(iel)]);
      if(tempDR < minLepDR) minLepDR = tempDR;
      if(tempDR > maxLepDR) maxLepDR = tempDR;
    }
    for(size_t imu(0); imu<sig_mu_nano_idx.size(); imu++) {
      double tempDR = dR(eta, nano.Muon_eta()[sig_mu_nano_idx.at(imu)],
                         phi, nano.Muon_phi()[sig_mu_nano_idx.at(imu)]);
      if(tempDR < minLepDR) minLepDR = tempDR;
      if(tempDR > maxLepDR) maxLepDR = tempDR;
    }
    photon_drmin.push_back(minLepDR);
    photon_drmax.push_back(maxLepDR);
  }

  //first, determine ordering based on signal and pt
  std::vector<NanoOrderEntry> nano_entries;
  for(int iph(0); iph<nano.nPhoton(); ++iph){
    NanoOrderEntry nano_entry;
    nano_entry.nano_idx = iph;
    nano_entry.pt = nano.Photon_pt()[iph]*scaleres_corr[iph];
    nano_entry.is_sig = IsSignal(nano, pico, iph, scaleres_corr[iph], 
                                 photon_drmin[iph], photon_el_pico_idx);
    nano_entries.push_back(nano_entry);
  }
  std::sort(nano_entries.begin(),nano_entries.end(), 
      [](NanoOrderEntry a, NanoOrderEntry b) {
        if (a.is_sig && !b.is_sig) return true;
        if (b.is_sig && !a.is_sig) return false;
        return (a.pt>b.pt);
      });
  std::vector<int> ordered_nano_indices;
  for (NanoOrderEntry nano_entry : nano_entries)
    ordered_nano_indices.push_back(nano_entry.nano_idx);

  //then, add branches to pico
  int pico_idx = 0;
  for(int iph : ordered_nano_indices) {
    float raw_pt = nano.Photon_pt()[iph];
    float pt = raw_pt*scaleres_corr[iph];
    float eta = nano.Photon_eta()[iph];
    float phi = nano.Photon_phi()[iph];
    float mva = nano.Photon_mvaID()[iph];
    bool eVeto = nano.Photon_electronVeto()[iph];

    if (pt <= PhotonPtCut) continue;
    if (!(nano.Photon_isScEtaEB()[iph] || nano.Photon_isScEtaEE()[iph])) 
      continue;

    bool isSignal = IsSignal(nano, pico, iph, scaleres_corr[iph], 
                             photon_drmin[iph], photon_el_pico_idx);

    //some logic to save highest pt photon passing HIG-19-014 selection
    bool isSignalhig019014 = (((nano.Photon_isScEtaEB()[iph] && mva > -0.4f) 
        || (nano.Photon_isScEtaEE()[iph] && mva > -0.58f)) && eVeto 
        && photon_drmin[iph] > 0.4f && pt > SignalPhotonPtCut);
    if (isSignalhig019014) {
      if (hig019014_photon_idx==-1 
          || pt > pico.out_photon_pt()[hig019014_photon_idx]) {
        pico.out_photon_idx_hig019014() = pico_idx;
        hig019014_photon_idx = pico_idx;
      }
    }

    //calculate eta w.r.t. origin (SCeta)
    float origin_eta = 0.0;
    if (nano.Photon_isScEtaEB()[iph]) {
      float pv_tan_theta_over_2 = exp(-1.0*eta);
      float pv_tan_theta = 2.0*pv_tan_theta_over_2/(1.0-pv_tan_theta_over_2*pv_tan_theta_over_2);
      float photon_unit_x = cos(nano.Photon_phi()[iph]);
      float photon_unit_y = sin(nano.Photon_phi()[iph]);
      float pv_ecal_dr = 130.0 - (nano.PV_x()*photon_unit_x+nano.PV_y()*photon_unit_y);
      float pv_ecal_dz = pv_ecal_dr/pv_tan_theta;
      float origin_theta = atan(130.0/(nano.PV_z()+pv_ecal_dz));
      if (origin_theta < 0) origin_theta += M_PI;
      origin_eta = -1.0*log(tan(origin_theta/2.0));
    }
    else { //if (nano.Photon_isScEtaEE()[iph])
      float pv_tan_theta_over_2 = exp(-1.0*eta);
      float pv_tan_theta = 2.0*pv_tan_theta_over_2/(1.0-pv_tan_theta_over_2*pv_tan_theta_over_2);
      float photon_unit_x = cos(nano.Photon_phi()[iph]);
      float photon_unit_y = sin(nano.Photon_phi()[iph]);
      float pv_ecal_dz = 310.0-nano.PV_z(); //+ endcap
      if (eta < 0) pv_ecal_dz = 310.0+nano.PV_z(); //- endcap
      float pv_ecal_dr = pv_ecal_dz*pv_tan_theta;
      float origin_theta = atan(((photon_unit_x*nano.PV_x()+photon_unit_y*nano.PV_y())+pv_ecal_dr)/310.0);
      if (origin_theta < 0) origin_theta += M_PI;
      origin_eta = -1.0*log(tan(origin_theta/2.0));
    }

    //find nearest jet and get PUID
    float photon_puid_disc = -1.0;
    float min_jet_dr = 999.0;
    if (nanoaod_version < 9.99) { //PUID not available for run 3
      for (int ijet = 0; ijet < nano.nJet(); ijet++) {
        float ph_jet_dr = dR(eta, nano.Jet_eta()[ijet], phi, nano.Jet_phi()[ijet]);
        if (ph_jet_dr < 0.4) {
          if (ph_jet_dr < min_jet_dr) {
            min_jet_dr = ph_jet_dr;
            photon_puid_disc = nano.Jet_puIdDisc()[ijet];
          }
        }
      }
    }

    //check for nearby gen particles with pT>15 GeV and fromHardprocess
    bool hardprocess = false;
    if (!isData) {
      for (int imc = 0; imc < nano.nGenPart(); imc++) {
        if (nano.GenPart_pt()[imc] > 15
            && ((GenPart_statusFlags[imc] & 0x100) != 0)) {
          //skip promptly decaying particles (W Z t H)
          int abs_mc_id = abs(nano.GenPart_pdgId()[imc]);
          if (abs_mc_id == 6 || abs_mc_id == 23 || abs_mc_id == 24
              || abs_mc_id ==25) continue;
          //use large radius to capture ex. fragmentation photons in jets
          if (dR(eta,nano.GenPart_eta()[imc],phi,nano.GenPart_phi()[imc])<0.4) {
            hardprocess = true;
          }
        }
      }
    }
    
    pico.out_photon_pt().push_back(pt);
    pico.out_photon_pt_raw().push_back(raw_pt);
    pico.out_photon_eta().push_back(eta);
    pico.out_photon_phi().push_back(phi);
    pico.out_photon_r9().push_back(nano.Photon_r9()[iph]);
    pico.out_photon_sieie().push_back(nano.Photon_sieie()[iph]);
    pico.out_photon_hoe().push_back(nano.Photon_hoe()[iph]);
    pico.out_photon_energyErr().push_back(nano.Photon_energyErr()[iph]);
    pico.out_photon_elveto().push_back(eVeto);
    pico.out_photon_isScEtaEB().push_back(nano.Photon_isScEtaEB()[iph]);
    pico.out_photon_isScEtaEE().push_back(nano.Photon_isScEtaEE()[iph]);
    pico.out_photon_sig().push_back(isSignal);
    pico.out_photon_drmin().push_back(photon_drmin[iph]);
    pico.out_photon_drmax().push_back(photon_drmax[iph]);
    pico.out_photon_elidx().push_back(photon_el_pico_idx[iph]);
    pico.out_photon_origin_eta().push_back(origin_eta);
    pico.out_photon_pixelseed().push_back(nano.Photon_pixelSeed()[iph]);
    pico.out_photon_idmva().push_back(mva);
    pico.out_photon_idCutBasedBitMap().push_back(nano.Photon_vidNestedWPBitmap()[iph]);
    pico.out_photon_idCutBased().push_back(Photon_cutBased[iph]);
    pico.out_photon_id().push_back(nano.Photon_mvaID_WP90()[iph]);
    pico.out_photon_id80().push_back(nano.Photon_mvaID_WP80()[iph]);
    if (!isData) {
      pico.out_photon_pflavor().push_back(nano.Photon_genPartFlav()[iph]);
      pico.out_photon_hardprocess().push_back(hardprocess);
      pico.out_sys_photon_pt_resup().push_back(raw_pt*smear_syst_up[iph]);
      pico.out_sys_photon_pt_resdn().push_back(raw_pt*smear_syst_dn[iph]);
      pico.out_sys_photon_pt_scaleup().push_back(pt*scale_syst_up[iph]);
      pico.out_sys_photon_pt_scaledn().push_back(pt*scale_syst_dn[iph]);
    }
    if (year=="2016APV"||year=="2016"||year=="2017"||year=="2018") {
      pico.out_photon_ecorr().push_back(nano.Photon_eCorr()[iph]);
      pico.out_photon_reliso().push_back(nano.Photon_pfRelIso03_all()[iph]);
      pico.out_photon_pudisc().push_back(photon_puid_disc);
    }
    else if (year=="2022"||year=="2022EE"||year=="2023"||year=="2023BPix") {
      pico.out_photon_reliso().push_back(nano.Photon_pfRelIso03_all_quadratic()[iph]);
      pico.out_photon_phiso().push_back(nano.Photon_pfPhoIso03()[iph]);
      pico.out_photon_chiso().push_back(nano.Photon_pfChargedIsoPFPV()[iph]);
      pico.out_photon_chiso_worst().push_back(nano.Photon_pfChargedIsoWorstVtx()[iph]);
      pico.out_photon_s4().push_back(nano.Photon_s4()[iph]);
      pico.out_photon_sieip().push_back(nano.Photon_sieip()[iph]);
      pico.out_photon_etawidth().push_back(nano.Photon_etaWidth()[iph]);
      pico.out_photon_phiwidth().push_back(nano.Photon_phiWidth()[iph]);
      pico.out_photon_energy_raw().push_back(nano.Photon_energyRaw()[iph]);
      pico.out_photon_esoversc().push_back(nano.Photon_esEnergyOverRawE()[iph]);
      pico.out_photon_essigmarr().push_back(nano.Photon_esEffSigmaRR()[iph]);
    }
    else {
      std::cout<<"Need code for new year in getZGammaPhBr (in photon_producer.cpp)"<<endl;
      exit(1);
    }
    if (nanoaod_version > 9.49 && nanoaod_version < 9.51) { //custom NanoAOD production
      pico.out_photon_phiso().push_back(nano.Photon_pfPhoIso03()[iph]);
      pico.out_photon_chiso().push_back(nano.Photon_pfChargedIso()[iph]); //note different name from NanoAODv12
      pico.out_photon_chiso_worst().push_back(nano.Photon_pfChargedIsoWorstVtx()[iph]);
      pico.out_photon_s4().push_back(nano.Photon_s4()[iph]);
      pico.out_photon_sieip().push_back(nano.Photon_sieip()[iph]);
      pico.out_photon_etawidth().push_back(nano.Photon_etaWidth()[iph]);
      pico.out_photon_phiwidth().push_back(nano.Photon_phiWidth()[iph]);
      pico.out_photon_energy_raw().push_back(nano.Photon_energyRaw()[iph]);
      pico.out_photon_esoversc().push_back(nano.Photon_esEnergyOverRawE()[iph]);
      pico.out_photon_essigmarr().push_back(nano.Photon_esEffSigmaRR()[iph]);
    }

    nphotons++;

    // All photons with pt > 15 GeV are considered for creating the ZGamma candidates
    if(isSignal) {
      pico.out_nphoton()++;
      sig_photon_nano_idx.push_back(iph);
      // save indices of matching jets
      if (Photon_jetIdx[iph]>=0)
        jet_isphoton_nano_idx.push_back(Photon_jetIdx[iph]);
      for (int ijet(0); ijet<nano.nJet(); ijet++)
        if (dR(eta, nano.Jet_eta()[ijet], phi, nano.Jet_phi()[ijet])<0.4f)
          jet_isphoton_nano_idx.push_back(ijet);
    }
    pico_idx++;
  }

  //Saves indices for FSR photons
  pico.out_nfsrphoton() = 0;
  for(int iph(0); iph < nano.nFsrPhoton(); ++iph){
    //These set of selections require the photon be close to the corresponding muon and relatively soft
    if (nano.FsrPhoton_pt()[iph] <= FsrPhotonPtCut) continue;
    if (fabs(nano.FsrPhoton_eta()[iph]) > FsrPhotonEtaCut) continue;
    if (nano.FsrPhoton_relIso03()[iph] > FsrPhotonIsoCut) continue;
    if (nano.FsrPhoton_dROverEt2()[iph] > FsrPhotondRCut) continue;

    //Check for separation between fsrphoton and first photon candidate
    if ((pico.out_nphoton() > 0) && dR(pico.out_photon_eta()[0],nano.FsrPhoton_eta()[iph],pico.out_photon_phi()[0],nano.FsrPhoton_phi()[iph]) < FsrSeparationReq) continue;
    
    //Add the values to the pico trees
    pico.out_fsrphoton_pt().push_back(nano.FsrPhoton_pt()[iph]);
    pico.out_fsrphoton_eta().push_back(nano.FsrPhoton_eta()[iph]);
    pico.out_fsrphoton_phi().push_back(nano.FsrPhoton_phi()[iph]);
    pico.out_fsrphoton_reliso().push_back(nano.FsrPhoton_relIso03()[iph]);
    pico.out_fsrphoton_muonidx().push_back(FsrPhoton_muonIdx[iph]);
    pico.out_fsrphoton_droveret2().push_back(nano.FsrPhoton_dROverEt2()[iph]);
    if (nanoaod_version +0.01 > 12) pico.out_fsrphoton_electronidx().push_back(nano.FsrPhoton_electronIdx()[iph]);
    pico.out_nfsrphoton()++;
  } 

  return sig_photon_nano_idx;
}

bool PhotonProducer::idPhoton(int bitmap, int level){
  // decision for each cut represented by 1 bit
  //0 - MinPtCut
  //1 - PhoSCEtaMultiRangeCut
  //2 - PhoSingleTowerHadOverEmCut
  //3 - PhoFull5x5SigmaIEtaIEtaCut
  //4 - PhoGenericRhoPtScaledCut
  //5 - PhoGenericRhoPtScaledCut
  //6 - PhoGenericRhoPtScaledCut
  bool pass = true;
  for (int i(0); i<7; i++){
    if (((bitmap >> i*2) & 0x3) < level) pass = false;
  }
  return pass;
}

