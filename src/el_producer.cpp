#include "el_producer.hpp"

#include "correction.hpp"
#include "utilities.hpp"

#include "TRandom3.h"

#include <algorithm>
#include <iomanip>
#include <memory>
#include <string>

using namespace std;

ElectronProducer::ElectronProducer(string year_, bool isData_, float nanoaod_version_){
  year = year_;
  isData = isData_;
  rng_ = TRandom3(4357);
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
        "data/zgamma/2022/electronSS.json");
    map_scale_ = cs_scale_syst_->at("2022Re-recoBCD_ScaleJSON");
    map_smearing_ = cs_scale_syst_->at("2022Re-recoBCD_SmearingJSON");
  }
  else if (year=="2022EE") {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2022EE/electronSS.json");
    map_scale_ = cs_scale_syst_->at("2022Re-recoE+PromptFG_ScaleJSON");
    map_smearing_ = cs_scale_syst_->at("2022Re-recoE+PromptFG_SmearingJSON");
  }
  else {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2022EE/electronSS.json");
    map_scale_ = cs_scale_syst_->at("2022Re-recoE+PromptFG_ScaleJSON");
    map_smearing_ = cs_scale_syst_->at("2022Re-recoE+PromptFG_SmearingJSON");
    std::cout << "WARNING: No dedicated EGM scale/smearing JSONs, defaulting to 2022EE" << std::endl;
  }
  nanoaod_version = nanoaod_version_;
}

ElectronProducer::~ElectronProducer(){
}

bool ElectronProducer::IsSignal(nano_tree &nano, int nano_idx, bool isZgamma, float scaleres_corr) {
  float pt = nano.Electron_pt()[nano_idx]*scaleres_corr;
  float eta = nano.Electron_eta()[nano_idx];
  float etasc = nano.Electron_deltaEtaSC()[nano_idx] + nano.Electron_eta()[nano_idx];
  float dz = nano.Electron_dz()[nano_idx];
  float dxy = nano.Electron_dxy()[nano_idx];
  float miniiso = nano.Electron_miniPFRelIso_all()[nano_idx];
  if (isZgamma) {
    if (pt <= ZgElectronPtCut) return false;
    if (fabs(etasc) > ElectronEtaCut) return false;
    if (fabs(dz) > dzCut) return false;
    if (fabs(dxy) > dxyCut) return false; 
    if (year=="2016APV"||year=="2016"||year=="2017"||year=="2018") {
      return nano.Electron_mvaFall17V2Iso_WPL()[nano_idx];
    }
    else if (year=="2022"||year=="2022EE"||year=="2023"||year=="2023BPix") {
       return HzzId_WP2022(pt, etasc, nano.Electron_mvaHZZIso()[nano_idx]);
    }
    else {
      std::cout << "WARNING: year value not found in cases. Defaulting Electron MVA to WP90" << std::endl;
      return nano.Electron_mvaIso_WP90()[nano_idx];
    }
  }
  else {
    pt = nano.Electron_pt()[nano_idx]/nano.Electron_eCorr()[nano_idx];
    int bitmap = nano.Electron_vidNestedWPBitmap()[nano_idx];
    bool isBarrel = fabs(eta) <= 1.479;
    bool id = idElectron_noIso(bitmap,3);
    if (pt <= SignalElectronPtCut) return false;
    if (fabs(eta) > ElectronEtaCut) return false;
    if (!idElectron_noIso(bitmap, 1)) return false;
    if ((isBarrel && fabs(dz)>0.10f) || (!isBarrel && fabs(dz)>0.20f)) return false;
    if ((isBarrel && fabs(dxy)>0.05f) || (!isBarrel && fabs(dxy)>0.10f)) return false; 
    return (id && miniiso < ElectronMiniIsoCut);
  }
  return false;
}

vector<int> ElectronProducer::WriteElectrons(nano_tree &nano, pico_tree &pico, vector<int> &jet_islep_nano_idx, vector<int> &jet_isvlep_nano_idx, vector<int> &sig_el_pico_idx, vector<int> &photon_el_pico_idx, bool isZgamma, bool isFastsim){
  vector<float> Jet_pt, Jet_mass;
  getJetWithJEC(nano, isFastsim, Jet_pt, Jet_mass);
  vector<int> Electron_photonIdx;
  getElectron_photonIdx(nano, nanoaod_version, Electron_photonIdx);
  vector<int> Photon_electronIdx;
  getPhoton_electronIdx(nano, nanoaod_version, Photon_electronIdx);

  for (int iph(0); iph<nano.nPhoton(); ++iph)
    photon_el_pico_idx.push_back(-1);

  //first, determine ordering based on signal and pt
  std::vector<NanoOrderEntry> nano_entries;
  for(int iel(0); iel<nano.nElectron(); ++iel){
    NanoOrderEntry nano_entry;
    nano_entry.nano_idx = iel;
    nano_entry.pt = nano.Electron_pt()[iel];
    nano_entry.is_sig = IsSignal(nano, iel, isZgamma);
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
  vector<int> sig_el_nano_idx;
  pico.out_nel() = 0; pico.out_nvel() = 0;
  int pico_idx = 0;
  for(int iel : ordered_nano_indices) {
    float pt = nano.Electron_pt()[iel];
    float eta = nano.Electron_eta()[iel];
    float etasc = nano.Electron_deltaEtaSC()[iel] + nano.Electron_eta()[iel];
    float phi = nano.Electron_phi()[iel];
    float dz = nano.Electron_dz()[iel];
    float dxy = nano.Electron_dxy()[iel];
    float miniiso = nano.Electron_miniPFRelIso_all()[iel];
    bool isSignal = false;
    bool id = false;
    float scaleres_corr = 1.0f;
    if(isZgamma) { // For Zgamma productions
      if (fabs(etasc) > ElectronEtaCut) continue;
      if (fabs(dz) > dzCut)  continue;
      if (fabs(dxy) > dxyCut) continue; 
      float scale_syst_up = 1.0f;
      float scale_syst_dn = 1.0f;
      float smearing_syst_up = 1.0f;
      float smearing_syst_dn = 1.0f;
      //deal with scale/smearing (systematics only for NanoAODv9 [run 2], full
      //correction for NanoAODv10+ [run3])
      if ((year=="2016APV"||year=="2016"||year=="2017"||year=="2018") 
          && !isData) {
        scale_syst_up = map_scale_syst_->evaluate({str_scale_syst_,"scaleup",
            etasc,nano.Electron_seedGain()[iel]});
        scale_syst_dn = map_scale_syst_->evaluate({str_scale_syst_,"scaledown",
            etasc,nano.Electron_seedGain()[iel]});
        smearing_syst_up = 1.0f+nano.Electron_dEsigmaUp()[iel];
        smearing_syst_dn = 1.0f+nano.Electron_dEsigmaDown()[iel];
      }
      else if (year=="2022"||year=="2022EE") {
        //TODO add 2023(BPix) when available
        float run = static_cast<float>(nano.run());
        if (isData) {
          //scale corrections applied to data
          scaleres_corr = map_scale_->evaluate({"total_correction",
              nano.Electron_seedGain()[iel],run,etasc,
              nano.Electron_r9()[iel],pt});
        }
        else {
          //smearing corrections applied to MC, syst.s also calculated
          float rho = map_smearing_->evaluate({"rho",etasc,
              nano.Electron_r9()[iel]});
          float err_rho = map_smearing_->evaluate({"err_rho",etasc,
              nano.Electron_r9()[iel]});
          float scale_unc = map_scale_->evaluate({"total_uncertainty",
              nano.Electron_seedGain()[iel],run,etasc,
              nano.Electron_r9()[iel],pt});
          float rand = rng_.Gaus();
          scaleres_corr = 1.0f+rand*rho;
          smearing_syst_up = 1.0f+rand*(rho+err_rho);
          smearing_syst_dn = 1.0f+rand*(rho-err_rho);
          scale_syst_up = (1.0f+scale_unc);
          scale_syst_dn = (1.0f-scale_unc);
        }
      }
      if (scaleres_corr*pt <= ZgElectronPtCut) continue;
      isSignal = IsSignal(nano, iel, isZgamma, scaleres_corr);

      pico.out_el_sip3d().push_back(nano.Electron_sip3d()[iel]);
      pico.out_el_phidx().push_back(Electron_photonIdx[iel]);
      pico.out_el_etPt().push_back(nano.Electron_scEtOverPt()[iel]);
      pico.out_el_eminusp().push_back(nano.Electron_eInvMinusPInv()[iel]);
      if (!isData) {
        pico.out_sys_el_pt_resup().push_back(pt*smearing_syst_up);
        pico.out_sys_el_pt_resdn().push_back(pt*smearing_syst_dn);
        pico.out_sys_el_pt_scaleup().push_back(pt*scaleres_corr*scale_syst_up);
        pico.out_sys_el_pt_scaledn().push_back(pt*scaleres_corr*scale_syst_dn);
      }
      if (year=="2016APV"||year=="2016"||year=="2017"||year=="2018") {
        pico.out_el_idmva().push_back(nano.Electron_mvaFall17V2Iso()[iel]);
        pico.out_el_id80().push_back(nano.Electron_mvaFall17V2Iso_WP80()[iel]);
        pico.out_el_id90().push_back(nano.Electron_mvaFall17V2Iso_WP90()[iel]);
        pico.out_el_idLoose().push_back(nano.Electron_mvaFall17V2Iso_WPL()[iel]);
      }
      else if (year=="2022"||year=="2022EE"||year=="2023"||year=="2023BPix") {
        bool hzz_wp2022 = HzzId_WP2022(scaleres_corr*pt,etasc,
                                       nano.Electron_mvaHZZIso()[iel]);
        pico.out_el_idmva().push_back(nano.Electron_mvaIso()[iel]);
        pico.out_el_idmvaHZZ().push_back(nano.Electron_mvaHZZIso()[iel]);
        pico.out_el_id80().push_back(nano.Electron_mvaIso_WP80()[iel]);
        pico.out_el_id90().push_back(nano.Electron_mvaIso_WP90()[iel]);
        pico.out_el_idLoose().push_back(hzz_wp2022);
        pico.out_el_fsrphotonidx().push_back(nano.Electron_fsrPhotonIdx()[iel]);
      }
      else {
        std::cout<<"Need code for new year in getZGammaElBr (in el_producer.cpp)"<<endl;
        exit(1);
      }
      int bitmap = nano.Electron_vidNestedWPBitmapHEEP()[iel];
      pico.out_el_ecal().push_back(EcalDriven(bitmap));
    }
    else {
      // Redefine pt and eta to match RA2B ntuples
      pt = nano.Electron_pt()[iel]/nano.Electron_eCorr()[iel];
      if (pt <= VetoElectronPtCut) continue;
      if (fabs(eta) > ElectronEtaCut) continue;
      int bitmap = nano.Electron_vidNestedWPBitmap()[iel];
      if (!idElectron_noIso(bitmap, 1)) continue;
      bool isBarrel = fabs(eta) <= 1.479f;
      if ((isBarrel && fabs(dz)>0.10f) || (!isBarrel && fabs(dz)>0.20f)) continue;
      if ((isBarrel && fabs(dxy)>0.05f) || (!isBarrel && fabs(dxy)>0.10f)) continue; 
      isSignal = IsSignal(nano, iel, isZgamma);
      id = idElectron_noIso(bitmap,3);
    }
    pico.out_el_pt().push_back(scaleres_corr*pt);
    pico.out_el_pt_raw().push_back(pt);
    pico.out_el_energyErr().push_back(nano.Electron_energyErr()[iel]);
    pico.out_el_eta().push_back(eta);
    pico.out_el_etasc().push_back(etasc);
    pico.out_el_phi().push_back(phi);
    pico.out_el_miniso().push_back(miniiso);
    pico.out_el_reliso().push_back(nano.Electron_pfRelIso03_all()[iel]);
    pico.out_el_dz().push_back(dz);
    pico.out_el_dxy().push_back(dxy);
    pico.out_el_ip3d().push_back(nano.Electron_ip3d()[iel]);
    pico.out_el_id().push_back(id);
    pico.out_el_sig().push_back(isSignal);
    pico.out_el_ispf().push_back(nano.Electron_isPFcand()[iel]);
    pico.out_el_charge().push_back(nano.Electron_charge()[iel]);
    if (!isData) {
      pico.out_el_pflavor().push_back(nano.Electron_genPartFlav()[iel]);
    }
    
    // veto electron selection
    if (miniiso < ElectronMiniIsoCut) {
      pico.out_nvel()++;
      pico.out_nvlep()++;
      // save indices of matching jets
      for (int ijet(0); ijet<nano.nJet(); ijet++) {
        if (dR(eta, nano.Jet_eta()[ijet], phi, nano.Jet_phi()[ijet])<0.4f &&
            fabs(Jet_pt[ijet] - nano.Electron_pt()[iel])/nano.Electron_pt()[iel] < 1.0f)
          jet_isvlep_nano_idx.push_back(ijet);
      }
    }

    for (int iph(0); iph<nano.nPhoton(); ++iph) {
      if (Photon_electronIdx[iph]==iel)
        photon_el_pico_idx[iph] = pico.out_el_pt().size()-1;
    }

    if (isSignal) {
      pico.out_nel()++;
      pico.out_nlep()++;
      sig_el_nano_idx.push_back(iel);
      sig_el_pico_idx.push_back(pico_idx);
      // save indices of matching jets
      for (int ijet(0); ijet<nano.nJet(); ijet++) {
        if (dR(eta, nano.Jet_eta()[ijet], phi, nano.Jet_phi()[ijet])<jetDRCut &&
            fabs(Jet_pt[ijet] - nano.Electron_pt()[iel])/nano.Electron_pt()[iel] < jetpTCut)
          jet_islep_nano_idx.push_back(ijet);
      }
    }
    pico_idx++;
  }
  return sig_el_nano_idx;
}

bool ElectronProducer::idElectron_noIso(int bitmap, int level){
  // decision for each cut represented by 3 bits (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
  // Electron_vidNestedWPBitmap 
  //0 - MinPtCut
  //1 - GsfEleSCEtaMultiRangeCut
  //2 - GsfEleDEtaInSeedCut
  //3 - GsfEleDPhiInCut
  //4 - GsfEleFull5x5SigmaIEtaIEtaCut
  //5 - GsfEleHadronicOverEMEnergyScaledCut
  //6 - GsfEleEInverseMinusPInverseCut
  //7 - GsfEleRelPFIsoScaledCut
  //8 - GsfEleConversionVetoCut
  //9 - GsfEleMissingHitsCut
  bool pass = true;
  // cout<<std::bitset<8*sizeof(bitmap)>(bitmap)<<endl;
  for (int i(0); i<10; i++){
    if (i==7) continue;
    if ( ((bitmap >> i*3) & 0x7) < level) pass = false;
  }
  return pass;
}

float ElectronProducer::ConvertMVA(float mva_mini) {
  // 2.0 / (1.0 + exp(-2.0 * response)) - 1)
  float mva_nano = 2.0 / (1.0 + exp(-2.0 * mva_mini)) - 1;
  return mva_nano;
}

bool ElectronProducer::HzzId_WP2022(float pt, float etasc, float hzzmvaid) {
  //2022 WPs for 2018 ID training taken from https://indico.cern.ch/event/1429005/contributions/6039535/attachments/2891374/5077286/240712_H4lrun3_Approval.pdf
  if (pt < 10.0f) {
    if (fabs(etasc) < 0.8f) {
     return (hzzmvaid > ConvertMVA(1.6339));
    }
    else if (fabs(etasc) < 1.479f) {
      return (hzzmvaid > ConvertMVA(1.5499));
    }
    else {
      return (hzzmvaid > ConvertMVA(2.0629));
    }
  }
  else {
    if (fabs(etasc) < 0.8f) {
      return (hzzmvaid > ConvertMVA(0.3685));
    }
    else if (fabs(etasc) < 1.479f) {
      return (hzzmvaid > ConvertMVA(0.2662));
    }
    else {
      return (hzzmvaid > ConvertMVA(-0.5444));
    }
  }
}

bool ElectronProducer::EcalDriven(int bitmap){
  // decision for each cut represented by 1 bit
  //0 - MinPtCut
  //1 - GsfEleSCEtaMultiRangeCut
  //2 - GsfEleDEtaInSeedCut
  //3 - GsfEleDPhiInCut
  //4 - GsfEleFull5x5SigmaIEtaIEtaWithSatCut
  //5 - GsfEleFull5x5E2x5OverE5x5WithSatCut
  //6 - GsfEleHadronicOverEMLinearCut
  //7 - GsfEleDxyCut
  //8 - GsfEleMissingHitsCut
  //9 - GsfEleEcalDrivenCut
  bool ecaldriven = bitmap >> 9;
  return ecaldriven;
}
//Branches change in run3

