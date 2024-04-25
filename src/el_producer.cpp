#include "el_producer.hpp"

#include "correction.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <iomanip>
#include <memory>
#include <string>

using namespace std;

ElectronProducer::ElectronProducer(int year_, bool isData_, bool preVFP, float nanoaod_version_){
  year = year_;
  isData = isData_;
  if (year==2016 && preVFP) {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2016preVFP_UL/EGM_ScaleUnc.json");
    str_scale_syst_ = "2016preVFP";
  }
  else if (year==2016) {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2016postVFP_UL/EGM_ScaleUnc.json");
    str_scale_syst_ = "2016postVFP";
  }
  else if (year==2017) {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2017_UL/EGM_ScaleUnc.json");
    str_scale_syst_ = "2017";
  }
  else if (year==2018) {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2018_UL/EGM_ScaleUnc.json");
    str_scale_syst_ = "2018";
  }
  else {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2018_UL/EGM_ScaleUnc.json");
    str_scale_syst_ = "2018";
    std::cout << "WARNING: No dedicated EGM scale/smearing JSONs, defaulting to 2018" << std::endl;
  }
  map_scale_syst_ = cs_scale_syst_->at("UL-EGM_ScaleUnc");
  nanoaod_version = nanoaod_version_;
}

ElectronProducer::~ElectronProducer(){
}

bool ElectronProducer::IsSignal(nano_tree &nano, int nano_idx, bool isZgamma) {
  float pt = nano.Electron_pt()[nano_idx];
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
    switch(year) {
      case 2016:
      case 2017:
      case 2018:
        return nano.Electron_mvaFall17V2Iso_WPL()[nano_idx];
      case 2022:
      case 2023:
      default:
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
    if(isZgamma) { // For Zgamma productions
      if (pt <= ZgElectronPtCut) continue;
      if (fabs(etasc) > ElectronEtaCut) continue;
      if (fabs(dz) > dzCut)  continue;
      if (fabs(dxy) > dxyCut) continue; 
      isSignal = IsSignal(nano, iel, isZgamma);
      float scale_syst_up = 1.0;
      float scale_syst_dn = 1.0;
      if (year <= 2018) {
        scale_syst_up = map_scale_syst_->evaluate({str_scale_syst_,"scaleup",etasc,
            nano.Electron_seedGain()[iel]});
        scale_syst_dn = map_scale_syst_->evaluate({str_scale_syst_,"scaledown",etasc,
            nano.Electron_seedGain()[iel]});
      }

      switch(year) {
        case 2016:
        case 2017:
        case 2018:
          pico.out_el_idmva().push_back(nano.Electron_mvaFall17V2Iso()[iel]);
          pico.out_el_sip3d().push_back(nano.Electron_sip3d()[iel]);
          pico.out_el_phidx().push_back(Electron_photonIdx[iel]);
          pico.out_el_id80().push_back(nano.Electron_mvaFall17V2Iso_WP80()[iel]);
          pico.out_el_id90().push_back(nano.Electron_mvaFall17V2Iso_WP90()[iel]);
          pico.out_el_idLoose().push_back(nano.Electron_mvaFall17V2Iso_WPL()[iel]);
          pico.out_el_etPt().push_back(nano.Electron_scEtOverPt()[iel]);
          pico.out_el_eminusp().push_back(nano.Electron_eInvMinusPInv()[iel]);
          pico.out_sys_el_pt_resup().push_back(pt*(1.0+nano.Electron_dEsigmaUp()[iel]));
          pico.out_sys_el_pt_resdn().push_back(pt*(1.0+nano.Electron_dEsigmaDown()[iel]));
          pico.out_sys_el_pt_scaleup().push_back(pt*scale_syst_up);
          pico.out_sys_el_pt_scaledn().push_back(pt*scale_syst_dn);
          break;
        case 2022:
        case 2023:
          pico.out_el_idmva().push_back(nano.Electron_mvaIso()[iel]);
          pico.out_el_sip3d().push_back(nano.Electron_sip3d()[iel]);
          pico.out_el_phidx().push_back(Electron_photonIdx[iel]);
          pico.out_el_id80().push_back(nano.Electron_mvaIso_WP80()[iel]);
          pico.out_el_id90().push_back(nano.Electron_mvaIso_WP90()[iel]);
          pico.out_el_idLoose().push_back(nano.Electron_mvaIso_WP90()[iel]);
          pico.out_el_etPt().push_back(nano.Electron_scEtOverPt()[iel]);
          pico.out_el_eminusp().push_back(nano.Electron_eInvMinusPInv()[iel]);
          break;
        default:
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
    pico.out_el_pt().push_back(pt);
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
            fabs(Jet_pt[ijet] - nano.Electron_pt()[iel])/nano.Electron_pt()[iel] < 1.f)
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

