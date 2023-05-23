#include "el_producer.hpp"

#include "utilities.hpp"

#include<iomanip>

using namespace std;

ElectronProducer::ElectronProducer(int year_, bool isData_){
    year = year_;
    isData = isData_;
}

ElectronProducer::~ElectronProducer(){
}

vector<int> ElectronProducer::WriteElectrons(nano_tree &nano, pico_tree &pico, vector<int> &jet_islep_nano_idx, vector<int> &jet_isvlep_nano_idx, vector<int> &sig_el_pico_idx, bool isZgamma, bool isFastsim){
  vector<float> Jet_pt, Jet_mass;
  getJetWithJEC(nano, isFastsim, Jet_pt, Jet_mass);

  vector<int> sig_el_nano_idx;
  pico.out_nel() = 0; pico.out_nvel() = 0;
  int pico_idx = 0;
  for(int iel(0); iel<nano.nElectron(); ++iel){
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
      if (fabs(dz) > 1.0)  continue;
      if (fabs(dxy) > 0.5) continue; 
      switch(year) {
        case 2016:
        case 2017:
        case 2018:
          isSignal = nano.Electron_mvaFall17V2Iso_WPL()[iel];
          pico.out_el_idmva().push_back(nano.Electron_mvaFall17V2Iso()[iel]);
          pico.out_el_sip3d().push_back(nano.Electron_sip3d()[iel]);
          pico.out_el_phidx().push_back(nano.Electron_photonIdx()[iel]);
          pico.out_el_id80().push_back(nano.Electron_mvaFall17V2Iso_WP80()[iel]);
          pico.out_el_id90().push_back(nano.Electron_mvaFall17V2Iso_WP90()[iel]);
          pico.out_el_idLoose().push_back(nano.Electron_mvaFall17V2Iso_WPL()[iel]);
          pico.out_el_etPt().push_back(nano.Electron_scEtOverPt()[iel]);
          pico.out_el_eminusp().push_back(nano.Electron_eInvMinusPInv()[iel]);
          break;
        case 2022:
          isSignal = nano.Electron_mvaIso_WP80()[iel];
          pico.out_el_idmva().push_back(nano.Electron_mvaIso()[iel]);
          pico.out_el_sip3d().push_back(nano.Electron_sip3d()[iel]);
          pico.out_el_phidx().push_back(nano.Electron_photonIdx()[iel]);
          pico.out_el_id80().push_back(nano.Electron_mvaIso_WP80()[iel]);
          pico.out_el_id90().push_back(nano.Electron_mvaIso_WP90()[iel]);
          pico.out_el_idLoose().push_back(nano.Electron_mvaIso_WP80()[iel]);
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
      bool isBarrel = fabs(eta) <= 1.479;
      if ((isBarrel && fabs(dz)>0.10) || (!isBarrel && fabs(dz)>0.20)) continue;
      if ((isBarrel && fabs(dxy)>0.05) || (!isBarrel && fabs(dxy)>0.10)) continue; 
      id = idElectron_noIso(bitmap,3);
      if (id && 
          miniiso < ElectronMiniIsoCut &&
          pt > SignalElectronPtCut)
        isSignal = true;
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
        if (dR(eta, nano.Jet_eta()[ijet], phi, nano.Jet_phi()[ijet])<0.4 &&
            fabs(Jet_pt[ijet] - nano.Electron_pt()[iel])/nano.Electron_pt()[iel] < 1)
          jet_isvlep_nano_idx.push_back(ijet);
      }
    }

    if (isSignal) {
      pico.out_nel()++;
      pico.out_nlep()++;
      sig_el_nano_idx.push_back(iel);
      sig_el_pico_idx.push_back(pico_idx);
      // save indices of matching jets
      for (int ijet(0); ijet<nano.nJet(); ijet++) {
        if (dR(eta, nano.Jet_eta()[ijet], phi, nano.Jet_phi()[ijet])<0.4 &&
            fabs(Jet_pt[ijet] - nano.Electron_pt()[iel])/nano.Electron_pt()[iel] < 1)
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

