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

vector<int> ElectronProducer::WriteElectrons(nano_tree &nano, pico_tree &pico, vector<int> &jet_islep_nano_idx, vector<int> &sig_el_pico_idx, bool isZgamma){
  vector<int> sig_el_nano_idx;
  pico.out_nel() = 0; pico.out_nvel() = 0;
  int pico_idx = 0;
  for(int iel(0); iel<nano.nElectron(); ++iel){
    float pt = nano.Electron_pt()[iel];///nano.Electron_eCorr()[iel]; 
    float etasc = nano.Electron_deltaEtaSC()[iel] + nano.Electron_eta()[iel];
    bool isSignal = false;
    bool id = false;
    if(isZgamma) { // For Zgamma productions
      if (pt <= ZgElectronPtCut) continue;
      if (fabs(etasc) > ElectronEtaCut) continue;
      if (fabs(nano.Electron_dz()[iel])>1.0)  continue;
      if (fabs(nano.Electron_dxy()[iel])>0.5) continue; 
      id = nano.Electron_mvaFall17V2Iso_WP90()[iel];
      double wp[2][3] = {{-0.145237, -0.0315746, -0.032173},
                         { 0.604775,  0.628743,   0.896462}};
      int ipt(1), ieta(2);
      if(pt>10) ipt = 0;
      if(fabs(etasc) < 0.8) ieta = 0;
      else if(fabs(etasc) < 1.479) ieta = 1;
      double mva = nano.Electron_mvaFall17V2Iso()[iel];
      isSignal = mva > wp[ipt][ieta];
      pico.out_el_idmva().push_back(mva);
      pico.out_el_sip3d().push_back(nano.Electron_sip3d()[iel]);
    }
    else {
      if (pt <= VetoElectronPtCut) continue;
      if (fabs(etasc) > ElectronEtaCut) continue;
      int bitmap = nano.Electron_vidNestedWPBitmap()[iel];
      if (!idElectron_noIso(bitmap, 1)) continue;
      bool isBarrel = fabs(etasc) <= 1.479;
      if ((isBarrel && fabs(nano.Electron_dz()[iel])>0.10) || (!isBarrel && fabs(nano.Electron_dz()[iel])>0.20)) continue;
      if ((isBarrel && fabs(nano.Electron_dxy()[iel])>0.05) || (!isBarrel && fabs(nano.Electron_dxy()[iel])>0.10)) continue; 
      id = idElectron_noIso(bitmap,3);
      if (id && 
          nano.Electron_miniPFRelIso_all()[iel] < ElectronMiniIsoCut &&
          pt > SignalElectronPtCut)
        isSignal = true;
    }
    pico.out_el_pt().push_back(pt);
    pico.out_el_eta().push_back(nano.Electron_eta()[iel]);
    pico.out_el_phi().push_back(nano.Electron_phi()[iel]);
    pico.out_el_miniso().push_back(nano.Electron_miniPFRelIso_all()[iel]);
    pico.out_el_reliso().push_back(nano.Electron_pfRelIso03_all()[iel]);
    pico.out_el_dz().push_back(nano.Electron_dz()[iel]);
    pico.out_el_dxy().push_back(nano.Electron_dxy()[iel]);
    pico.out_el_ip3d().push_back(nano.Electron_ip3d()[iel]);
    pico.out_el_id().push_back(id);
    pico.out_el_sig().push_back(isSignal);
    pico.out_el_ispf().push_back(nano.Electron_isPFcand()[iel]);
    pico.out_el_charge().push_back(nano.Electron_charge()[iel]);
    if (!isData) {
      pico.out_el_pflavor().push_back(nano.Electron_genPartFlav()[iel]);
    }
    
    if (nano.Electron_miniPFRelIso_all()[iel] < ElectronMiniIsoCut) {
      pico.out_nvel()++;
      pico.out_nvlep()++;
    }
    // cout<<"Adding to"<<(isSignal ? "signal":"veto")<<" leptons"<<iel<<": pt = "<<setw(10)<<nano.Electron_pt()[iel]
    //                                 <<" eta = "<<setw(10)<<nano.Electron_eta()[iel]
    //                                 <<" phi = "<<setw(10)<<nano.Electron_phi()[iel]
    //                                 <<" m = "<<setw(10)<<nano.Electron_mass()[iel];
    if (isSignal) {
      pico.out_nel()++;
      pico.out_nlep()++;
      sig_el_nano_idx.push_back(iel);
      sig_el_pico_idx.push_back(pico_idx);
      // save indices of matching jets
      if (nano.Electron_isPFcand()[iel] && nano.Electron_jetIdx()[iel]>=0) {
        jet_islep_nano_idx.push_back(nano.Electron_jetIdx()[iel]);
      } else {
        for (int ijet(0); ijet<nano.nJet(); ijet++) {
          if (dR(nano.Electron_eta()[iel], nano.Jet_eta()[ijet], nano.Electron_phi()[iel], nano.Jet_phi()[ijet])<0.4)
            jet_islep_nano_idx.push_back(ijet);
        }
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
