#include "el_producer.hpp"

#include "utilities.hpp"

using namespace std;

ElectronProducer::ElectronProducer(int year_){
    year = year_;
}

ElectronProducer::~ElectronProducer(){
}

vector<int> ElectronProducer::WriteElectrons(nano_tree &nano, pico_tree &pico, vector<int> &jet_islep_nano_idx){
    
  vector<int> sig_el_nano_idx;
  pico.out_nel() = 0; pico.out_nvel() = 0;
  for(int iel(0); iel<nano.nElectron(); ++iel){
      
    float pt = nano.Electron_pt()[iel];///nano.Electron_eCorr()[iel]; 
    float etasc = nano.Electron_deltaEtaSC()[iel] + nano.Electron_eta()[iel];

    if (pt <= VetoElectronPtCut) continue;
    if (fabs(etasc) > ElectronEtaCut) continue;

    int bitmap = nano.Electron_vidNestedWPBitmap()[iel];
    if (!idElectron_noIso(bitmap, 1)) continue;

    bool isBarrel = fabs(etasc) <= 1.479;
    if ((isBarrel && fabs(nano.Electron_dz()[iel])>0.10) || (!isBarrel && fabs(nano.Electron_dz()[iel])>0.20)) continue;
    if ((isBarrel && fabs(nano.Electron_dxy()[iel])>0.05) || (!isBarrel && fabs(nano.Electron_dxy()[iel])>0.10)) continue; 

    bool isSignal = false;
    if (idElectron_noIso(bitmap, 3) && 
        nano.Electron_miniPFRelIso_all()[iel] < ElectronMiniIsoCut &&
        pt > SignalElectronPtCut)
      isSignal = true;

    pico.out_el_pt().push_back(pt);
    pico.out_el_eta().push_back(nano.Electron_eta()[iel]);
    pico.out_el_phi().push_back(nano.Electron_phi()[iel]);
    pico.out_el_miniso().push_back(nano.Electron_miniPFRelIso_all()[iel]);
    pico.out_el_reliso().push_back(nano.Electron_pfRelIso03_all()[iel]);
    pico.out_el_dz().push_back(nano.Electron_dz()[iel]);
    pico.out_el_dxy().push_back(nano.Electron_dxy()[iel]);
    pico.out_el_ip3d().push_back(nano.Electron_ip3d()[iel]);
    pico.out_el_sig().push_back(isSignal);
    pico.out_el_ispf().push_back(nano.Electron_isPFcand()[iel]);
    pico.out_el_charge().push_back(nano.Electron_charge()[iel]);
    
    pico.out_el_tm().push_back(nano.Electron_genPartIdx()[iel]!=-1);
    
    pico.out_nvel()++;
    if (isSignal) {
      pico.out_nel()++;
      sig_el_nano_idx.push_back(iel);

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
    if (i==5) continue;
    if ( ((bitmap >> i*3) & 0x7) < level) pass = false;
  }
  return pass;
}
