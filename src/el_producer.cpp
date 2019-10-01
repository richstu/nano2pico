#include "el_producer.hpp"

#include "utilities.hpp"

using namespace std;

ElectronProducer::ElectronProducer(int year_){
    year = year_;
}

ElectronProducer::~ElectronProducer(){
}

vector<int> ElectronProducer::WriteElectrons(nano_tree &nano, pico_tree &pico){
    
  vector<int> sig_el_nano_idx;
  pico.out_nel() = 0; pico.out_nvel() = 0;
  for(int iel(0); iel<nano.nElectron(); ++iel){
      //cut-based ID Fall17 V2 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
    if (nano.Electron_cutBased()[iel] == 0) continue;
    if (nano.Electron_pt()[iel] <= VetoElectronPtCut) continue;
    if (fabs(nano.Electron_deltaEtaSC()[iel] + nano.Electron_eta()[iel]) > ElectronEtaCut) continue;
    if (nano.Electron_miniPFRelIso_all()[iel]==ElectronMiniIsoCut) continue;

    bool isSig = nano.Electron_cutBased()[iel]>=3 && nano.Electron_pt()[iel] > SignalElectronPtCut;

    pico.out_el_pt().push_back(nano.Electron_pt()[iel]);
    pico.out_el_eta().push_back(nano.Electron_eta()[iel]);
    pico.out_el_phi().push_back(nano.Electron_phi()[iel]);
    pico.out_el_miniso().push_back(nano.Electron_miniPFRelIso_all()[iel]);
    pico.out_el_reliso().push_back(nano.Electron_pfRelIso03_all()[iel]);
    pico.out_el_dz().push_back(nano.Electron_dz()[iel]);
    pico.out_el_d0().push_back(nano.Electron_dxy()[iel]);
    pico.out_el_ip3d().push_back(nano.Electron_ip3d()[iel]);
    pico.out_el_sig().push_back(isSig);
    pico.out_el_ispf().push_back(nano.Electron_isPFcand()[iel]);
    pico.out_el_charge().push_back(nano.Electron_charge()[iel]);

      //this will be filled in mc_producer??
    pico.out_el_tm().push_back(false);
    
    pico.out_nvel()++;
    if (isSig) {
      pico.out_nel()++;
      sig_el_nano_idx.push_back(iel);
    }
  }
  return sig_el_nano_idx;
}
