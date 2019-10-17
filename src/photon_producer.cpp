#include "photon_producer.hpp"

#include "utilities.hpp"

using namespace std;

PhotonProducer::PhotonProducer(int year_){
    year = year_;
}

PhotonProducer::~PhotonProducer(){
}

vector<int> PhotonProducer::WritePhotons(nano_tree &nano, pico_tree &pico, vector<int> &jet_isphoton_nano_idx){
  pico.out_nphoton() = 0; 
  vector<int> sig_photon_nano_idx;
  for(int iph(0); iph<nano.nPhoton(); ++iph){
    if (!nano.Photon_mvaID_WP90()[iph]) continue;
    if (!nano.Photon_electronVeto()[iph]) continue;
    if (nano.Photon_pt()[iph] <= PhotonPtCut) continue;
    if (fabs(nano.Photon_eta()[iph]) > PhotonEtaCut) continue;
    if ((1.4442 <= fabs(nano.Photon_eta()[iph])) && (fabs(nano.Photon_eta()[iph]) <= 1.566)) continue; // exclude gap
    // if (nano.Photon_pfRelIso03_all()[iph]==PhotonRelIsoCut) continue; // no isolation cut in 2016...?

    pico.out_photon_pt().push_back(nano.Photon_pt()[iph]);
    pico.out_photon_eta().push_back(nano.Photon_eta()[iph]);
    pico.out_photon_phi().push_back(nano.Photon_phi()[iph]);
    pico.out_photon_reliso().push_back(nano.Photon_pfRelIso03_all()[iph]);
    pico.out_photon_r9().push_back(nano.Photon_r9()[iph]);
    pico.out_photon_pflavor().push_back(nano.Photon_genPartFlav()[iph]);
    
    pico.out_nphoton()++;
    sig_photon_nano_idx.push_back(iph);
    // save indices of matching jets
    if (nano.Photon_jetIdx()[iph]>=0) {
      jet_isphoton_nano_idx.push_back(nano.Photon_jetIdx()[iph]);
    } else {
      for (int ijet(0); ijet<nano.nJet(); ijet++) {
        if (dR(nano.Photon_eta()[iph], nano.Jet_eta()[ijet], nano.Photon_phi()[iph], nano.Jet_phi()[ijet])<0.4)
          jet_isphoton_nano_idx.push_back(ijet);
      }
    }
  }

  return sig_photon_nano_idx;
}
