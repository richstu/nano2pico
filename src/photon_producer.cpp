#include "photon_producer.hpp"

#include "utilities.hpp"

using namespace std;

PhotonProducer::PhotonProducer(int year_){
    year = year_;
}

PhotonProducer::~PhotonProducer(){
}

vector<int> PhotonProducer::WritePhotons(nano_tree &nano, pico_tree &pico, vector<int> &jet_isphoton_nano_idx, 
                                         vector<int> &sig_el_nano_idx, vector<int> &sig_mu_nano_idx){
  pico.out_nphoton() = 0; 
  vector<int> sig_photon_nano_idx;
  for(int iph(0); iph<nano.nPhoton(); ++iph){
    if (nano.Photon_pt()[iph] <= PhotonPtCut) continue;
    if (fabs(nano.Photon_eta()[iph]) > PhotonEtaCut) continue;
    if (!(nano.Photon_isScEtaEB()[iph] || nano.Photon_isScEtaEE()[iph])) continue;
    // if (nano.Photon_pfRelIso03_all()[iph]==PhotonRelIsoCut) continue; // no isolation cut in 2016...?

    pico.out_photon_pt().push_back(nano.Photon_pt()[iph]);
    pico.out_photon_eta().push_back(nano.Photon_eta()[iph]);
    pico.out_photon_phi().push_back(nano.Photon_phi()[iph]);
    pico.out_photon_reliso().push_back(nano.Photon_pfRelIso03_all()[iph]);
    pico.out_photon_r9().push_back(nano.Photon_r9()[iph]);
    pico.out_photon_pflavor().push_back(nano.Photon_genPartFlav()[iph]);
    pico.out_photon_elveto().push_back(nano.Photon_electronVeto()[iph]);
    pico.out_photon_id().push_back(nano.Photon_mvaID_WP90()[iph]);
    bool isSig = nano.Photon_mvaID_WP90()[iph] && nano.Photon_electronVeto()[iph];
    pico.out_photon_sig().push_back(isSig);
    // Find min(dR) between photon and signal lepton
    double minLepDR(999.);
    for(size_t iel(0); iel<sig_el_nano_idx.size(); iel++) {
      double tempDR = dR(nano.Photon_eta()[iph], nano.Electron_eta()[sig_el_nano_idx.at(iel)], 
                         nano.Photon_phi()[iph], nano.Electron_phi()[sig_el_nano_idx.at(iel)]);
      if(tempDR < minLepDR) minLepDR = tempDR;
    }
    for(size_t imu(0); imu<sig_mu_nano_idx.size(); imu++) {
      double tempDR = dR(nano.Photon_eta()[iph], nano.Muon_eta()[sig_mu_nano_idx.at(imu)], 
                         nano.Photon_phi()[iph], nano.Muon_phi()[sig_mu_nano_idx.at(imu)]);
      if(tempDR < minLepDR) minLepDR = tempDR;
    }
    pico.out_photon_drmin().push_back(minLepDR);
    if(isSig) {
      pico.out_nphoton()++;
      sig_photon_nano_idx.push_back(iph);
      // save indices of matching jets
      if (nano.Photon_jetIdx()[iph]>=0) 
        jet_isphoton_nano_idx.push_back(nano.Photon_jetIdx()[iph]);
      else 
        for (int ijet(0); ijet<nano.nJet(); ijet++) 
          if (dR(nano.Photon_eta()[iph], nano.Jet_eta()[ijet], nano.Photon_phi()[iph], nano.Jet_phi()[ijet])<0.4)
            jet_isphoton_nano_idx.push_back(ijet);
    }
  }
  return sig_photon_nano_idx;
}
