#include "ph_producer.hpp"

#include "utilities.hpp"

using namespace std;

PhotonProducer::PhotonProducer(int year_){
    year = year_;
}

PhotonProducer::~PhotonProducer(){
}

void PhotonProducer::WritePhotons(nano_tree &nano, pico_tree &pico){
    pico.out_nph() = 0; 
    for(int iph(0); iph<nano.nPhoton(); ++iph){
      if (!nano.Photon_mvaID_WP90()[iph]) continue;
      if (!nano.Photon_electronVeto()[iph]) continue;
      if (nano.Photon_pt()[iph] <= PhotonPtCut) continue;
      if (fabs(nano.Photon_eta()[iph]) > PhotonEtaCut) continue;
      if ((1.4442 <= fabs(nano.Photon_eta()[iph])) && (fabs(nano.Photon_eta()[iph]) <= 1.566)) continue; // exclude gap
      // if (nano.Photon_pfRelIso03_all()[iph]==PhotonRelIsoCut) continue; // no isolation cut in 2016...

      pico.out_ph_pt().push_back(nano.Photon_pt()[iph]);
      pico.out_ph_eta().push_back(nano.Photon_eta()[iph]);
      pico.out_ph_phi().push_back(nano.Photon_phi()[iph]);
      pico.out_ph_reliso().push_back(nano.Photon_pfRelIso03_all()[iph]);
      pico.out_ph_r9().push_back(nano.Photon_r9()[iph]);

      //this will be filled in mc_producer??
      pico.out_ph_tm().push_back(false);
      
      pico.out_nph()++;
    }

    return;
}
