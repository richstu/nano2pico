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
      if (nano.Photon_pt().at(iph) > SignalPhotonPtCut) {
        pico.out_ph_pt().push_back(nano.Photon_pt().at(iph));
        pico.out_nph()++;
      }
    }

    return;
}
