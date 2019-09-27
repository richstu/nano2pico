#include "mu_producer.hpp"

#include "utilities.hpp"

using namespace std;

MuonProducer::MuonProducer(int year_){
    year = year_;
}

MuonProducer::~MuonProducer(){
}

void MuonProducer::WriteMuons(nano_tree &nano, pico_tree &pico){
    pico.out_nmu() = 0;
    for(int imu(0); imu<nano.nMuon(); ++imu){
      if (nano.Muon_pt().at(imu) > SignalMuonPtCut) {
        pico.out_mu_pt().push_back(nano.Muon_pt().at(imu));
        pico.out_nmu()++;
      }
    }

    return;
}
