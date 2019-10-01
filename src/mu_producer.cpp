#include "mu_producer.hpp"

#include "utilities.hpp"

using namespace std;

MuonProducer::MuonProducer(int year_){
    year = year_;
}

MuonProducer::~MuonProducer(){
}

vector<int> MuonProducer::WriteMuons(nano_tree &nano, pico_tree &pico){

  vector<int> sig_mu_nano_idx;
  pico.out_nmu() = 0; pico.out_nvmu() = 0;
  for(int imu(0); imu<nano.nMuon(); ++imu){
    if (!nano.Muon_mediumId()[imu]) continue;
    if (nano.Muon_pt()[imu] <= VetoMuonPtCut) continue;
    if (fabs(nano.Muon_eta()[imu]) > MuonEtaCut) continue;
    if (nano.Muon_miniPFRelIso_all()[imu]==MuonMiniIsoCut) continue;

    bool isSig = nano.Muon_pt()[imu] > SignalMuonPtCut;

    pico.out_mu_pt().push_back(nano.Muon_pt()[imu]);
    pico.out_mu_eta().push_back(nano.Muon_eta()[imu]);
    pico.out_mu_phi().push_back(nano.Muon_phi()[imu]);
    pico.out_mu_miniso().push_back(nano.Muon_miniPFRelIso_all()[imu]);
    pico.out_mu_reliso().push_back(nano.Muon_pfRelIso03_all()[imu]);
    pico.out_mu_dz().push_back(nano.Muon_dz()[imu]);
    pico.out_mu_d0().push_back(nano.Muon_dxy()[imu]);
    pico.out_mu_ip3d().push_back(nano.Muon_ip3d()[imu]);
    pico.out_mu_sig().push_back(isSig);
    pico.out_mu_charge().push_back(nano.Muon_charge()[imu]);

      //this will be filled in mc_producer??
    pico.out_mu_tm().push_back(false);

    pico.out_nvmu()++;
    if (isSig) {
      pico.out_nmu()++;
      sig_mu_nano_idx.push_back(imu);
    }
  }
  return sig_mu_nano_idx;
}
