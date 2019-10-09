#include "tk_producer.hpp"

#include "utilities.hpp"

using namespace std;

IsoTrackProducer::IsoTrackProducer(int year_){
    year = year_;
}

IsoTrackProducer::~IsoTrackProducer(){
}

void IsoTrackProducer::WriteIsoTracks(nano_tree &nano, pico_tree &pico, 
                                      vector<int> &sig_el_nano_idx, vector<int> &sig_mu_nano_idx){
  
  pico.out_ntk() = 0;
  // N.B. Objects that end up in the slimmedElecrtons or slimmedMuons collections are not stored 
  // as IsoTrack in Nano so we have to loop over all three collections
  for (int itk(0); itk < nano.nIsoTrack(); itk++) {
    IsGoodTk(pico, nano.IsoTrack_pdgId()[itk],
             nano.IsoTrack_pt()[itk], nano.IsoTrack_eta()[itk], nano.IsoTrack_phi()[itk], 
             nano.IsoTrack_miniPFRelIso_chg()[itk], nano.IsoTrack_pfRelIso03_chg()[itk],
             nano.IsoTrack_dxy()[itk], nano.IsoTrack_dz()[itk],
             GetMT(nano.MET_pt(), nano.MET_phi(),  nano.IsoTrack_pt()[itk], nano.IsoTrack_phi()[itk]));
  }

  // collect relevant Electrons
  for (int iel(0); iel < nano.nElectron(); iel++) {
    if (find(sig_el_nano_idx.begin(), sig_el_nano_idx.end(), iel)!=sig_el_nano_idx.end()) continue;
    IsGoodTk(pico, nano.Electron_pdgId()[iel],
             nano.Electron_pt()[iel], nano.Electron_eta()[iel], nano.Electron_phi()[iel], 
             nano.Electron_miniPFRelIso_chg()[iel], nano.Electron_pfRelIso03_chg()[iel],
             nano.Electron_dxy()[iel], nano.Electron_dz()[iel],
             GetMT(nano.MET_pt(), nano.MET_phi(),  nano.Electron_pt()[iel], nano.Electron_phi()[iel]));
  }
  
  // collect relevant Muons
  for (int imu(0); imu < nano.nMuon(); imu++) {
    if (find(sig_mu_nano_idx.begin(), sig_mu_nano_idx.end(), imu)!=sig_mu_nano_idx.end()) continue;
    IsGoodTk(pico, nano.Muon_pdgId()[imu],
             nano.Muon_pt()[imu], nano.Muon_eta()[imu], nano.Muon_phi()[imu], 
             nano.Muon_miniPFRelIso_chg()[imu], nano.Muon_pfRelIso03_chg()[imu],
             nano.Muon_dxy()[imu], nano.Muon_dz()[imu],
             GetMT(nano.MET_pt(), nano.MET_phi(),  nano.Muon_pt()[imu], nano.Muon_phi()[imu]));
  }
  return;
}

bool IsoTrackProducer::IsGoodTk(pico_tree &pico, int pdgid, float pt, float eta, float phi, 
                                float miniso_chg, float reliso_chg, float dxy, float dz, float mt){
  
  // re-applying cuts used in NanoAOD so that they would also apply to "tracks" in the lepton collections
  // ((pt>5 && (abs(pdgId) == 11 || abs(pdgId) == 13)) || pt > 10) && 
  // (abs(pdgId) < 15 || abs(eta) < 2.5) && 
  // abs(dxy) < 0.2 && 
  // abs(dz) < 0.1 && 
  // ((pfIsolationDR03().chargedHadronIso < 5 && pt < 25) || pfIsolationDR03().chargedHadronIso/pt < 0.2)
  
  pdgid = abs(pdgid);

  if (pdgid!=11 && pdgid!=13 && pdgid!=15 && pdgid!=211) return false; // allow also taus compared to babymaker
  
  if (pdgid==11 || pdgid==13) {
    if (pt < 5) return false;
    if ((pt < 25 && reliso_chg*pt > 5) || reliso_chg < 0.2) return false; //combination of our old cuts & Nano baseline
  } else {
    if (pt < 10) return false;
    if ((pt < 25 && reliso_chg*pt > 5) || reliso_chg < 0.1) return false; //combination of our old cuts & Nano baseline
  }

  if (fabs(eta) > 2.5) return false; // not applied to all tracks in Nano
  if (fabs(dxy)  > 0.2) return false; //applied to tracks but not leptons in Nano
  if (fabs(dz)  > 0.1) return false; //applied to tracks but not leptons in Nano
  // if (mt > 100) return false; // we should revisit whether this is useful

  pico.out_tk_pdgid().push_back(pdgid);
  pico.out_tk_pt().push_back(pt);
  pico.out_tk_eta().push_back(eta);
  pico.out_tk_phi().push_back(phi);
  pico.out_tk_dxy().push_back(dxy);
  pico.out_tk_miniso_chg().push_back(miniso_chg);
  pico.out_tk_reliso_chg().push_back(reliso_chg);
  pico.out_tk_dxy().push_back(dxy);
  pico.out_tk_dz().push_back(dz);
  pico.out_tk_mt().push_back(mt);

  // filled by mc_producer
  pico.out_tk_tm().push_back(false);    

  pico.out_ntk()++;
  return true;
}
