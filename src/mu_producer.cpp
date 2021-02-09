#include "mu_producer.hpp"

#include "utilities.hpp"

using namespace std;

MuonProducer::MuonProducer(int year_, bool isData_){
    year = year_;
    isData = isData_;
}

MuonProducer::~MuonProducer(){
}

vector<int> MuonProducer::WriteMuons(nano_tree &nano, pico_tree &pico, vector<int> &jet_islep_nano_idx, vector<int> &jet_isvlep_nano_idx, vector<int> &sig_mu_pico_idx, bool isZgamma, bool isFastsim){
  vector<float> Jet_pt, Jet_mass;
  getJetWithJEC(nano, isFastsim, Jet_pt, Jet_mass);

  vector<int> sig_mu_nano_idx;
  pico.out_nmu() = 0; pico.out_nvmu() = 0;
  int pico_idx = 0;
  for(int imu(0); imu<nano.nMuon(); ++imu){
    float pt = nano.Muon_pt()[imu];
    float eta = nano.Muon_eta()[imu];
    bool isSignal = false;
    if(isZgamma) { // For Zgamma productions
      if (pt <= ZgMuonPtCut) continue;
      if (fabs(eta) > MuonEtaCut) continue;
      if (fabs(nano.Muon_dz()[imu])>1.0)  continue;
      if (fabs(nano.Muon_dxy()[imu])>0.5) continue; 
      if ((nano.Muon_looseId()[imu] || (pt > 200 && nano.Muon_highPtId()[imu])) && 
           nano.Muon_pfRelIso03_all()[imu] < MuonRelIsoCut &&
           nano.Muon_sip3d()[imu] < 4)
        isSignal = true;
      pico.out_mu_sip3d().push_back(nano.Muon_sip3d()[imu]);
    }
    else {
      if (!nano.Muon_mediumId()[imu]) continue;
      if (pt <= VetoMuonPtCut) continue;
      if (fabs(eta) > MuonEtaCut) continue;
      if (pt > SignalMuonPtCut &&
        nano.Muon_miniPFRelIso_all()[imu] < MuonMiniIsoCut &&
        fabs(nano.Muon_dz()[imu])<=0.5 &&
        fabs(nano.Muon_dxy()[imu])<=0.2)
        isSignal = true;
    }
    pico.out_mu_pt().push_back(pt);
    pico.out_mu_eta().push_back(eta);
    pico.out_mu_phi().push_back(nano.Muon_phi()[imu]);
    pico.out_mu_miniso().push_back(nano.Muon_miniPFRelIso_all()[imu]);
    pico.out_mu_reliso().push_back(nano.Muon_pfRelIso03_all()[imu]);
    pico.out_mu_dz().push_back(nano.Muon_dz()[imu]);
    pico.out_mu_dxy().push_back(nano.Muon_dxy()[imu]);
    pico.out_mu_ip3d().push_back(nano.Muon_ip3d()[imu]);
    pico.out_mu_id().push_back(nano.Muon_looseId()[imu]);
    pico.out_mu_sig().push_back(isSignal);
    pico.out_mu_charge().push_back(nano.Muon_charge()[imu]);
    if (!isData)
      pico.out_mu_pflavor().push_back(nano.Muon_genPartFlav()[imu]);

    // veto muon selection
    if (nano.Muon_miniPFRelIso_all()[imu] < MuonMiniIsoCut && 
        fabs(nano.Muon_dz()[imu])<=0.5 &&
        fabs(nano.Muon_dxy()[imu])<=0.2) {
      pico.out_nvmu()++;
      pico.out_nvlep()++;
      // save indices of matching jets
      for (int ijet(0); ijet<nano.nJet(); ijet++) {
        if (dR(nano.Muon_eta()[imu], nano.Jet_eta()[ijet], nano.Muon_phi()[imu], nano.Jet_phi()[ijet])<0.4 &&
          fabs(Jet_pt[ijet] - nano.Muon_pt()[imu])/nano.Muon_pt()[imu] < 1)
          jet_isvlep_nano_idx.push_back(ijet);
      }
    }
    if (isSignal) {
      pico.out_nmu()++;
      pico.out_nlep()++;
      sig_mu_nano_idx.push_back(imu);
      sig_mu_pico_idx.push_back(pico_idx);

      // save indices of matching jets
      for (int ijet(0); ijet<nano.nJet(); ijet++) {
        if (dR(nano.Muon_eta()[imu], nano.Jet_eta()[ijet], nano.Muon_phi()[imu], nano.Jet_phi()[ijet])<0.4 &&
          fabs(Jet_pt[ijet] - nano.Muon_pt()[imu])/nano.Muon_pt()[imu] < 1)
          jet_islep_nano_idx.push_back(ijet);
      }
    }
    pico_idx++;
  }
  return sig_mu_nano_idx;
}
