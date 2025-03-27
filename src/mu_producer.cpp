#include "mu_producer.hpp"


#include "utilities.hpp"

#include <algorithm>

using namespace std;

MuonProducer::MuonProducer(int year_, bool isData_, float nanoaod_version_, std::string rocco_file) :
  year(year_),
  isData(isData_),
  rc(rocco_file),
  rng(4357),
  nanoaod_version(nanoaod_version_){
    (void) year;
}

MuonProducer::~MuonProducer(){
}

bool MuonProducer::IsSignal(nano_tree &nano, int nano_idx, bool isZgamma) {
  float pt = nano.Muon_pt()[nano_idx];
  float eta = nano.Muon_eta()[nano_idx];
  if(isZgamma) { // For Zgamma productions
    if (pt <= ZgMuonPtCut) return false;
    if (fabs(eta) > MuonEtaCut) return false;
    if (fabs(nano.Muon_dz()[nano_idx])>dzCut)  return false;
    if (fabs(nano.Muon_dxy()[nano_idx])>dxyCut) return false; 
    if (nano.Muon_looseId()[nano_idx] && 
         nano.Muon_pfRelIso03_all()[nano_idx] < MuonRelIsoCut &&
         nano.Muon_sip3d()[nano_idx] < MuonSip3dCut)
      return true;
    return false;
  }
  else {
    if (!nano.Muon_mediumId()[nano_idx]) return false;
    if (pt <= VetoMuonPtCut) return false;
    if (fabs(eta) > MuonEtaCut) return false;
    if (pt > SignalMuonPtCut &&
      nano.Muon_miniPFRelIso_all()[nano_idx] < MuonMiniIsoCut &&
      fabs(nano.Muon_dz()[nano_idx])<=0.5f && 
      fabs(nano.Muon_dxy()[nano_idx])<=0.2f)
      return true;
    return false;
  }
  return false;
}

vector<int> MuonProducer::WriteMuons(nano_tree &nano, pico_tree &pico, vector<int> &jet_islep_nano_idx, vector<int> &jet_isvlep_nano_idx, vector<int> &sig_mu_pico_idx, bool isZgamma, bool isFastsim){
  vector<float> Jet_pt, Jet_mass;
  getJetWithJEC(nano, isFastsim, Jet_pt, Jet_mass);
  vector<int> Muon_fsrPhotonIdx;
  getMuon_fsrPhotonIdx(nano, nanoaod_version, Muon_fsrPhotonIdx);
  vector<int> Muon_nTrackerLayers;
  getMuon_nTrackerLayers(nano, nanoaod_version, Muon_nTrackerLayers);
  vector<int> Muon_genPartIdx;;
  if (!isData) getMuon_genPartIdx(nano, nanoaod_version, Muon_genPartIdx);

  //first, determine ordering based on signal and pt
  std::vector<NanoOrderEntry> nano_entries;
  for(int imu(0); imu<nano.nMuon(); ++imu){
    NanoOrderEntry nano_entry;
    nano_entry.nano_idx = imu;
    nano_entry.pt = nano.Muon_pt()[imu];
    nano_entry.is_sig = IsSignal(nano, imu, isZgamma);
    nano_entries.push_back(nano_entry);
  }
  std::sort(nano_entries.begin(),nano_entries.end(), 
      [](NanoOrderEntry a, NanoOrderEntry b) {
        if (a.is_sig && !b.is_sig) return true;
        if (b.is_sig && !a.is_sig) return false;
        return (a.pt>b.pt);
      });
  std::vector<int> ordered_nano_indices;
  for (NanoOrderEntry nano_entry : nano_entries)
    ordered_nano_indices.push_back(nano_entry.nano_idx);

  vector<int> sig_mu_nano_idx;
  pico.out_nmu() = 0; pico.out_nvmu() = 0;
  int pico_idx = 0;
  for(int imu : ordered_nano_indices){
    float pt = nano.Muon_pt()[imu];
    float eta = nano.Muon_eta()[imu];
    bool isSignal = false;
    if(isZgamma) { // For Zgamma productions
      if (pt <= ZgMuonPtCut) continue;
      if (fabs(eta) > MuonEtaCut) continue;
      if (fabs(nano.Muon_dz()[imu])>dzCut)  continue;
      if (fabs(nano.Muon_dxy()[imu])>dxyCut) continue; 
      isSignal = IsSignal(nano, imu, isZgamma);
      pico.out_mu_sip3d().push_back(nano.Muon_sip3d()[imu]);
      pico.out_mu_mediumid().push_back(nano.Muon_mediumId()[imu]);
      pico.out_mu_tightid().push_back(nano.Muon_tightId()[imu]);
      pico.out_mu_highptid().push_back(nano.Muon_highPtId()[imu]);
      pico.out_mu_fsrphotonid().push_back(Muon_fsrPhotonIdx[imu]);
    }
    else {
      if (!nano.Muon_mediumId()[imu]) continue;
      if (pt <= VetoMuonPtCut) continue;
      if (fabs(eta) > MuonEtaCut) continue;
      isSignal = IsSignal(nano, imu, isZgamma);
    }
    pico.out_mu_raw_pt().push_back(pt);
    pico.out_mu_raw_ptErr().push_back(nano.Muon_ptErr()[imu]);
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
    if (year < 2022) {
      //Rochester corrections, see https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/common/muonScaleResProducer.py
      if (isData) {
        pico.out_mu_pt().push_back(pt*rc.kScaleDT(nano.Muon_charge()[imu],pt,eta,nano.Muon_phi()[imu]));
        pico.out_mu_ptErr().push_back(pt*rc.kScaleDTerror(nano.Muon_charge()[imu],pt,eta,nano.Muon_phi()[imu]));
      }
      else {
        if (Muon_genPartIdx[imu] > 0 && Muon_genPartIdx[imu] < nano.nGenPart()) {
          float gen_pt = nano.GenPart_pt()[Muon_genPartIdx[imu]];
          pico.out_mu_pt().push_back(pt*rc.kSpreadMC(nano.Muon_charge()[imu],pt,eta,nano.Muon_phi()[imu],gen_pt));
          pico.out_mu_ptErr().push_back(pt*rc.kSpreadMCerror(nano.Muon_charge()[imu],pt,eta,nano.Muon_phi()[imu],gen_pt));
        }
        else {
          float unif_rand = rng.Uniform();
          pico.out_mu_pt().push_back(pt*rc.kSmearMC(nano.Muon_charge()[imu],pt,eta,nano.Muon_phi()[imu],Muon_nTrackerLayers[imu],unif_rand));
          pico.out_mu_ptErr().push_back(pt*rc.kSmearMCerror(nano.Muon_charge()[imu],pt,eta,nano.Muon_phi()[imu],Muon_nTrackerLayers[imu],unif_rand));
        }
      }
    }
    else if (nanoaod_version > 11.95) {
      pico.out_mu_pt().push_back(nano.Muon_bsConstrainedPt()[imu]);
      pico.out_mu_ptErr().push_back(nano.Muon_bsConstrainedPtErr()[imu]);
    }
    else {
      pico.out_mu_pt().push_back(pt);
      pico.out_mu_ptErr().push_back(nano.Muon_ptErr()[imu]);
    }
    if (!isData)
      pico.out_mu_pflavor().push_back(nano.Muon_genPartFlav()[imu]);

    // veto muon selection
    if (nano.Muon_miniPFRelIso_all()[imu] < MuonMiniIsoCut && 
        fabs(nano.Muon_dz()[imu])<=0.5f &&
        fabs(nano.Muon_dxy()[imu])<=0.2f) {
      pico.out_nvmu()++;
      pico.out_nvlep()++;
      // save indices of matching jets
      for (int ijet(0); ijet<nano.nJet(); ijet++) {
        if (dR(nano.Muon_eta()[imu], nano.Jet_eta()[ijet], nano.Muon_phi()[imu], nano.Jet_phi()[ijet])<0.4f &&
          fabs(Jet_pt[ijet] - nano.Muon_pt()[imu])/nano.Muon_pt()[imu] < 1.0f)
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
        if (dR(nano.Muon_eta()[imu], nano.Jet_eta()[ijet], nano.Muon_phi()[imu], nano.Jet_phi()[ijet])<0.4f &&
          fabs(Jet_pt[ijet] - nano.Muon_pt()[imu])/nano.Muon_pt()[imu] < 1.0f)
          jet_islep_nano_idx.push_back(ijet);
      }
    }
    pico_idx++;
  }
  return sig_mu_nano_idx;
}
