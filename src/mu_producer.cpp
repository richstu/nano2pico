#include "mu_producer.hpp"


#include "correction.hpp"
#include "muon_scare.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <memory>

using namespace std;

//TODO update year to string in process_nano
MuonProducer::MuonProducer(string year_, bool isData_, float nanoaod_version_, std::string rocco_file) :
  year(year_),
  isData(isData_),
  rc(rocco_file),
  rng(4357),
  nanoaod_version(nanoaod_version_),
  run3(false){
  if (year=="2022") {
    cs_scare_ = correction::CorrectionSet::from_file(
        "data/zgamma/2022/2022_Summer22.json");
    run3 = true;
  }
  else if (year=="2022EE") {
    cs_scare_ = correction::CorrectionSet::from_file(
        "data/zgamma/2022EE/2022_Summer22EE.json");
    run3 = true;
  }
  else if (year=="2023") {
    cs_scare_ = correction::CorrectionSet::from_file(
        "data/zgamma/2023/2023_Summer23.json");
    run3 = true;
  }
  else if (year=="2023BPix") {
    cs_scare_ = correction::CorrectionSet::from_file(
        "data/zgamma/2023BPix/2023_Summer23BPix.json");
    run3 = true;
  }
}

MuonProducer::~MuonProducer(){
}

bool MuonProducer::IsSignal(nano_tree &nano, int nano_idx, bool isZgamma, float pt, bool skip_pt) {
  float eta = nano.Muon_eta()[nano_idx];
  if(isZgamma) { // For Zgamma productions
    if (pt <= ZgMuonPtCut && !skip_pt) return false;
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
    if (pt <= VetoMuonPtCut && !skip_pt) return false;
    if (fabs(eta) > MuonEtaCut) return false;
    if ((pt > SignalMuonPtCut || skip_pt) &&
      nano.Muon_miniPFRelIso_all()[nano_idx] < MuonMiniIsoCut &&
      fabs(nano.Muon_dz()[nano_idx])<=0.5f && 
      fabs(nano.Muon_dxy()[nano_idx])<=0.2f)
      return true;
    return false;
  }
  return false;
}

vector<int> MuonProducer::WriteMuons(nano_tree &nano, pico_tree &pico, vector<int> &jet_islep_nano_idx, vector<int> &jet_isvlep_nano_idx, vector<int> &sig_mu_pico_idx, bool isZgamma, bool isSignal_sample, bool isFastsim){
  vector<float> Jet_pt, Jet_mass;
  getJetWithJEC(nano, isFastsim, Jet_pt, Jet_mass);
  vector<int> Muon_fsrPhotonIdx;
  getMuon_fsrPhotonIdx(nano, nanoaod_version, Muon_fsrPhotonIdx);
  vector<int> Muon_nTrackerLayers;
  getMuon_nTrackerLayers(nano, nanoaod_version, Muon_nTrackerLayers);
  vector<int> Muon_genPartIdx;
  if (!isData) getMuon_genPartIdx(nano, nanoaod_version, Muon_genPartIdx);

  //scale+resolution corrections
  vector<float> muon_pt_corr;
  vector<float> muon_pt_scaleup;
  vector<float> muon_pt_scaledn;
  vector<float> muon_pt_resup;
  vector<float> muon_pt_resdn;
  for(int imu(0); imu<nano.nMuon(); ++imu){
    float eta = nano.Muon_eta()[imu];
    float phi = nano.Muon_phi()[imu];
    int charge = nano.Muon_charge()[imu];
    int nTrackerLayers = Muon_nTrackerLayers[imu];
    if (!run3) {
      //Rochester corrections, see https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/common/muonScaleResProducer.py
      float pt = nano.Muon_pt()[imu];
      float scale_sf = rc.kScaleDT(charge,pt,eta,phi);
      if (isData) {
        muon_pt_corr.push_back(pt*scale_sf);
      }
      else {
        float resolution_sf = 1.0;
        float scale_unc = rc.kScaleDTerror(charge,pt,eta,phi);
        float resolution_unc = 1.0;
        if (Muon_genPartIdx[imu] > 0 
            && Muon_genPartIdx[imu] < nano.nGenPart()) {
          float gen_pt = nano.GenPart_pt()[Muon_genPartIdx[imu]];
          resolution_sf = rc.kSpreadMC(charge,pt,eta,phi,gen_pt);
          resolution_unc = rc.kSpreadMCerror(charge,pt,eta,phi,gen_pt)
            *(resolution_sf > 1.0 ? 1.0 : -1.0);
        }
        else {
          float unif_rand = rng.Uniform();
          resolution_sf = rc.kSmearMC(charge,pt,eta,phi,nTrackerLayers,
              unif_rand);
          resolution_unc = rc.kSmearMCerror(charge,pt,eta,phi,nTrackerLayers,
              unif_rand)*(resolution_sf > 1.0 ? 1.0 : -1.0);
        }
        muon_pt_corr.push_back(pt*resolution_sf);
        muon_pt_scaleup.push_back(pt*resolution_sf
                                  *(scale_sf+scale_unc)/scale_sf);
        muon_pt_scaledn.push_back(pt*resolution_sf
                                  *(scale_sf-scale_unc)/scale_sf);
        muon_pt_resup.push_back(pt*(resolution_sf+resolution_unc));
        muon_pt_resdn.push_back(pt*(resolution_sf-resolution_unc));
      }
    }
    else {
      float pt = nano.Muon_bsConstrainedPt()[imu];
      if (isData) {
        muon_pt_corr.push_back(scarekit::pt_scale(1, pt, eta, phi,
            charge, cs_scare_));
      }
      else {
        float sca_pt = scarekit::pt_scale(0, pt, eta, phi, charge, cs_scare_);
        float re_pt = scarekit::pt_resol(sca_pt, eta, 
            static_cast<float>(nTrackerLayers), cs_scare_);
        muon_pt_corr.push_back(re_pt);
        muon_pt_scaleup.push_back(scarekit::pt_scale_var(re_pt, eta, phi, 
            charge, "up", cs_scare_));
        muon_pt_scaledn.push_back(scarekit::pt_scale_var(re_pt, eta, phi, 
            charge, "dn", cs_scare_));
        muon_pt_resup.push_back(scarekit::pt_resol_var(sca_pt, re_pt, eta, 
            "up", cs_scare_));
        muon_pt_resdn.push_back(scarekit::pt_resol_var(sca_pt, re_pt, eta, 
            "dn", cs_scare_));
      }
    }
  }

  //first, determine ordering based on signal and pt
  std::vector<NanoOrderEntry> nano_entries;
  for(int imu(0); imu<nano.nMuon(); ++imu){
    NanoOrderEntry nano_entry;
    nano_entry.nano_idx = imu;
    nano_entry.pt = muon_pt_corr[imu];
    nano_entry.is_sig = IsSignal(nano, imu, isZgamma, muon_pt_corr[imu]);
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
    float pt = muon_pt_corr[imu];
    float eta = nano.Muon_eta()[imu];
    bool isSignal = false;
    bool isSignal_nopt = false;
    if(isZgamma) { // For Zgamma productions
      if (pt <= PicoMuonPtCut) continue;
      if (fabs(eta) > MuonEtaCut) continue;
      if (fabs(nano.Muon_dz()[imu])>dzCut)  continue;
      if (fabs(nano.Muon_dxy()[imu])>dxyCut) continue; 
      isSignal = IsSignal(nano, imu, isZgamma, pt);
      isSignal_nopt = IsSignal(nano, imu, isZgamma, pt, true);
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
      isSignal = IsSignal(nano, imu, isZgamma, pt);
    }
    pico.out_mu_raw_pt().push_back(nano.Muon_pt()[imu]);
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
    if (!isData && isSignal_sample) {
      pico.out_mu_pflavor().push_back(nano.Muon_genPartFlav()[imu]);
      pico.out_sys_mu_pt_scaleup().push_back(muon_pt_scaleup[imu]);
      pico.out_sys_mu_pt_scaledn().push_back(muon_pt_scaledn[imu]);
      pico.out_sys_mu_pt_resup().push_back(muon_pt_resup[imu]);
      pico.out_sys_mu_pt_resdn().push_back(muon_pt_resdn[imu]);
      pico.out_sys_mu_sig_scaleup().push_back(isSignal_nopt 
          && (muon_pt_scaleup[imu] > ZgMuonPtCut));
      pico.out_sys_mu_sig_scaledn().push_back(isSignal_nopt 
          && (muon_pt_scaledn[imu] > ZgMuonPtCut));
      pico.out_sys_mu_sig_resup().push_back(isSignal_nopt 
          && (muon_pt_resup[imu] > ZgMuonPtCut));
      pico.out_sys_mu_sig_resdn().push_back(isSignal_nopt 
          && (muon_pt_resdn[imu] > ZgMuonPtCut));
    }
    if (nanoaod_version > 11.95) {
      pico.out_mu_ptErr().push_back(nano.Muon_bsConstrainedPtErr()[imu]);
    }
    else {
      pico.out_mu_ptErr().push_back(nano.Muon_ptErr()[imu]);
    }

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
