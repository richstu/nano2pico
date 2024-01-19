#ifndef H_EVENT_WEIGHTER
#define H_EVENT_WEIGHTER

#include <vector>
#include <utility>

#include "TH2D.h"
#include "TH2F.h"

#include "pico_tree.hpp"

#include "correction.hpp"

class EventWeighter{
public:
  EventWeighter(int year, bool preVFP, const std::vector<float> &btag_wpts);

  void ElectronSF(pico_tree &pico);

  void PhotonIDSF(pico_tree &pico, float &w_photon_id);

  void PhotonCSEVSF(pico_tree &pico, float &w_photon_csev, std::vector<float> &sys_photon_csev);

  void MuonSF(pico_tree &pico);

  void PileupSF(pico_tree &pico);

  void bTaggingSF(pico_tree &pico);

private:
  std::string in_file_electron_;
  std::string in_file_photon_;
  std::string in_file_photon_mceff_;
  std::string in_file_muon_;
  std::string in_file_muon_lowpt_reco_;
  std::string in_file_muon_lowpt_id_;
  std::string in_file_muon_mceff_;
  std::string in_file_pu_;
  std::string in_file_btag_;
  std::string in_file_btag_mceff_;
  std::string key_;
  std::string puName_;
  std::unique_ptr<correction::CorrectionSet> cs_electron_;
  std::unique_ptr<correction::CorrectionSet> cs_photon_;
  std::unique_ptr<correction::CorrectionSet> cs_photon_mceff_;
  std::unique_ptr<correction::CorrectionSet> cs_muon_lowpt_reco_;
  std::unique_ptr<correction::CorrectionSet> cs_muon_lowpt_id_;
  std::unique_ptr<correction::CorrectionSet> cs_muon_mceff_;
  std::unique_ptr<correction::CorrectionSet> cs_muon_;
  std::unique_ptr<correction::CorrectionSet> cs_pileup_;
  std::unique_ptr<correction::CorrectionSet> cs_btag_;
  std::unique_ptr<correction::CorrectionSet> cs_btag_mceff_;
  correction::Correction::Ref map_electron_;
  correction::Correction::Ref map_photon_id_;
  correction::Correction::Ref map_photon_csev_;
  correction::Correction::Ref map_photon_csev_mceff_;
  correction::Correction::Ref map_muon_looseid_;
  correction::Correction::Ref map_muon_highptid_;
  correction::Correction::Ref map_muon_iso_;
  correction::Correction::Ref map_muon_lowpt_reco_;
  correction::Correction::Ref map_muon_lowpt_id_;
  correction::Correction::Ref map_muon_mceff_;
  correction::Correction::Ref map_pileup_;
  correction::Correction::Ref map_btag_;
  correction::Correction::Ref map_udsgtag_;
  float btag_wp_loose_;
  float btag_wp_medium_;
  float btag_wp_tight_;
};

#endif
