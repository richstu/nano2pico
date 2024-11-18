#ifndef H_EVENT_WEIGHTER
#define H_EVENT_WEIGHTER

#include <string>
#include <vector>
#include <utility>

#include "TH2D.h"
#include "TH2F.h"

#include "pico_tree.hpp"

#include "correction.hpp"

class EventWeighter{
public:
  EventWeighter(std::string year, const std::vector<float> &btag_wpts);

  void ElectronSF(pico_tree &pico);

  void PhotonSF(pico_tree &pico);

  void MuonSF(pico_tree &pico);

  void PileupSF(pico_tree &pico);

  void bTaggingSF(pico_tree &pico);

private:
  std::string in_file_electron_;
  std::string in_file_photon_;
  std::string in_file_photon_mceff_;
  std::string in_file_muon_;
  std::string in_file_pu_;
  std::string in_file_btag_;
  std::string in_file_btag_mceff_;
  std::string key_;
  std::string puName_;
  std::string year_;
  std::unique_ptr<correction::CorrectionSet> cs_electron_;
  std::unique_ptr<correction::CorrectionSet> cs_photon_;
  std::unique_ptr<correction::CorrectionSet> cs_photon_mceff_;
  std::unique_ptr<correction::CorrectionSet> cs_muon_;
  std::unique_ptr<correction::CorrectionSet> cs_pileup_;
  std::unique_ptr<correction::CorrectionSet> cs_btag_;
  std::unique_ptr<correction::CorrectionSet> cs_btag_mceff_;
  correction::Correction::Ref map_photon_id_;
  correction::Correction::Ref map_photon_csev_;
  correction::Correction::Ref map_photon_mceff_;
  correction::Correction::Ref map_photon_mcunc_;
  correction::Correction::Ref map_electron_id_pass_;
  correction::Correction::Ref map_electron_id_pass_unc_;
  correction::Correction::Ref map_electron_id_fail_;
  correction::Correction::Ref map_electron_id_fail_unc_;
  correction::Correction::Ref map_muon_id_pass_;
  correction::Correction::Ref map_muon_id_pass_unc_;
  correction::Correction::Ref map_muon_id_fail_;
  correction::Correction::Ref map_muon_id_fail_unc_;
  correction::Correction::Ref map_pileup_;
  correction::Correction::Ref map_btag_;
  correction::Correction::Ref map_udsgtag_;
  float btag_wp_loose_;
  float btag_wp_medium_;
  float btag_wp_tight_;
};

#endif
