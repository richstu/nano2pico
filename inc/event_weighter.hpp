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
  EventWeighter(int year, bool preVFP = false);

  void ElectronIDSF(pico_tree &pico, float &w_el_id);

  void PhotonIDSF(pico_tree &pico, float &w_photon_id);

  void PhotonCSEVSF(pico_tree &pico, float &w_photon_csev);

  void MuonIDSF(pico_tree &pico, float &w_muon_id);

  void MuonIsoSF(pico_tree &pico, float &w_muon_iso);

  void PileupSF(pico_tree &pico, float &w_pu);

  void bTaggingSF(pico_tree &pico, float &w_btag);

private:
  std::string in_file_electron_;
  std::string in_file_photon_;
  std::string in_file_muon_;
  std::string in_file_pu_;
  std::string in_file_btag_;
  std::string key_;
  std::string puName_;
  std::unique_ptr<correction::CorrectionSet> cs_electron_;
  std::unique_ptr<correction::CorrectionSet> cs_photon_;
  std::unique_ptr<correction::CorrectionSet> cs_muon_;
  std::unique_ptr<correction::CorrectionSet> cs_pileup_;
  std::unique_ptr<correction::CorrectionSet> cs_btag_;
  correction::Correction::Ref map_electron_;
  correction::Correction::Ref map_photon_id_;
  correction::Correction::Ref map_photon_csev_;
  correction::Correction::Ref map_muon_looseid_;
  correction::Correction::Ref map_muon_highptid_;
  correction::Correction::Ref map_muon_iso_;
  correction::Correction::Ref map_pileup_;
  correction::Correction::Ref map_btag_;
};

#endif
