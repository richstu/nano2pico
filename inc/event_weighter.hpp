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

  void ElectronIDSF(pico_tree &pico, float &w_el_id, std::vector<float> &sys_lep);

  void PhotonIDSF(pico_tree &pico, float &w_photon_id);

  void PhotonCSEVSF(pico_tree &pico, float &w_photon_csev);

  void MuonTotalSF(pico_tree &pico, float &w_muon_tot, std::vector<float> &sys_lep);

  void PileupSF(pico_tree &pico, float &w_pu, float &sys_pu_up, float &sys_pu_down);

  void bTaggingSF(pico_tree &pico, float &w_btag);

private:
  std::string in_file_electron_;
  std::string in_file_photon_;
  std::string in_file_muon_;
  std::string in_file_muon_lowpt_reco_;
  std::string in_file_muon_lowpt_id_;
  std::string in_file_muon_mceff_;
  std::string in_file_pu_;
  std::string in_file_btag_;
  std::string key_;
  std::string puName_;
  std::unique_ptr<correction::CorrectionSet> cs_electron_;
  std::unique_ptr<correction::CorrectionSet> cs_photon_;
  std::unique_ptr<correction::CorrectionSet> cs_muon_lowpt_reco_;
  std::unique_ptr<correction::CorrectionSet> cs_muon_lowpt_id_;
  std::unique_ptr<correction::CorrectionSet> cs_muon_mceff_;
  std::unique_ptr<correction::CorrectionSet> cs_muon_;
  std::unique_ptr<correction::CorrectionSet> cs_pileup_;
  std::unique_ptr<correction::CorrectionSet> cs_btag_;
  correction::Correction::Ref map_electron_;
  correction::Correction::Ref map_photon_id_;
  correction::Correction::Ref map_photon_csev_;
  correction::Correction::Ref map_muon_looseid_;
  correction::Correction::Ref map_muon_highptid_;
  correction::Correction::Ref map_muon_iso_;
  correction::Correction::Ref map_muon_lowpt_reco_;
  correction::Correction::Ref map_muon_lowpt_id_;
  correction::Correction::Ref map_muon_mceff_;
  correction::Correction::Ref map_pileup_;
  correction::Correction::Ref map_btag_;
};

#endif
