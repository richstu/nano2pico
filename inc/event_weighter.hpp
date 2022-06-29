#ifndef H_EVENT_WEIGHTER
#define H_EVENT_WEIGHTER

#include <vector>
#include <utility>

#include "TH2D.h"
#include "TH2F.h"

#include "pico_tree.hpp"

class EventWeighter{
public:
  EventWeighter(int year, bool preVFP = false);

  void ElectronIDSF(pico_tree &pico, float &w_el_id);

  void PhotonIDSF(pico_tree &pico, float &w_photon_id);

  void PhotonCSEVSF(pico_tree &pico, float &w_photon_csev);

  void MuonIDSF(pico_tree &pico, float &w_muon_id);

  void MuonIsoSF(pico_tree &pico, float &w_muon_iso);

  void PileupSF(pico_tree &pico, float &w_pu);

private:
  std::string in_file_electron_;
  std::string in_file_photon_;
  std::string in_file_muon_;
  std::string in_file_pu_;
  std::string key_;
  std::string puName_;
};

#endif
