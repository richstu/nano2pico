#ifndef H_LEPTON_WEIGHTER
#define H_LEPTON_WEIGHTER

#include <vector>
#include <utility>

#include "TH2D.h"
#include "TH2F.h"

#include "pico_tree.hpp"

class LeptonWeighter{
public:
  LeptonWeighter(int year, bool isZgamma, bool doGH16 = false);

  void FullSim(pico_tree &pico, float &w_lep, std::vector<float> &sys_lep);
  void FastSim(pico_tree &pico, float &w_fs_lep, std::vector<float> &sys_fs_lep);
  
private:
  std::string in_full_mu_med_, hist_full_mu_med_;
  std::string in_full_mu_iso_, hist_full_mu_iso_;
  std::string in_full_mu_trk_, hist_full_mu_trk_;
  std::string in_full_el_med_, hist_full_el_med_;
  std::string in_full_el_iso_, hist_full_el_iso_;
  std::string in_full_el_trk_, hist_full_el_trk_;
  std::string in_fast_mu_med_, hist_fast_mu_med_;
  std::string in_fast_mu_iso_, hist_fast_mu_iso_;
  std::string in_fast_el_med_, hist_fast_el_med_;
  std::string in_fast_el_iso_, hist_fast_el_iso_;
  std::string in_full_el_id_,  hist_full_el_id_;
  std::string in_full_el_tk_,  hist_full_el_tk_;
  std::string in_full_el_lo_,  hist_full_el_lo_;
  std::string in_full_mu_id_,  hist_full_mu_id_;
  std::string in_full_mu_is_,  hist_full_mu_is_;

  bool do_full_el_med_, do_full_el_iso_, do_full_el_trk_;
  bool do_full_mu_med_, do_full_mu_iso_, do_full_mu_trk_;

  bool do_full_el_id_, do_full_el_tk_, do_full_el_lo_;
  bool do_full_mu_id_, do_full_mu_is_;

  bool do_fast_el_med_, do_fast_el_iso_;
  bool do_fast_mu_med_, do_fast_mu_iso_;

  TH2F sf_full_mu_med_;
  TH2F sf_full_mu_iso_;
  TH2F sf_full_mu_trk_;

  TH2F sf_full_el_med_;
  TH2F sf_full_el_iso_;
  TH2F sf_full_el_trk_;

  TH2F sf_full_el_id_;
  TH2F sf_full_el_tk_;
  TH2F sf_full_el_lo_;

  TH2F sf_full_mu_id_;
  TH2F sf_full_mu_is_;

  TH2D sf_fast_mu_med_;
  TH2D sf_fast_mu_iso_;
  
  TH2D sf_fast_el_med_;
  TH2D sf_fast_el_iso_;

  std::pair<double, double> GetMuonScaleFactor(pico_tree &pico, std::size_t ilep);
  std::pair<double, double> GetElectronScaleFactor(pico_tree &pico, std::size_t ilep);

  std::pair<double, double> GetMuonScaleFactorFS(pico_tree &pico, std::size_t ilep);
  std::pair<double, double> GetElectronScaleFactorFS(pico_tree &pico, std::size_t ilep);
};

#endif
