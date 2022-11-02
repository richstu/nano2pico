#ifndef H_PHOTON_WEIGHTER
#define H_PHOTON_WEIGHTER

#include <vector>
#include <utility>

#include "TH2D.h"
#include "TH2F.h"

#include "pico_tree.hpp"

class PhotonWeighter{
public:
  PhotonWeighter(int year, bool isZgamma);

  void FullSim(pico_tree &pico, float &w_photon, std::vector<float> &sys_photon);
  
private:
  std::string in_full_photon_id_, hist_full_photon_id_;
  std::string in_full_photon_ev_, hist_full_photon_ev_;

  bool do_full_photon_id_, do_full_photon_ev_;

  int year_;

  TH2F sf_full_photon_id_;
  TH2F sf_full_photon_ev_;

  float GetRegion(float const & eta, float const & r9);
  std::pair<double, double> GetPhotonScaleFactor(pico_tree &pico, std::size_t iph);
};

#endif
