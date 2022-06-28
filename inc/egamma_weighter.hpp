#ifndef H_EGAMMA_WEIGHTER
#define H_EGAMMA_WEIGHTER

#include <vector>
#include <utility>

#include "TH2D.h"
#include "TH2F.h"

#include "pico_tree.hpp"

class EgammaWeighter{
public:
  EgammaWeighter(int year, bool preVFP = false);

  void GetSF(pico_tree &pico, float &w_lep);

  void IDSF(pico_tree &pico, float &w_photon_id);

  void CSEVSF(pico_tree &pico, float &w_photon_csev);

private:
  std::string in_file_, in_file_photon_, key_;
};

#endif
