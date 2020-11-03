#ifndef H_PREFIRE_WEIGHTER
#define H_PREFIRE_WEIGHTER

#include <vector>
#include <utility>

#include "TH2D.h"
#include "TH2F.h"

#include "nano_tree.hpp"

class PrefireWeighter{
public:
  PrefireWeighter(int year, bool use_jet_empt=true);

  void EventWeight(nano_tree & nano, float & w_prefire, std::vector<float> & sys_prefire);
  
private:
  bool use_jet_empt_;

  std::string prefire_jet_empt_filename_, prefire_jet_empt_histname_;
  std::string prefire_jet_filename_, prefire_jet_histname_;
  std::string prefire_photon_filename_, prefire_photon_histname_;

  TH2F sf_hist_prefire_jet_, sf_hist_prefire_jet_empt_, sf_hist_prefire_photon_;

  bool do_prefire_jet_, do_prefire_jet_empt_, do_prefire_photon_;
};

#endif
