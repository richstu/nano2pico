#ifndef H_MU_PRODUCER
#define H_MU_PRODUCER

#include <string>

#include "TRandom3.h"

#include "nano_tree.hpp"
#include "pico_tree.hpp"
#include "RoccoR.hpp"

class MuonProducer{
public:

  explicit MuonProducer(int year, bool isData, float nanoaod_version, std::string rocco_file);
  ~MuonProducer();

  const float SignalMuonPtCut  = 20.0;
  const float VetoMuonPtCut    = 10.0;
  const float ZgMuonPtCut      =  5.0;
  const float MuonEtaCut        = 2.4;
  const float MuonMiniIsoCut    = 0.2;
  const float MuonRelIsoCut     = 0.35;
  const float dzCut         = 0.1;
  const float dxyCut        = 0.5;
  const float MuonHighPt        = 200;
  const float MuonSip3dCut      = 4.0;

  std::vector<int> WriteMuons(nano_tree &nano, pico_tree &pico, std::vector<int> &jet_islep_nano_idx, std::vector<int> &jet_isvlep_nano_idx, std::vector<int> &sig_mu_pico_idx, bool isZgamma, bool isFastsim);

private:
  int year;
  bool isData;
  RoccoR rc;
  TRandom3 rng;
  float nanoaod_version;

  bool IsSignal(nano_tree &nano, int nano_idx, bool isZgamma);
  
};

#endif
