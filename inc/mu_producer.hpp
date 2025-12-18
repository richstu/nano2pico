#ifndef H_MU_PRODUCER
#define H_MU_PRODUCER

#include <memory>
#include <string>

#include "TRandom3.h"

#include "correction.hpp"
#include "nano_tree.hpp"
#include "pico_tree.hpp"
#include "RoccoR.hpp"

class MuonProducer{
public:

  explicit MuonProducer(std::string year, bool isData, float nanoaod_version, std::string rocco_file);
  ~MuonProducer();

  const float SignalMuonPtCut  = 20.0;
  const float VetoMuonPtCut    = 10.0;
  const float ZgMuonPtCut      =  5.0;
  const float PicoMuonPtCut    =  3.0;
  const float MuonEtaCut        = 2.4;
  const float MuonMiniIsoCut    = 0.2;
  const float MuonRelIsoCut     = 0.35;
  const float dzCut         = 1.0;
  const float dxyCut        = 0.5;
  const float MuonHighPt        = 200;
  const float MuonSip3dCut      = 4.0;

  std::vector<int> WriteMuons(nano_tree &nano, pico_tree &pico, std::vector<int> &jet_islep_nano_idx, std::vector<int> &jet_isvlep_nano_idx, std::vector<int> &sig_mu_pico_idx, bool isZgamma, bool is_signal_sample, bool isFastsim);

private:
  std::string year;
  bool isData;
  RoccoR rc;
  TRandom3 rng;
  float nanoaod_version;
  bool run3;
  std::unique_ptr<correction::CorrectionSet> cs_scare_;

  bool IsSignal(nano_tree &nano, int nano_idx, bool isZgamma, float pt, 
                bool skip_pt=false);
  
};

#endif
