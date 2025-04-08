#ifndef H_PHOTON_PRODUCER
#define H_PHOTON_PRODUCER

#include <string>
#include <vector>

#include "TRandom3.h"

#include "correction.hpp"
#include "nano_tree.hpp"
#include "pico_tree.hpp"

class PhotonProducer{
public:

  explicit PhotonProducer(std::string year, bool isData, float nanoaod_version);
  ~PhotonProducer();

  // check what these should be in a relevant AN
  const float PhotonPtCut       = 10.0;
  const float SignalPhotonPtCut = 15.0;
  const float PhotonEtaCut      = 2.5;
  const float PhotonRelIsoCut   = 0.1; 

  const float FsrPhotonPtCut    = 2.0;
  const float FsrPhotonEtaCut   = 2.4;
  const float FsrPhotonIsoCut   = 1.8;
  const float FsrPhotondRCut    = 0.012;
  const float FsrSeparationReq  = 0.2;

  std::vector<int> WritePhotons(nano_tree &nano, pico_tree &pico, 
                                std::vector<int> &jet_isphoton_nano_idx, 
                                std::vector<int> &sig_el_nano_idx, 
                                std::vector<int> &sig_mu_nano_idx,
                                std::vector<int> &photon_el_pico_idx);

private:
  std::string year;
  bool isData;
  float nanoaod_version;
  TRandom3 rng_;

  bool idPhoton(int bitmap, int level);
  std::unique_ptr<correction::CorrectionSet> cs_scale_syst_;
  correction::Correction::Ref map_scale_syst_;
  correction::CompoundCorrection::Ref map_scale_;
  correction::Correction::Ref map_smearing_;
  std::string str_scale_syst_;
};

#endif
