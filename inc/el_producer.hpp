#ifndef H_EL_PRODUCER
#define H_EL_PRODUCER

#include <string>
#include <memory>

#include "TRandom3.h"

#include "correction.hpp"
#include "nano_tree.hpp"
#include "pico_tree.hpp"

class ElectronProducer{
public:

  explicit ElectronProducer(std::string year, bool isData, float nanoaod_version);
  ~ElectronProducer();

  const float SignalElectronPtCut  = 20.0;
  const float VetoElectronPtCut    = 10.0;
  const float ZgElectronPtCut      =  7.0;
  const float ElectronEtaCut     = 2.5;
  const float ElectronMiniIsoCut = 0.1;
  const float ElectronRelIsoCut = 0.35;

  const float dxyCut = 0.5;
  const float dzCut = 1.0;
  const float jetDRCut = 0.4;
  const float jetpTCut = 1.0;

  std::vector<int> WriteElectrons(nano_tree &nano, pico_tree &pico, 
                                  std::vector<int> &jet_islep_nano_idx, 
                                  std::vector<int> &jet_isvlep_nano_idx, 
                                  std::vector<int> &sig_el_pico_idx, 
                                  std::vector<int> &photon_el_pico_idx, 
                                  bool isZgamma, bool isFastsim);

  float ConvertMVA(float mva_mini);

private:
  std::string year;
  bool isData;
  std::unique_ptr<correction::CorrectionSet> cs_scale_syst_;
  correction::Correction::Ref map_scale_syst_; //run2, just has uncertainties
  correction::Correction::Ref map_scale_; //run3 has scale correction and uncertainty
  correction::Correction::Ref map_smearing_; //run3 has smearing correction and uncertainty
  std::string str_scale_syst_;
  TRandom3 rng_;
  float nanoaod_version;

  bool IsSignal(nano_tree& nano, int nano_idx, bool isZgamma);
  bool idElectron_noIso(int bitmap, int level);
  bool EcalDriven(int bitmap);
};

#endif
