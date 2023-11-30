#ifndef H_EL_PRODUCER
#define H_EL_PRODUCER

#include <string>
#include <memory>

#include "correction.hpp"
#include "nano_tree.hpp"
#include "pico_tree.hpp"

class ElectronProducer{
public:

  explicit ElectronProducer(int year, bool isData, bool preVFP);
  ~ElectronProducer();

  const float SignalElectronPtCut  = 20.0;
  const float VetoElectronPtCut    = 10.0;
  const float ZgElectronPtCut      =  7.0;
  const float ElectronEtaCut     = 2.5;
  const float ElectronMiniIsoCut = 0.1;
  const float ElectronRelIsoCut = 0.35;

  std::vector<int> WriteElectrons(nano_tree &nano, pico_tree &pico, 
                                  std::vector<int> &jet_islep_nano_idx, 
                                  std::vector<int> &jet_isvlep_nano_idx, 
                                  std::vector<int> &sig_el_pico_idx, 
                                  std::vector<int> &photon_el_pico_idx, 
                                  bool isZgamma, bool isFastsim);

private:
  int year;
  bool isData;
  std::unique_ptr<correction::CorrectionSet> cs_scale_syst_;
  correction::Correction::Ref map_scale_syst_;
  std::string str_scale_syst_;

  bool IsSignal(nano_tree& nano, int nano_idx, bool isZgamma);
  bool idElectron_noIso(int bitmap, int level);
  bool EcalDriven(int bitmap);
};

#endif
