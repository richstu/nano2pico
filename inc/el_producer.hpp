#ifndef H_EL_PRODUCER
#define H_EL_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class ElectronProducer{
public:

  explicit ElectronProducer(int year, bool isData);
  ~ElectronProducer();

  const float SignalElectronPtCut  = 20.0;
  const float VetoElectronPtCut    = 10.0;
  const float ZgElectronPtCut      =  7.0;
  const float ElectronEtaCut     = 2.5;
  const float ElectronMiniIsoCut = 0.1;
  const float ElectronRelIsoCut = 0.35;

  std::vector<int> WriteElectrons(nano_tree &nano, pico_tree &pico, std::vector<int> &jet_islep_nano_idx, bool isZgamma, bool isTTZ);

private:
  int year;
  bool isData;

  bool idElectron_noIso(int bitmap, int level);
  
};

#endif
