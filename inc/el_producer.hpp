#ifndef H_EL_PRODUCER
#define H_EL_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class ElectronProducer{
public:

  explicit ElectronProducer(int year);
  ~ElectronProducer();

  const float SignalElectronPtCut  = 20.0;
  const float VetoElectronPtCut    = 10.0;
  const float ElectronEtaCut     = 2.5;
  const float ElectronMiniIsoCut = 0.1;

  void WriteElectrons(nano_tree &nano, pico_tree &pico);

private:
  int year;
  
};

#endif