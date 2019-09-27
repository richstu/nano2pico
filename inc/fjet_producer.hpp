#ifndef H_FATJET_PRODUCER
#define H_FATJET_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class FatJetProducer{
public:

  explicit FatJetProducer(int year);
  ~FatJetProducer();

  const float FatJetPtCut      = 200.0;
  const float FatJetEtaCut    = 2.4;

  void WriteFatJets(nano_tree &nano, pico_tree &pico);

private:
  int year;
  
};

#endif