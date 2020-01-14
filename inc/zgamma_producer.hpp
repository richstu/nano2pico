#ifndef H_ZGAMMA_PRODUCER
#define H_ZGAMMA_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class ZGammaVarProducer{
public:

  explicit ZGammaVarProducer(int year);
  ~ZGammaVarProducer();

  void WriteZGammaVars(nano_tree &nano, pico_tree &pico, std::vector<int> sig_jet_nano_idx);

private:
  int year;
  
};

#endif
