#ifndef H_GAMMAGAMMA_PRODUCER
#define H_GAMMAGAMMA_PRODUCER

#include "pico_tree.hpp"

class GammaGammaVarProducer{
public:

  explicit GammaGammaVarProducer(int year);
  ~GammaGammaVarProducer();

  void WriteGammaGammaVars(pico_tree &pico);

private:
  int year;
  
};

#endif
