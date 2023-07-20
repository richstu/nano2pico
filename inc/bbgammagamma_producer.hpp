#ifndef H_BBGAMMAGAMMA_PRODUCER
#define H_BBGAMMAGAMMA_PRODUCER

#include "pico_tree.hpp"

class BBGammaGammaVarProducer{
public:

  explicit BBGammaGammaVarProducer(int year);
  ~BBGammaGammaVarProducer();

  void WriteBBGammaGammaVars(pico_tree &pico);

private:
  int year;
  
};

#endif
