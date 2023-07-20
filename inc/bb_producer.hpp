#ifndef H_BB_PRODUCER
#define H_BB_PRODUCER

#include "pico_tree.hpp"

class BBVarProducer{
public:

  explicit BBVarProducer(int year);
  ~BBVarProducer();

  void WriteBBVars(pico_tree &pico, bool doDeelFlav);

private:
  int year;
  
};

#endif
