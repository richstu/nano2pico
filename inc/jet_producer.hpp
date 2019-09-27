#ifndef H_JET_PRODUCER
#define H_JET_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class JetProducer{
public:

  explicit JetProducer(int year);
  ~JetProducer();

  const float JetPtCut      = 30.0;
  const float JetEtaCut    = 2.4;

  void WriteJets(nano_tree &nano, pico_tree &pico);

private:
  int year;
  
};

#endif