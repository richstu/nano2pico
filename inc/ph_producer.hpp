#ifndef H_PH_PRODUCER
#define H_PH_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class PhotonProducer{
public:

  explicit PhotonProducer(int year);
  ~PhotonProducer();

  // check what these should be in a relevant AN
  const float SignalPhotonPtCut  = 20.0;
  const float PhotonEtaCut     = 2.5;
  const float PhotonMiniIsoCut = 0.1;

  void WritePhotons(nano_tree &nano, pico_tree &pico);

private:
  int year;
  
};

#endif