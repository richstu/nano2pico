#ifndef H_TK_PRODUCER
#define H_TK_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class IsoTrackProducer{
public:

  explicit IsoTrackProducer(int year);
  ~IsoTrackProducer();

  const float LeptonIsoTrackPtCut  = 10.0;
  const float LeptonIsoTrackEtaCut  = 2.4;
  const float HadronIsoTrackPtCut    = 5.0;
  const float HadronIsoTrackEtaCut    = 2.4;

  void WriteIsoTracks(nano_tree &nano, pico_tree &pico);

private:
  int year;
  
};

#endif