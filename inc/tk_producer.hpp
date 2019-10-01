#ifndef H_TK_PRODUCER
#define H_TK_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class IsoTrackProducer{
public:

  explicit IsoTrackProducer(int year);
  ~IsoTrackProducer();

  const float IsoTrackEtaCut  = 2.5;
  const float IsoTrackDzCut  = 0.1;
  const float IsoTrackMtCut  = 0.1;
  //dependent on track type
  const float LeptonIsoTrackPtCut  = 10.0;
  const float HadronIsoTrackPtCut    = 5.0;
  const float LeptonIsoTrackRelIsoCut  = 0.2;
  const float HadronIsoTrackRelIsoCut    = 0.1;

  std::vector<int> WriteIsoTracks(nano_tree &nano, pico_tree &pico);

private:
  int year;
  
};

#endif