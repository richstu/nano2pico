#ifndef H_TK_PRODUCER
#define H_TK_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class IsoTrackProducer{
public:

  explicit IsoTrackProducer(int year);
  ~IsoTrackProducer();

  // track cuts a bit involved so see directly IsGoodTk for readability

  void WriteIsoTracks(nano_tree &nano, pico_tree &pico, 
                      std::vector<int> &sig_el_nano_idx,
                      std::vector<int> &sig_mu_nano_idx, bool isFastsim, bool isUL);

  bool IsGoodTk(pico_tree &pico, bool isNanoElectron, bool isNanoMuon, int pdgid, float pt, float eta, float phi, 
                float miniso, float reliso, float d0, float dz, float mt);

private:
  int year;
  
};

#endif
