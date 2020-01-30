#ifndef H_PHOTON_PRODUCER
#define H_PHOTON_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class PhotonProducer{
public:

  explicit PhotonProducer(int year);
  ~PhotonProducer();

  // check what these should be in a relevant AN
  const float PhotonPtCut       = 10.0;
  const float SignalPhotonPtCut = 15.0;
  const float PhotonEtaCut      = 2.5;
  const float PhotonRelIsoCut   = 0.1; 

 std::vector<int> WritePhotons(nano_tree &nano, pico_tree &pico, std::vector<int> &jet_isphoton_nano_idx, 
                                               std::vector<int> &sig_el_nano_idx, std::vector<int> &sig_mu_nano_idx);

private:
  int year;
  
};

#endif
