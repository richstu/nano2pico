#ifndef H_JET_PRODUCER
#define H_JET_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class JetProducer{
public:

  explicit JetProducer(int year);
  ~JetProducer();

  const float JetPtCut     = 30.0;
  const float JetEtaCut    = 2.4;

  std::vector<int> WriteJets(nano_tree &nano, pico_tree &pico, std::vector<int> jet_islep_nano_idx,
                 const std::vector<float> &btag_wpts, const std::vector<float> &btag_df_wpts);
  void WriteJetSys(nano_tree &nano, pico_tree &pico, 
                              std::vector<int> &sig_jet_nano_idx, const float &btag_wpt);
private:
  int year;
};

#endif