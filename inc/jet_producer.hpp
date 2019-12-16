#ifndef H_JET_PRODUCER
#define H_JET_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class JetProducer{
public:

  explicit JetProducer(int year, bool verbose=false);
  ~JetProducer();

  const float JetPtCut     = 30.0;
  const float JetEtaCut    =  2.4;
  const float ZgJetEtaCut  =  4.7;

  void SetVerbose(bool verbose_){ verbose = verbose_; };

  std::vector<int> WriteJets(nano_tree &nano, pico_tree &pico, 
                             std::vector<int> jet_islep_nano_idx, std::vector<int> jet_isphoton_nano_idx,
                             const std::vector<float> &btag_wpts, const std::vector<float> &btag_df_wpts,
                             bool isZgamma);
  void WriteFatJets(nano_tree &nano, pico_tree &pico);
  void WriteSubJets(nano_tree &nano, pico_tree &pico);
  void WriteJetSystemPt(nano_tree &nano, pico_tree &pico, 
                              std::vector<int> &sig_jet_nano_idx, const float &btag_wpt);
private:
  int year;
  bool verbose;
};

#endif
