#ifndef H_JET_PRODUCER
#define H_JET_PRODUCER

#include "hig_producer.hpp"
#include "nano_tree.hpp"
#include "pico_tree.hpp"

class JetProducer{
public:

  explicit JetProducer(int year, float min_jet_pt, float max_jet_eta, bool isData, bool verbose=false);
  ~JetProducer();

  void SetVerbose(bool verbose_){ verbose = verbose_; };

  std::vector<int> WriteJets(nano_tree &nano, pico_tree &pico, 
                             std::vector<int> jet_islep_nano_idx, 
                             std::vector<int> jet_isvlep_nano_idx, 
                             std::vector<int> jet_isphoton_nano_idx,
                             const std::vector<float> &btag_wpts, 
                             const std::vector<float> &btag_df_wpts, 
                             bool isFastsim,
                             std::vector<HiggsConstructionVariables> &sys_higvars);
  void WriteFatJets(nano_tree &nano, pico_tree &pico);
  void WriteSubJets(nano_tree &nano, pico_tree &pico);
  void WriteJetSystemPt(nano_tree &nano, pico_tree &pico, 
                              std::vector<int> &sig_jet_nano_idx, const float &btag_wpt, bool isFastsim);
private:
  int year;
  bool verbose;
  float min_jet_pt;
  float max_jet_eta;
  bool isData;
};

#endif
