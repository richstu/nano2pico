#ifndef H_JET_PRODUCER
#define H_JET_PRODUCER

#include <string>

#include "correction.hpp"
#include "hig_producer.hpp"
#include "met_producer.hpp"
#include "nano_tree.hpp"
#include "pico_tree.hpp"
#include "correction.hpp"

#include "TLorentzVector.h"
#include "TRandom3.h"

class JetMetProducer{
public:

  explicit JetMetProducer(int year, std::string year_string, 
                          float nanoaod_version, float min_jet_pt, 
                          float max_jet_eta, bool isData, bool is_preUL, 
                          bool verbose=false);
  ~JetMetProducer();

  void SetVerbose(bool verbose_){ verbose = verbose_; };

  void WriteMet(nano_tree &nano, pico_tree &pico);
  std::vector<int> WriteJetMet(nano_tree &nano, pico_tree &pico, 
                               std::vector<int> jet_islep_nano_idx, 
                               std::vector<int> jet_isvlep_nano_idx, 
                               std::vector<int> jet_isphoton_nano_idx,
                               const std::vector<float> &btag_wpts, 
                               const std::vector<float> &btag_df_wpts, 
                               bool isFastsim, 
                               bool isSignal,
                               bool is2022preEE,
                               std::vector<HiggsConstructionVariables> &sys_higvars);
  void WriteFatJets(nano_tree &nano, pico_tree &pico);
  void WriteSubJets(nano_tree &nano, pico_tree &pico);
  void WriteJetSystemPt(nano_tree &nano, pico_tree &pico, 
                        std::vector<int> &sig_jet_nano_idx, 
                        const float &btag_wpt, bool isFastsim);
private:

  void GetJetUncertainties(nano_tree &nano, pico_tree &pico, 
                           std::vector<float> &jer_nm_factor, 
                           std::vector<float> &jer_up_factor,
                           std::vector<float> &jer_dn_factor,
                           std::vector<float> &jes_up_factor,
                           std::vector<float> &jes_dn_factor);

  int year;
  std::string year_string;
  float nanoaod_version;
  bool verbose;
  float min_jet_pt;
  float max_jet_eta;
  bool isData;
  bool is_preUL;
  MetProducer met_producer;
  TRandom3 rng_;
  std::unique_ptr<correction::CorrectionSet> cs_jerc_;
  correction::Correction::Ref map_jes_;
  correction::Correction::Ref map_jersf_;
  correction::Correction::Ref map_jermc_;
  correction::CompoundCorrection::Ref map_jec_;
  std::string in_file_jetveto_;
  std::unique_ptr<correction::CorrectionSet> cs_jetveto_;
  correction::Correction::Ref map_jetveto_;
};

#endif
