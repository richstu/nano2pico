#ifndef H_JET_PRODUCER
#define H_JET_PRODUCER

#include "correction.hpp"
#include "hig_producer.hpp"
#include "nano_tree.hpp"
#include "pico_tree.hpp"

#include "TLorentzVector.h"
#include "TRandom3.h"

class JetProducer{
public:

  explicit JetProducer(int year, float nanoaod_version, float min_jet_pt, 
                       float max_jet_eta, bool isData, bool preVFP, bool verbose=false);
  ~JetProducer();

  void SetVerbose(bool verbose_){ verbose = verbose_; };

  std::vector<int> WriteJets(nano_tree &nano, pico_tree &pico, 
                             std::vector<int> jet_islep_nano_idx, 
                             std::vector<int> jet_isvlep_nano_idx, 
                             std::vector<int> jet_isphoton_nano_idx,
                             const std::vector<float> &btag_wpts, 
                             const std::vector<float> &btag_df_wpts, 
                             bool isFastsim, 
                             bool isSignal,
                             bool is_preUL,
                             std::vector<HiggsConstructionVariables> &sys_higvars);
  void WriteFatJets(nano_tree &nano, pico_tree &pico);
  void WriteSubJets(nano_tree &nano, pico_tree &pico);
  void WriteJetSystemPt(nano_tree &nano, pico_tree &pico, 
                        std::vector<int> &sig_jet_nano_idx, 
                        const float &btag_wpt, bool isFastsim);
private:

  void GetJetWithJEC(nano_tree &nano, std::vector<float> &Jet_pt, 
                     std::vector<float> &Jet_pt_jerUp, 
                     std::vector<float> &Jet_pt_jerDn, 
                     std::vector<float> &Jet_mass, 
                     std::vector<float> &Jet_mass_jerUp, 
                     std::vector<float> &Jet_mass_jerDn);

  int year;
  float nanoaod_version;
  bool verbose;
  float min_jet_pt;
  float max_jet_eta;
  bool isData;
  std::unique_ptr<correction::CorrectionSet> cs_jerc_;
  correction::Correction::Ref map_jes_;
  correction::Correction::Ref map_jersf_;
  correction::Correction::Ref map_jermc_;
  TRandom3 rng_;
};

#endif
