#ifndef H_EVENT_TOOLS
#define H_EVENT_TOOLS


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class EventTools{
public:

  explicit EventTools(const std::string &name, int year, bool isData);
  ~EventTools();

  enum Dataset {EGamma, SingleElectron, SingleMuon, DoubleEG, DoubleMuon, MET, JetHT, Muon, JetMET, MuonEG};

  void WriteStitch(nano_tree &nano, pico_tree &pico);
  void WriteDataQualityFilters(nano_tree& nano, pico_tree& pico, std::vector<int> sig_jet_nano_idx,
                               float min_jet_pt, bool isFastsim, bool is_preUL);
  bool SaveTriggerDecisions(nano_tree& nano, pico_tree& pico, bool isZgamma);
  void WriteTriggerEfficiency(pico_tree &pico);
  int GetEventType();


private:
  const std::string name;
  int year;
  bool isTTJets_LO_Incl;
  bool isTTJets_LO_MET;
  bool isTTJets_LO_HT;
  bool isWJets_LO;
  bool isDYJets_LO;
  bool isWW;
  bool isWZ;
  bool isZZ;
  bool isFastSim;
  bool isData;
  int dataset;
};

#endif
