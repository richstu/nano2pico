#ifndef H_EVENT_TOOLS
#define H_EVENT_TOOLS


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class EventTools{
public:

  explicit EventTools(const std::string &name, int year);
  ~EventTools();

  enum Dataset {EGamma, SingleElectron, SingleMuon, DoubleEG, DoubleMuon, MET, JetHT};

  void WriteStitch(nano_tree &nano, pico_tree &pico);
  void WriteDataQualityFilters(nano_tree& nano, pico_tree& pico, std::vector<int> sig_jet_nano_idx,
                               float min_jet_pt, float max_jet_pt, bool isData, bool isFastsim);
  bool SaveTriggerDecisions(nano_tree& nano, pico_tree& pico, bool isZgamma);
  void WriteTriggerEfficiency(pico_tree &pico);
  int GetEventType();


private:
  const std::string name;
  int year;
  bool isTTJets_LO_MET;
  bool isTTJets_LO_Incl;
  bool isWJets_LO;
  bool isDYJets_LO;
  int dataset;
};

#endif