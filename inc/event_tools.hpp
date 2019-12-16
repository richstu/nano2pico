#ifndef H_EVENT_TOOLS
#define H_EVENT_TOOLS


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class EventTools{
public:

  explicit EventTools(const std::string &name, int year);
  ~EventTools();

  void WriteStitch(nano_tree &nano, pico_tree &pico);
  void WriteDataQualityFilters(nano_tree& nano, pico_tree& pico, std::vector<int> &sig_jet_nano_idx,
                               bool isData, bool isFastsim);
  void CopyTriggerDecisions(nano_tree& nano, pico_tree& pico);
  void WriteTriggerEfficiency(pico_tree &pico);
  int GetEventType();


private:
  const std::string name;
  int year;
  bool isTTJets_LO;
  bool isWJets_LO;
  bool isDYJets_LO;

};

#endif