#ifndef H_ISR_TOOLS
#define H_ISR_TOOLS


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class ISRTools{
public:

  explicit ISRTools(const std::string &name, int year, float nanoaod_version);
  ~ISRTools();

  bool IsLastCopyBeforeFSR_or_LastCopy(nano_tree &nano, int mc_index);
  void WriteISRSystemPt(nano_tree &nano, pico_tree &pico);
  // Makes child_map
  // child_map[mom_index] = [child_index]
  std::map<int, std::vector<int> > GetChildMap(nano_tree &nano);
  void WriteISRJetMultiplicity(nano_tree &nano, pico_tree &pico);
  void WriteISRWeights(pico_tree &pico);

private:
  const std::string name;
  int year;
  bool isTTJets_LO;
  bool isGluino;
  bool isTChi;
  float nanoaod_version;

};

#endif
