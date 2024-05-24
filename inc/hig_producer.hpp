#ifndef H_HIG_PRODUCER
#define H_HIG_PRODUCER

#include <vector>

#include "TLorentzVector.h"

#include "nano_tree.hpp"
#include "pico_tree.hpp"

struct HiggsConstructionVariables {
  std::vector<TLorentzVector> jet_lv;
  std::vector<float> jet_deepcsv;
};

class HigVarProducer{
public:

  explicit HigVarProducer(int year);
  ~HigVarProducer();

  void WriteHigVars(pico_tree& pico, bool doDeepFlav, bool isSignal,
                    std::vector<HiggsConstructionVariables> sys_higvars, float nanoaod_version);

private:
  int year;
  
};

#endif
