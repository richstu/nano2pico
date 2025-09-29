#ifndef H_ZGAMMA_PRODUCER
#define H_ZGAMMA_PRODUCER

#include <vector>

#include "TLorentzVector.h"

#include "nano_tree.hpp"
#include "pico_tree.hpp"
#include "KinZfitter.hpp"

class ZGammaVarProducer{
public:

  explicit ZGammaVarProducer(int year);
  ~ZGammaVarProducer();

  void WriteZGammaVars(nano_tree &nano, pico_tree &pico, 
      std::vector<int> sig_jet_nano_idx, bool is_signal);

  struct RefitResults {
    TLorentzVector l1;
    TLorentzVector l2;
    int status;
    int covmatstatus;
    float minnll;
  };

  RefitResults PerformKinematicRefit(pico_tree &pico, int ll_variation);

  // calculates kinematic angles {cosTheta, costheta, phi}
  std::vector<double> CalculateAngles(TLorentzVector lplus, 
      TLorentzVector lminus, TLorentzVector ph);

private:
  int year;
  KinZfitter *kinZfitter; 
  const float el_m = 0.000511f;
  const float mu_m = 0.10566f;
};

#endif
