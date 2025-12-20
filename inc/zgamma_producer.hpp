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

  void WriteZGammaVars(pico_tree &pico);

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
