// MC cross sections

#ifndef H_CROSS_SECTIONS
#define H_CROSS_SECTIONS

#include "TString.h"

namespace xsec{

  float crossSection(const TString &file, int year);
  void gluinoCrossSection(int glu_mass, double &xsec, double &xsec_unc, int year);
  void higgsinoCrossSection(int hig_mass, double &xsec, double &xsec_unc, int year);
  float fractionNegWeights(const TString &file);
}

#endif
