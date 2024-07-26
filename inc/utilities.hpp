//----------------------------------------------------------------------------
// utilities - Various functions used accross the code
//----------------------------------------------------------------------------

#ifndef H_UTILITIES
#define H_UTILITIES

#include <cstddef>
#include <cstdio>
#include <cmath>

#include <iostream>
#include <string>
#include <vector>
#include <set>

#include <unistd.h>

#include "TString.h"
#include "TTree.h"
#include "TGraph.h"
#include "TLorentzVector.h"

#include "nano_tree.hpp"

#define ERROR(x) do{throw std::runtime_error(std::string("Error in file ")+__FILE__+" at line "+std::to_string(__LINE__)+" (in "+__func__+"): "+x);}while(false)
#define DBG(x) do{std::cerr << "In " << __FILE__ << " at line " << __LINE__ << " (in function " << __func__ << "): " << x << std::endl;}while(false)

class Variable{
public:
  Variable():
    type_(""),
    name_(""){
  }

  Variable( std::string &type,
            std::string &name):
    type_(type),
    name_(name){
  }

  bool operator<( Variable& var) {
    return type_<var.type_ || (type_==var.type_ && name_<var.name_);
  }

  std::string type_, name_;
};

/*!\brief struct used to sort nano entries
 */
struct NanoOrderEntry {
  int nano_idx;
  float pt;
  bool is_sig;
};

const long double PI = acos(-1.L);
long double DeltaPhi(long double phi1, long double phi2);
long double SignedDeltaPhi(long double phi1, long double phi2);
float dR(float eta1, float eta2, float phi1, float phi2);
double cosThetaJeff(TLorentzVector lminus, TLorentzVector lplus, TLorentzVector photon);
TString roundNumber(double num, int decimals, double denom=1.);
TString addCommas(double num);
long double AddInQuadrature(long double x, long double y);
long double GetMass(long double e, long double px, long double py, long double pz);
long double GetMT(long double m1, long double pt1, long double phi1,
                  long double m2, long double pt2, long double phi2);
long double GetMT(long double pt1, long double phi1,
                  long double pt2, long double phi2);
bool Contains(const std::string& text, const std::string& pattern);

std::vector<std::string> Tokenize(const std::string& input,
                                  const std::string& tokens=" ");
void get_count_and_uncertainty(TTree& tree,
                               const std::string& cut,
                               double& count,
                               double& uncertainty);
void AddPoint(TGraph& graph, const double x, const double y);
TString hoursMinSec(long seconds);

template<class T>
T noInfNan(T val, T defval=1.){
  if(isnan(val)) {
    std::cout<<"Value is NaN. Returning "<<defval<<std::endl;
    return defval;
  } else if(isinf(val)) {
    std::cout<<"Value is Inf. Returning "<<defval<<std::endl;
    return defval;
  } else return val;
}

template<class T>
short Sign(T val){
  return (T(0) < val) - (val < T(0));
}

std::string execute(const std::string &cmd);

void ReplaceAll(std::string &str, const std::string &orig, const std::string &rep);
std::string CopyReplaceAll(const std::string str, const std::string &orig, const std::string &rep);

void SplitFilePath(const std::string &path, std::string &dir_name, std::string &base_name);

void propagate_uncertainty_product(float a, float a_unc, float b, float b_unc, float& prod, float& prod_unc);
void propagate_uncertainty_ratio(float num, float num_unc, float den, float den_unc, float& ratio, float& ratio_unc);

void getMETWithJEC(nano_tree & nano, int year, bool isFastsim, float & MET_pt, float & MET_phi, bool is_preUL);
void getJetWithJEC(nano_tree & nano, bool isFastsim, std::vector<float> & Jet_pt, std::vector<float> & Jet_mass);
void getJetId(nano_tree & nano, float nanoaod_version, std::vector<int> & Jet_jetId);
void getFatJet_btagDDBvL(nano_tree & nano, float nanoaod_version, std::vector<float> & FatJet_btagDDBvL);
void getFatJet_particleNetWithMass_WvsQCD(nano_tree & nano, float nanoaod_version, std::vector<float> & FatJet_particleNetWithMass_WvsQCD);
void getFatJet_particleNetWithMass_ZvsQCD(nano_tree & nano, float nanoaod_version, std::vector<float> & FatJet_particleNetWithMass_ZvsQCD);
void getFatJet_particleNetWithMass_TvsQCD(nano_tree & nano, float nanoaod_version, std::vector<float> & FatJet_particleNetWithMass_TvsQCD);
void getFatJet_particleNet_mass(nano_tree & nano, float nanoaod_version, std::vector<float> & FatJet_particleNet_mass);
void getPhoton_electronIdx(nano_tree & nano, float nanoaod_version, std::vector<int> & Photon_electronIdx);
void getMuon_fsrPhotonIdx(nano_tree & nano, float nanoaod_version, std::vector<int> & Muon_fsrPhotonIdx);
void getElectron_photonIdx(nano_tree & nano, float nanoaod_version, std::vector<int> & Electron_photonIdx);
void getFsrPhoton_muonIdx(nano_tree & nano, float nanoaod_version, std::vector<int> & FsrPhoton_muonIdx);
void getPhoton_jetIdx(nano_tree & nano, float nanoaod_version, std::vector<int> & Photon_jetIdx);
void getPhoton_cutBased(nano_tree & nano, float nanoaod_version, std::vector<int> & Photon_cutBased);
void getFatJet_subJetIdx1(nano_tree & nano, float nanoaod_version, std::vector<int> & FatJet_subJetIdx1);
void getFatJet_subJetIdx2(nano_tree & nano, float nanoaod_version, std::vector<int> & FatJet_subJetIdx2);
void getMuon_nTrackerLayers(nano_tree & nano, float nanoaod_version, std::vector<int> & Muon_nTrackerLayers);
void getMuon_genPartIdx(nano_tree & nano, float nanoaod_version, std::vector<int> & Muon_genPartIdx);
void getJet_genJetIdx(nano_tree & nano, float nanoaod_version, std::vector<int> & Jet_genJetIdx);
void getJet_hadronFlavour(nano_tree & nano, float nanoaod_version, std::vector<int> & Jet_hadronFlavour);
void getJet_partonFlavour(nano_tree & nano, float nanoaod_version, std::vector<int> & Jet_partonFlavour);
void getGenPart_genPartIdxMother(nano_tree & nano, float nanoaod_version, std::vector<int> & GenPart_genPartIdxMother);
void getGenPart_statusFlags(nano_tree & nano, float nanoaod_version, std::vector<int> & GenPart_statusFlags);
void getGenJet_partonFlavour(nano_tree & nano, float nanoaod_version, std::vector<int> & GenJet_partonFlavour);
#endif
