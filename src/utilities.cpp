//----------------------------------------------------------------------------
// utilities - Various functions used accross the code
//----------------------------------------------------------------------------

#include "utilities.hpp"

#include <cmath>

#include <deque>
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <iomanip>   // setw

#include <libgen.h>

#include "TCollection.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TList.h"
#include "TString.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TRegexp.h"
#include "TLorentzVector.h"

#include "nano_tree.hpp"

using namespace std;

long double DeltaPhi(long double phi1, long double phi2){
  long double dphi = fmod(fabs(phi2-phi1), 2.L*PI);
  return dphi>PI ? 2.L*PI-dphi : dphi;
}

long double SignedDeltaPhi(long double phi1, long double phi2){
  long double dphi = fmod(phi2-phi1, 2.L*PI);
  if(dphi>PI){
    return dphi-2.L*PI;
  }else if(dphi<-PI){
    return dphi+2.L*PI;
  }else{
    return dphi;
  }
}

float dR(float eta1, float eta2, float phi1, float phi2) {
  return AddInQuadrature(eta1-eta2, DeltaPhi(phi1,phi2));
}

double cosThetaJeff(TLorentzVector lminus, TLorentzVector lplus, TLorentzVector photon) {
  // Calculates the angle between the Zs spin and the lepton in the Zs rest frame
  TLorentzVector ll = lminus + lplus;
  lminus.Boost(-ll.BoostVector());
  photon.Boost(-ll.BoostVector());
  TVector3 l(lminus.Vect()), p(photon.Vect());
  double costj = l*p/(l.Mag()*p.Mag());
  return costj;
}

TString roundNumber(double num, int decimals, double denom){
  if(denom==0) return " - ";
  double neg = 1; if(num*denom<0) neg = -1;
  num /= neg*denom; num += 0.5*pow(10.,-decimals);
  long num_int = static_cast<long>(num);
  long num_dec = static_cast<long>((1+num-num_int)*pow(10.,decimals));
  TString s_dec = ""; s_dec += num_dec; s_dec.Remove(0,1);
  TString result="";
  if(neg<0) result+="-";
  result+= num_int;
  if(decimals>0) {
    result+="."; result+=s_dec;
  }

  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<decimals-afterdot.Length(); i++)
    result += "0";
  return result;
}

TString addCommas(double num){
  TString result(""); result += num;
  int posdot(result.First('.'));
  if(posdot==-1) posdot = result.Length();
  for(int ind(posdot-3); ind > 0; ind -= 3)
    result.Insert(ind, ",");
  return result;
}

long double AddInQuadrature(long double x, long double y){
  if(fabs(y)>fabs(x)){
    const long double temp = y;
    y=x;
    x=temp;
  }
  if(x==0.) return y;
  const long double rat=y/x;
  return fabs(x)*sqrt(1.0L+rat*rat);
}

long double GetMass(long double e, long double px, long double py, long double pz){
  px/=e; py/=e; pz/=e;
  return fabs(e)*sqrt(1.0L-px*px-py*py-pz*pz);
}

long double GetMT(long double m1, long double pt1, long double phi1,
                  long double m2, long double pt2, long double phi2){
  return sqrt(m1*m1+m2*m2+2.L*(sqrt((m1*m1+pt1*pt1)*(m2*m2+pt2*pt2))-pt1*pt2*cos(phi2-phi1)));
}

long double GetMT(long double pt1, long double phi1,
                  long double pt2, long double phi2){
  //Faster calculation in massless case
  return sqrt(2.L*pt1*pt2*(1.L-cos(phi2-phi1)));
}

bool Contains(const string& text, const string& pattern){
  return text.find(pattern) != string::npos;
}

vector<string> Tokenize(const string& input,
                        const string& tokens){
  char* ipt(new char[input.size()+1]);
  memcpy(ipt, input.data(), input.size());
  ipt[input.size()]=static_cast<char>(0);
  char* ptr(strtok(ipt, tokens.c_str()));
  vector<string> output(0);
  while(ptr!=NULL){
    output.push_back(ptr);
    ptr=strtok(NULL, tokens.c_str());
  }
  return output;
}

void get_count_and_uncertainty(TTree& tree,
                               const string& cut,
                               double& count,
                               double& uncertainty){
  const string hist_name("temp");
  TH1D temp(hist_name.c_str(), "", 1, -1.0, 1.0);
  tree.Project(hist_name.c_str(), "0.0", cut.c_str());
  count=temp.IntegralAndError(0,2,uncertainty);
}

void AddPoint(TGraph& graph, const double x, const double y){
  graph.SetPoint(graph.GetN(), x, y);
}

string execute(const string &cmd){
  FILE *pipe = popen(cmd.c_str(), "r");
  if(!pipe) throw runtime_error("Could not open pipe.");
  const size_t buffer_size = 128;
  char buffer[buffer_size];
  string result = "";
  while(!feof(pipe)){
    if(fgets(buffer, buffer_size, pipe) != NULL) result += buffer;
  }

  pclose(pipe);
  return result;
}

TString hoursMinSec(long seconds){
  int minutes((seconds/60)%60), hours(seconds/3600);
  TString hhmmss("");
  if(hours<10) hhmmss += "0";
  hhmmss += hours; hhmmss += ":";
  if(minutes<10) hhmmss += "0";
  hhmmss += minutes; hhmmss += ":";
  if((seconds%60)<10) hhmmss += "0";
  hhmmss += seconds%60; 

  return hhmmss;
}

void ReplaceAll(string &str, const string &orig, const string &rep){
  size_t loc = 0;
  while ((loc = str.find(orig, loc)) != string::npos) {
    str.replace(loc, orig.length(), rep);
    loc += rep.length();
  }
}

string CopyReplaceAll(string str, const string &orig, const string &rep){
  ReplaceAll(str, orig, rep);
  return str;
}

void SplitFilePath(const string &path, string &dir_name, string &base_name){
  vector<char> cstr(path.c_str(), path.c_str()+path.size()+1);
  dir_name = dirname(&cstr.at(0));
  cstr = vector<char>(path.c_str(), path.c_str()+path.size()+1);
  base_name = basename(&cstr.at(0));
}

//simple propagation of Gaussian uncertainties based on linear Taylor expansion
void propagate_uncertainty_product(float a, float a_unc, float b, float b_unc, float& prod, float& prod_unc) {
  prod = a*b;
  prod_unc = hypotf(a_unc*b,b_unc*a);
}

//simple propagation of Gaussian uncertainties based on linear Taylor expansion
void propagate_uncertainty_ratio(float num, float num_unc, float den, float den_unc, float& ratio, float& ratio_unc) {
  ratio = 1.0;
  ratio_unc = 1.0;
  if (den != 0.0) {
    ratio = num/den;
    ratio_unc = hypotf(num_unc/den, den_unc*num/den/den);
  }
  //else if (den_unc != 0.0) {
  //  ratio = num/(den_unc/2.0);
  //  ratio_unc = fabs(ratio-num/den_unc);
  //}
}

void getMETWithJEC(nano_tree & nano, int year, bool isFastsim, float & MET_pt, float & MET_phi, bool is_preUL) {
  if (isFastsim) { 
    if (year==2017 && is_preUL) {
      MET_pt = nano.METFixEE2017_T1_pt(); 
      MET_phi = nano.METFixEE2017_T1_phi();
    } else {
      MET_pt = nano.MET_T1_pt();
      MET_phi = nano.MET_T1_phi();
    }
  } else {
    if (year == 2017 && is_preUL) {
      MET_pt = nano.METFixEE2017_pt();
      MET_phi = nano.METFixEE2017_phi();
    } else {
      MET_pt = nano.MET_pt();
      MET_phi = nano.MET_phi();
    }
  }
}
void getJetWithJEC(nano_tree & nano, bool isFastsim, vector<float> & Jet_pt, vector<float> & Jet_mass) {
  Jet_pt.resize(nano.nJet());
  Jet_mass.resize(nano.nJet());
  for(int ijet(0); ijet<nano.nJet(); ++ijet){
    if (isFastsim) {
      Jet_pt[ijet] = nano.Jet_pt_nom()[ijet];
      Jet_mass[ijet] = nano.Jet_mass_nom()[ijet];
    } else {
      Jet_pt[ijet] = nano.Jet_pt()[ijet];
      Jet_mass[ijet] = nano.Jet_mass()[ijet];
    }
  }
}

void getJetId(nano_tree & nano, float nanoaod_version, vector<int> & Jet_jetId) {
  Jet_jetId.resize(nano.nJet());
  for(int ijet(0); ijet<nano.nJet(); ++ijet){
    if (nanoaod_version+0.01 > 11.9) Jet_jetId[ijet] = nano.Jet_jetId_11p9()[ijet];
    else Jet_jetId[ijet] = nano.Jet_jetId()[ijet];
  }
}

void getFatJet_btagDDBvL(nano_tree & nano, float nanoaod_version, vector<float> & FatJet_btagDDBvL) {
  FatJet_btagDDBvL.resize(nano.nFatJet());
  for(int ijet(0); ijet<nano.nFatJet(); ++ijet){
    if (nanoaod_version+0.01 < 9) FatJet_btagDDBvL[ijet] = nano.FatJet_btagDDBvL()[ijet];
    else FatJet_btagDDBvL[ijet] = nano.FatJet_btagDDBvLV2()[ijet];
  }
}

void getFatJet_particleNetWithMass_WvsQCD(nano_tree & nano, float nanoaod_version, 
                                          vector<float> & FatJet_particleNetWithMass_WvsQCD) {
  FatJet_particleNetWithMass_WvsQCD.resize(nano.nFatJet());
  for(int ijet(0); ijet<nano.nFatJet(); ++ijet){
    if (nanoaod_version+0.01 < 11.9) 
      FatJet_particleNetWithMass_WvsQCD[ijet] = nano.FatJet_particleNet_WvsQCD()[ijet];
    else 
      FatJet_particleNetWithMass_WvsQCD[ijet] = nano.FatJet_particleNetWithMass_WvsQCD()[ijet];
  }
}

void getFatJet_particleNetWithMass_ZvsQCD(nano_tree & nano, float nanoaod_version, 
                                          vector<float> & FatJet_particleNetWithMass_ZvsQCD) {
  FatJet_particleNetWithMass_ZvsQCD.resize(nano.nFatJet());
  for(int ijet(0); ijet<nano.nFatJet(); ++ijet){
    if (nanoaod_version+0.01 < 11.9) 
      FatJet_particleNetWithMass_ZvsQCD[ijet] = nano.FatJet_particleNet_ZvsQCD()[ijet];
    else 
      FatJet_particleNetWithMass_ZvsQCD[ijet] = nano.FatJet_particleNetWithMass_ZvsQCD()[ijet];
  }
}

void getFatJet_particleNetWithMass_TvsQCD(nano_tree & nano, float nanoaod_version, 
                                          vector<float> & FatJet_particleNetWithMass_TvsQCD) {
  FatJet_particleNetWithMass_TvsQCD.resize(nano.nFatJet());
  for(int ijet(0); ijet<nano.nFatJet(); ++ijet){
    if (nanoaod_version+0.01 < 11.9) 
      FatJet_particleNetWithMass_TvsQCD[ijet] = nano.FatJet_particleNet_TvsQCD()[ijet];
    else 
      FatJet_particleNetWithMass_TvsQCD[ijet] = nano.FatJet_particleNetWithMass_TvsQCD()[ijet];
  }
}

void getFatJet_particleNet_mass(nano_tree & nano, float nanoaod_version, 
                                          vector<float> & FatJet_particleNet_mass) {
  FatJet_particleNet_mass.resize(nano.nFatJet());
  for(int ijet(0); ijet<nano.nFatJet(); ++ijet){
    if (nanoaod_version+0.01 < 11.9) 
      FatJet_particleNet_mass[ijet] = nano.FatJet_particleNet_mass()[ijet];
    else 
      FatJet_particleNet_mass[ijet] = nano.FatJet_particleNet_massCorr()[ijet];
  }
}

void getPhoton_electronIdx(nano_tree & nano, float nanoaod_version, vector<int> & Photon_electronIdx) {
  Photon_electronIdx.resize(nano.nPhoton());
  for(int iphoton(0); iphoton<nano.nPhoton(); ++iphoton){
    if (nanoaod_version+0.01 > 11.9) Photon_electronIdx[iphoton] = nano.Photon_electronIdx_11p9()[iphoton];
    else Photon_electronIdx[iphoton] = nano.Photon_electronIdx()[iphoton];
  }
}

void getMuon_fsrPhotonIdx(nano_tree & nano, float nanoaod_version, vector<int> & Muon_fsrPhotonIdx) {
  Muon_fsrPhotonIdx.resize(nano.nMuon());
  for(int imuon(0); imuon<nano.nMuon(); ++imuon){
    if (nanoaod_version+0.01 > 11.9) Muon_fsrPhotonIdx[imuon] = nano.Muon_fsrPhotonIdx_11p9()[imuon];
    else Muon_fsrPhotonIdx[imuon] = nano.Muon_fsrPhotonIdx()[imuon];
  }
}

void getElectron_photonIdx(nano_tree & nano, float nanoaod_version, vector<int> & Electron_photonIdx) {
  Electron_photonIdx.resize(nano.nElectron());
  //cout<<"nElectron: "<<nano.nElectron()<<endl;
  //cout<<"  "<<nano.Electron_photonIdx()[0]<<endl;
  ////cout<<"  "<<nano.Electron_photonIdx_short()[0]<<endl;
  //cout<<"nPhoton: "<<nano.nPhoton()<<endl;
  //cout<<"  "<<nano.Photon_cutBased()[0]<<endl;
  //cout<<"  "<<nano.Photon_cutBased_char()[0]<<endl;
  for(int iel(0); iel<nano.nElectron(); ++iel){
    if (nanoaod_version+0.01 > 11.9) Electron_photonIdx[iel] = nano.Electron_photonIdx_11p9()[iel];
    else Electron_photonIdx[iel] = nano.Electron_photonIdx()[iel];
  }
}

void getFsrPhoton_muonIdx(nano_tree & nano, float nanoaod_version, vector<int> & FsrPhoton_muonIdx) {
  FsrPhoton_muonIdx.resize(nano.nFsrPhoton());
  for(int ipart(0); ipart<nano.nFsrPhoton(); ++ipart){
    if (nanoaod_version+0.01 > 11.9) FsrPhoton_muonIdx[ipart] = nano.FsrPhoton_muonIdx_11p9()[ipart];
    else FsrPhoton_muonIdx[ipart] = nano.FsrPhoton_muonIdx()[ipart];
  }
}

void getPhoton_jetIdx(nano_tree & nano, float nanoaod_version, vector<int> & Photon_jetIdx) {
  Photon_jetIdx.resize(nano.nPhoton());
  for(int ipart(0); ipart<nano.nPhoton(); ++ipart){
    if (nanoaod_version+0.01 > 11.9) Photon_jetIdx[ipart] = nano.Photon_jetIdx_11p9()[ipart];
    else Photon_jetIdx[ipart] = nano.Photon_jetIdx()[ipart];
  }
}

void getPhoton_cutBased(nano_tree & nano, float nanoaod_version, vector<int> & Photon_cutBased) {
  Photon_cutBased.resize(nano.nPhoton());
  for(int ipart(0); ipart<nano.nPhoton(); ++ipart){
    if (nanoaod_version+0.01 > 11.9) Photon_cutBased[ipart] = nano.Photon_cutBased_11p9()[ipart];
    else Photon_cutBased[ipart] = nano.Photon_cutBased()[ipart];
  }
}

void getFatJet_subJetIdx1(nano_tree & nano, float nanoaod_version, vector<int> & FatJet_subJetIdx1) {
  FatJet_subJetIdx1.resize(nano.nFatJet());
  for(int ipart(0); ipart<nano.nFatJet(); ++ipart){
    if (nanoaod_version+0.01 > 11.9) FatJet_subJetIdx1[ipart] = nano.FatJet_subJetIdx1_11p9()[ipart];
    else FatJet_subJetIdx1[ipart] = nano.FatJet_subJetIdx1()[ipart];
  }
}

void getFatJet_subJetIdx2(nano_tree & nano, float nanoaod_version, vector<int> & FatJet_subJetIdx2) {
  FatJet_subJetIdx2.resize(nano.nFatJet());
  for(int ipart(0); ipart<nano.nFatJet(); ++ipart){
    if (nanoaod_version+0.01 > 11.9) FatJet_subJetIdx2[ipart] = nano.FatJet_subJetIdx2_11p9()[ipart];
    else FatJet_subJetIdx2[ipart] = nano.FatJet_subJetIdx2()[ipart];
  }
}

void getMuon_nTrackerLayers(nano_tree & nano, float nanoaod_version, vector<int> & Muon_nTrackerLayers) {
  Muon_nTrackerLayers.resize(nano.nMuon());
  for(int ipart(0); ipart<nano.nMuon(); ++ipart){
    if (nanoaod_version+0.01 > 11.9) Muon_nTrackerLayers[ipart] = nano.Muon_nTrackerLayers_11p9()[ipart];
    else Muon_nTrackerLayers[ipart] = nano.Muon_nTrackerLayers()[ipart];
  }
}

void getMuon_genPartIdx(nano_tree & nano, float nanoaod_version, vector<int> & Muon_genPartIdx) {
  Muon_genPartIdx.resize(nano.nMuon());
  for(int ipart(0); ipart<nano.nMuon(); ++ipart){
    if (nanoaod_version+0.01 > 12) Muon_genPartIdx[ipart] = nano.Muon_genPartIdx_12p0()[ipart];                                                                                                                                          
    else Muon_genPartIdx[ipart] = nano.Muon_genPartIdx()[ipart];                                                                                                                                                                         
  }
}

void getJet_genJetIdx(nano_tree & nano, float nanoaod_version, vector<int> & Jet_genJetIdx) {
  Jet_genJetIdx.resize(nano.nJet());
  for(int ipart(0); ipart<nano.nJet(); ++ipart){
    if (nanoaod_version+0.01 > 12) Jet_genJetIdx[ipart] = nano.Jet_genJetIdx_12p0()[ipart];
    else Jet_genJetIdx[ipart] = nano.Jet_genJetIdx()[ipart];
  }
}

void getJet_hadronFlavour(nano_tree & nano, float nanoaod_version, vector<int> & Jet_hadronFlavour) {
  Jet_hadronFlavour.resize(nano.nJet());
  for(int ipart(0); ipart<nano.nJet(); ++ipart){
    if (nanoaod_version+0.01 > 12) Jet_hadronFlavour[ipart] = nano.Jet_hadronFlavour_12p0()[ipart];
    else Jet_hadronFlavour[ipart] = nano.Jet_hadronFlavour()[ipart];
  }
}

void getJet_partonFlavour(nano_tree & nano, float nanoaod_version, vector<int> & Jet_partonFlavour) {
  Jet_partonFlavour.resize(nano.nJet());
  for(int ipart(0); ipart<nano.nJet(); ++ipart){
    if (nanoaod_version+0.01 > 12) Jet_partonFlavour[ipart] = nano.Jet_partonFlavour_12p0()[ipart];
    else Jet_partonFlavour[ipart] = nano.Jet_partonFlavour()[ipart];
  }
}

void getGenPart_genPartIdxMother(nano_tree & nano, float nanoaod_version, vector<int> & GenPart_genPartIdxMother) {
  GenPart_genPartIdxMother.resize(nano.nGenPart());
  for(int ipart(0); ipart<nano.nGenPart(); ++ipart){
    if (nanoaod_version+0.01 > 12) GenPart_genPartIdxMother[ipart] = nano.GenPart_genPartIdxMother_12p0()[ipart];
    else GenPart_genPartIdxMother[ipart] = nano.GenPart_genPartIdxMother()[ipart];
  }
}

void getGenPart_statusFlags(nano_tree & nano, float nanoaod_version, vector<int> & GenPart_statusFlags) {
  GenPart_statusFlags.resize(nano.nGenPart());                                                                                                                                                                                                    
  for(int ipart(0); ipart<nano.nGenPart(); ++ipart){
    if (nanoaod_version+0.01 > 12) GenPart_statusFlags[ipart] = nano.GenPart_statusFlags_12p0()[ipart];                                                                                                                                  
    else GenPart_statusFlags[ipart] = nano.GenPart_statusFlags()[ipart];                                                                                                                                                                 
  }
}

void getGenJet_partonFlavour(nano_tree & nano, float nanoaod_version, vector<int> & GenJet_partonFlavour) {
  GenJet_partonFlavour.resize(nano.nGenJet());
  for(int ipart(0); ipart<nano.nGenJet(); ++ipart){
    if (nanoaod_version+0.01 > 12) GenJet_partonFlavour[ipart] = nano.GenJet_partonFlavour_12p0()[ipart];
    else GenJet_partonFlavour[ipart] = nano.GenJet_partonFlavour()[ipart];
  }
}
