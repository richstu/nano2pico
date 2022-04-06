#include "photon_weighter.hpp"

#include <string>
#include <numeric>

#include "TFile.h"
#include "TGraphAsymmErrors.h"

#include "utilities.hpp"

using namespace std;

namespace{
  template<typename T>
    T LoadSF(const string &file_name, const string &item_name){
    string path = "data/"+file_name;
    TFile f(path.c_str(), "read");
    if(!f.IsOpen()) ERROR("Could not open "+file_name);
    T* item = static_cast<T*>(f.Get(item_name.c_str()));
    if(!item) ERROR("Could not find "+item_name+" in "+file_name);
    return *item;
  }
  template<typename T>
    pair<double, double> GetSF(const T &h, double x, double y, bool ignore_error = false){
    pair<double, double> sf;
    auto bin = h.FindFixBin(x, y);
    if((h.IsBinOverflow(bin) || h.IsBinUnderflow(bin))
       && h.GetBinContent(bin) == 0. && h.GetBinError(bin) == 0.){
      auto bin_x = h.GetXaxis()->FindFixBin(x);
      auto bin_y = h.GetYaxis()->FindFixBin(y);
      if(bin_x <= 0) bin_x = 1;
      if(bin_x > h.GetNbinsX()) bin_x = h.GetNbinsX();
      if(bin_y <= 0) bin_y = 1;
      if(bin_y > h.GetNbinsY()) bin_y = h.GetNbinsY();
      sf = {h.GetBinContent(bin_x, bin_y), h.GetBinError(bin_x, bin_y)};
    }else{
      sf = {h.GetBinContent(bin), h.GetBinError(bin)};
    }
    if(ignore_error) sf.second = 0.;
    return sf;
  }

  pair<double, double> MergeSF(pair<double, double> a,
                               pair<double, double> b){
    double sf = a.first * b.first;
    double err = hypot(a.first*b.second, b.first*a.second);
    return {sf, err};
  }
}

PhotonWeighter::PhotonWeighter(int year, bool isZgamma){
  if (isZgamma) {
    if (year==2016) {
      in_full_photon_id_ = "zgamma/Fall17V2_2016_MVAwp90_photons.root"; hist_full_photon_id_ = "EGamma_SF2D";
      in_full_photon_ev_ = "zgamma/Photon_ElVetoSFs_80X_Summer16.root"; hist_full_photon_ev_ = "Scaling_Factors_CSEV_R9 Inclusive";
    } else if (year==2017) {
      in_full_photon_id_ = "zgamma/2017_PhotonsMVAwp90.root";    hist_full_photon_id_ = "EGamma_SF2D";
      in_full_photon_ev_ = "zgamma/CSEV_ScaleFactors_2017.root"; hist_full_photon_ev_ = "MVA_ID";
    } else {
      in_full_photon_id_ = "zgamma/2018_PhotonsMVAwp90.root"; hist_full_photon_id_ = "EGamma_SF2D";
      in_full_photon_ev_ = "zgamma/CSEV_2018.root";           hist_full_photon_ev_ = "eleVeto_SF";
    }
    do_full_photon_id_ = (in_full_photon_id_!=""); do_full_photon_ev_ = (in_full_photon_ev_!="");

    if (do_full_photon_id_) sf_full_photon_id_ = LoadSF<TH2F>(in_full_photon_id_, hist_full_photon_id_);
    if (do_full_photon_ev_) sf_full_photon_ev_ = LoadSF<TH2F>(in_full_photon_ev_, hist_full_photon_ev_);
  }
}

void PhotonWeighter::FullSim(pico_tree &pico, float &w_photon, vector<float> &sys_photon){
  pair<double, double> sf(1., 0.);
  for(size_t i = 0; i < pico.out_photon_sig().size(); ++i){
    if(pico.out_photon_sig().at(i)){
      sf = MergeSF(sf, GetPhotonScaleFactor(pico, i));
    }
  }
  w_photon = sf.first;
  sys_photon = vector<float>{static_cast<float>(sf.first+sf.second),
                             static_cast<float>(sf.first-sf.second)};
}

std::pair<double, double> PhotonWeighter::GetPhotonScaleFactor(pico_tree &pico, size_t iph){
  //https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations#E_gamma_Energy_Corrections
  //ID and Electron Veto SFs applied
  double pt = pico.out_photon_pt().at(iph);
  double eta = pico.out_photon_eta().at(iph);
  vector<pair<double, double> > sfs;
  //Axes swapped, asymmetric in eta
  if (do_full_photon_id_) sfs.push_back(GetSF(sf_full_photon_id_, eta, pt));
  if (do_full_photon_ev_) sfs.push_back(GetSF(sf_full_photon_ev_, eta, pt));
  sfs.push_back(make_pair(1., pt<20. || pt >80. ? 0.01 : 0.));//Systematic uncertainty
  return accumulate(sfs.cbegin(), sfs.cend(), make_pair(1., 0.), MergeSF);
}

