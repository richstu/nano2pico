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
      //cout<<"  below's overflow/underflow: "<<sf.first<<endl;
    }else{
      sf = {h.GetBinContent(bin), h.GetBinError(bin)};
    }
    //cout<<"eta: "<<x<<" pt: "<<y<<" overflow: "<<h.IsBinUnderflow(bin)<<" underflow: "<<h.IsBinUnderflow(bin)<<" binContent: "<<h.GetBinContent(bin)<<endl;
    if(ignore_error) sf.second = 0.;
    return sf;
  }
  template<typename T>
    pair<double, double> GetSF(const T &h, double x, bool ignore_error = false){
    pair<double, double> sf;
    auto bin = h.FindFixBin(x);
    if((h.IsBinOverflow(bin) || h.IsBinUnderflow(bin))
       && h.GetBinContent(bin) == 0. && h.GetBinError(bin) == 0.){
      auto bin_x = h.GetXaxis()->FindFixBin(x);
      if(bin_x <= 0) bin_x = 1;
      if(bin_x > h.GetNbinsX()) bin_x = h.GetNbinsX();
      sf = {h.GetBinContent(bin_x), h.GetBinError(bin_x)};
      //cout<<"  [oneVar] below's overflow/underflow: "<<sf.first<<endl;
    }else{
      sf = {h.GetBinContent(bin), h.GetBinError(bin)};
    }
    //cout<<"[oneVar] eta: "<<x<<" overflow: "<<h.IsBinUnderflow(bin)<<" underflow: "<<h.IsBinUnderflow(bin)<<" binContent: "<<h.GetBinContent(bin)<<endl;
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
  year_ = year;
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
    //cout<<"Photon["<<i<<"] sf: "<<GetPhotonScaleFactor(pico, i).first<<endl;
    if(pico.out_photon_sig().at(i)){
      sf = MergeSF(sf, GetPhotonScaleFactor(pico, i));
    }
  }
  w_photon = sf.first;
  sys_photon = vector<float>{static_cast<float>(sf.first+sf.second),
                             static_cast<float>(sf.first-sf.second)};
}

// Returns 1.5(EB high R9), 2.5(EB low R9), 4.5(EE high R9), 5.5(EE low R9)
float PhotonWeighter::GetRegion(float const & eta, float const & r9) {
  if (fabs(eta) < 1.4442f) { //EB
    if (fabs(r9) > 0.94f) return 1.5;
    else return 2.5;
  } else if (fabs(eta) < 2.5f) {
    if (fabs(r9) > 0.94f) return 4.5;
    else return 5.5;
  } else return -0.5;
}

std::pair<double, double> PhotonWeighter::GetPhotonScaleFactor(pico_tree &pico, size_t iph){
  //https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations#E_gamma_Energy_Corrections
  //ID and Electron Veto SFs applied
  double pt = pico.out_photon_pt().at(iph);
  double eta = pico.out_photon_eta().at(iph);
  vector<pair<double, double> > sfs;
  //Axes swapped, asymmetric in eta
  if (do_full_photon_id_) sfs.push_back(GetSF(sf_full_photon_id_, eta, pt));
  if (do_full_photon_ev_) {
    if (year_ != 2017) sfs.push_back(GetSF(sf_full_photon_ev_, eta, pt));
    else sfs.push_back(GetSF(sf_full_photon_ev_, GetRegion(pico.out_photon_eta().at(iph), pico.out_photon_r9().at(iph))));
  }
  sfs.push_back(make_pair(1., pt<20. || pt >80. ? 0.01 : 0.));//Systematic uncertainty
  return accumulate(sfs.cbegin(), sfs.cend(), make_pair(1., 0.), MergeSF);
}

