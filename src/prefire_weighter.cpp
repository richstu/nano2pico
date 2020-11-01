#include "prefire_weighter.hpp"

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
    //note that weights are stored as prefire weights whereas we want Non-prefire weights
    sf.first = 1.0-sf.first;
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

PrefireWeighter::PrefireWeighter(int year, bool use_jet_empt){
  use_jet_empt_ = use_jet_empt;

  if (year == 2016) {
    prefire_jet_empt_filename_ = "L1prefiring_jetempt_2016BtoH.root";         prefire_jet_empt_histname_ = "L1prefiring_jetempt_2016BtoH";
    prefire_jet_filename_ = "L1prefiring_jetpt_2016BtoH.root";                prefire_jet_histname_ = "L1prefiring_jetpt_2016BtoH";
    prefire_photon_filename_ = "L1prefiring_photonpt_2016BtoH.root";          prefire_photon_histname_ = "L1prefiring_photonpt_2016BtoH";
  }
  else if (year == 2017) {
    prefire_jet_empt_filename_ = "L1prefiring_jetempt_2017BtoF.root";         prefire_jet_empt_histname_ = "L1prefiring_jetempt_2017BtoF";
    prefire_jet_filename_ = "L1prefiring_jetpt_2017BtoF.root";                prefire_jet_histname_ = "L1prefiring_jetpt_2017BtoF";
    prefire_photon_filename_ = "L1prefiring_photonpt_2017BtoF.root";          prefire_photon_histname_ = "L1prefiring_photonpt_2017BtoF";
  }
  //no prefire issue in 2018
  
  do_prefire_jet_ = (prefire_jet_filename_!="");
  do_prefire_jet_empt_ = (prefire_jet_empt_filename_!="");
  do_prefire_photon_ = (prefire_photon_filename_!="");

  if (do_prefire_jet_) sf_hist_prefire_jet_ = LoadSF<TH2F>(prefire_jet_filename_, prefire_jet_histname_);
  if (do_prefire_jet_empt_) sf_hist_prefire_jet_empt_ = LoadSF<TH2F>(prefire_jet_empt_filename_, prefire_jet_empt_histname_);
  if (do_prefire_photon_) sf_hist_prefire_photon_ = LoadSF<TH2F>(prefire_photon_filename_, prefire_photon_histname_);
}


void PrefireWeighter::EventWeight(nano_tree &nano, float & w_prefire, std::vector<float> & sys_prefire) {
  //prefire weight procedure described at https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe
  //which links here: https://lathomas.web.cern.ch/lathomas/TSGStuff/L1Prefiring/PrefiringMaps_2016and2017/
  
  std::vector<std::pair<double, double> > sfs;
  std::vector<std::pair<double, double> > photon_sfs;

  if (do_prefire_photon_ && (do_prefire_jet_ || do_prefire_jet_empt_)) {

    //get prefiring SFs for isolated photons
    //note: nanoAOD-tools also checks if the photon is in the electron list, we are following TreeMaker, which does not
    for (unsigned int ph_idx = 0; ph_idx < nano.nPhoton(); ph_idx++) {
      //note: TreeMaker seems to use a lower pt cut of 2. GeV while nanoAOD-tools uses a lower cut of 20. GeV
      if (nano.Photon_pt()[ph_idx] < 2. || fabs(nano.Photon_eta()[ph_idx]) < 2.0 || fabs(nano.Photon_eta()[ph_idx]) > 3.0 || photon_overlaps_with_jet[ph_idx]) {
	//insert dummy SFs for unaffected photons to make indexing easier for overlap removal
        photon_sfs.push_back(std::pair<double,double>{1.,0.});
      }
      else {
        photon_sfs.push_back(GetSF(sf_hist_prefire_photon_, nano.Photon_eta()[ph_idx], nano.Photon_pt()[ph_idx], false));
      }
    } // /photon loop

    //get prefiring SFs for jets and photons that overlap with jets
    for (unsigned int jet_idx = 0; jet_idx < nano.nJet(); jet_idx++) {
      if (nano.Jet_pt()[jet_idx] < 2. || fabs(nano.Jet_eta()[jet_idx]) < 2.0 || fabs(nano.Jet_eta()[jet_idx]) > 3.0) continue;
      float pt = nano.Jet_pt()[jet_idx];
      std::pair<double,double> jet_sf{1.,0.};
      if (use_jet_empt_) {
        pt = pt*(nano.Jet_neEmEF()[jet_idx]+nano.Jet_chEmEF()[jet_idx]);
        jet_sf = GetSF(sf_hist_prefire_jet_empt_, nano.Jet_eta()[jet_idx], pt, false);
      }
      else {
        jet_sf = GetSF(sf_hist_prefire_jet_, nano.Jet_eta()[jet_idx], pt, false);
      }
      std::vector<std::pair<double, double>> overlapping_photon_prefiring_sfs;
      bool overlapping_photons = false;
      //photon overlap removal
      for (unsigned int ph_idx = 0; ph_idx < nano.nPhoton(); ph_idx++) {
        if (nano.Photon_pt()[ph_idx] < 2. || fabs(nano.Photon_eta()[ph_idx]) < 2.0 || fabs(nano.Photon_eta()[ph_idx]) > 3.0) continue;
        if (dR(nano.Photon_eta()[ph_idx], nano.Jet_eta()[jet_idx], nano.Photon_phi()[ph_idx], nano.Jet_phi()[jet_idx])<0.4) {
          overlapping_photons = true;
	  if (jet_sf.first < photon_sfs[ph_idx].first) 
            //if jet non-prefire weight is lower, replace photon SF
            photon_sfs[ph_idx] = jet_sf;
        }
      }
      if (!overlapping_photon) {
        //if no overlapping photon directly add jet SF to overall SF
	sfs.push_back(jet_sf);
      }
    } // /jet loop

    //add all overlap-adjusted photon SFs to overall SF, then combine all SFs
    for (unsigned int ph_idx = 0; ph_idx < nano.nPhoton(); ph_idx++) {
      sfs.push_back(photon_sfs[ph_idx]);
    }

    std::pair<double,double> full_sf = accumulate(sfs.cbegin(), sfs.cend(), make_pair(1., 0.), MergeSF);

    w_prefire = full_sf.first;
    sys_prefire = std::vector<float>{static_cast<float>(full_sf.first+full_sf.second),static_cast<float>(full_sf.first-full_sf.second)};
  }
}

