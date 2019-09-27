#include "lepton_weighter.hpp"

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

LeptonWeighter::LeptonWeighter(int year){
  if (year==2016) {
    in_full_mu_med_ = "TnP_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root";  hist_full_mu_med_ = "SF";
    in_full_mu_iso_ = "TnP_NUM_MiniIsoTight_DENOM_MediumID_VAR_map_pt_eta.root";  hist_full_mu_iso_ = "SF";
    in_full_mu_trk_ = "TnP_NUM_MediumIP2D_DENOM_LooseID_VAR_map_pt_eta.root";  hist_full_mu_trk_ = "SF";
    in_full_el_med_ = "ElectronScaleFactors_Run2016.root";  hist_full_el_med_ = "Run2016_CutBasedMediumNoIso94XV2";
    in_full_el_iso_ = "ElectronScaleFactors_Run2016.root";  hist_full_el_iso_ = "Run2016_Mini";
    in_full_el_trk_ = "egammaEffi_EGM2D_ETge20_recoSF2016_19_02_09.root";  hist_full_el_trk_ = "EGamma_SF2D";

    in_fast_mu_med_ = "sf_fast_mu_mediumID_2016.root";  hist_fast_mu_med_ = "histo2D";
    in_fast_mu_iso_ = "sf_fast_mu_mini02_2016.root";  hist_fast_mu_iso_ = "histo2D";

    in_fast_el_med_ = "sf_fast_el_mediumID_2016.root";  hist_fast_el_med_ = "histo2D";
    in_fast_el_iso_ = "sf_fast_el_mini01_2016.root";  hist_fast_el_iso_ = "histo2D";
  } else if (year==2017) {
    in_full_mu_med_ = "Muon_Run2017_SF_ID.root";  hist_full_mu_med_ = "NUM_MediumID_DEN_genTracks_pt_abseta";
    in_full_mu_iso_ = "Muon_MinIso02_wrtMediumID_SF_Run2017.root";  hist_full_mu_iso_ = "TnP_MC_NUM_MiniIso02Cut_DEN_MediumID_PAR_pt_eta";
    in_full_mu_trk_ = "";  hist_full_mu_trk_ = ""; // not needed for 2017
    in_full_el_med_ = "ElectronScaleFactors_Run2017.root";  hist_full_el_med_ = "Run2017_CutBasedMediumNoIso94XV2";
    in_full_el_iso_ = "ElectronScaleFactors_Run2017.root";  hist_full_el_iso_ = "Run2017_MVAVLooseTightIP2DMini";
    in_full_el_trk_ = "egammaEffi_EGM2D_ETge20_recoSF2017_19_02_09.root";  hist_full_el_trk_ = "EGamma_SF2D";

    in_fast_mu_med_ = "detailed_mu_full_fast_sf_17.root";  hist_fast_mu_med_ = "miniIso02_MediumId_sf";
    in_fast_mu_iso_ = "";  hist_fast_mu_iso_ = ""; // included in the ID SF above

    in_fast_el_med_ = "detailed_ele_full_fast_sf_17.root";  hist_fast_el_med_ = "CutBasedMediumNoIso94XV2_sf";
    in_fast_el_iso_ = "detailed_ele_full_fast_sf_17.root";  hist_fast_el_iso_ = "MVAVLooseTightIP2DMini_sf";
  } else {
    in_full_mu_med_ = "Muon_Run2018_SF_ID.root";  hist_full_mu_med_ = "NUM_MediumID_DEN_TrackerMuons_pt_abseta";
    in_full_mu_iso_ = "Muon_MinIso02_wrtMediumID_SF_Run2017.root";  hist_full_mu_iso_ = "TnP_MC_NUM_MiniIso02Cut_DEN_MediumID_PAR_pt_eta"; // official recommendation is now to use 2017
    in_full_mu_trk_ = "";  hist_full_mu_trk_ = ""; // not needed for 2018
    in_full_el_med_ = "ElectronScaleFactors_Run2018.root";  hist_full_el_med_ = "Run2018_CutBasedMediumNoIso94XV2";
    in_full_el_iso_ = "ElectronScaleFactors_Run2018.root";  hist_full_el_iso_ = "Run2018_Mini";
    in_full_el_trk_ = "egammaEffi_EGM2D_ETge10_recoSF2018_19_04_04.root";  hist_full_el_trk_ = "EGamma_SF2D";

    in_fast_mu_med_ = "detailed_mu_full_fast_sf_18.root";  hist_fast_mu_med_ = "miniIso02_MediumId_sf";
    in_fast_mu_iso_ = "";  hist_fast_mu_iso_ = ""; // included in the ID SF above

    in_fast_el_med_ = "detailed_ele_full_fast_sf_18.root";  hist_fast_el_med_ = "CutBasedMediumNoIso94XV2_sf";
    in_fast_el_iso_ = "detailed_ele_full_fast_sf_18.root";  hist_fast_el_iso_ = "MVAVLooseTightIP2DMini_sf";
  }

  do_full_el_med_ = (in_full_el_med_!=""); do_full_el_iso_ = (in_full_el_iso_!=""); do_full_el_trk_ = (in_full_el_trk_!="");
  do_full_mu_med_ = (in_full_mu_med_!=""); do_full_mu_iso_ = (in_full_mu_iso_!=""); do_full_mu_trk_ = (in_full_mu_trk_!="");

  do_fast_el_med_ = (in_fast_el_med_!=""); do_fast_el_iso_ = (in_fast_el_iso_!="");
  do_fast_mu_med_ = (in_fast_mu_med_!=""); do_fast_mu_iso_ = (in_fast_mu_iso_!="");

  //https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#Muons_AN1
  if (do_full_mu_med_) sf_full_mu_med_ = LoadSF<TH2F>(in_full_mu_med_, hist_full_mu_med_);
  if (do_full_mu_iso_) sf_full_mu_iso_ = LoadSF<TH2F>(in_full_mu_iso_, hist_full_mu_iso_);
  if (do_full_mu_trk_) sf_full_mu_trk_ = LoadSF<TH2F>(in_full_mu_trk_, hist_full_mu_trk_);

  //https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#Electrons_AN1
  if (do_full_el_med_) sf_full_el_med_ = LoadSF<TH2F>(in_full_el_med_, hist_full_el_med_);
  if (do_full_el_iso_) sf_full_el_iso_ = LoadSF<TH2F>(in_full_el_iso_, hist_full_el_iso_);
  if (do_full_el_trk_) sf_full_el_trk_ = LoadSF<TH2F>(in_full_el_trk_, hist_full_el_trk_);

  if (do_fast_mu_med_) sf_fast_mu_med_ = LoadSF<TH2D>(in_fast_mu_med_, hist_fast_mu_med_);
  if (do_fast_mu_iso_) sf_fast_mu_iso_ = LoadSF<TH2D>(in_fast_mu_iso_, hist_fast_mu_iso_);

  if (do_fast_el_med_) sf_fast_el_med_ = LoadSF<TH2D>(in_fast_el_med_, hist_fast_el_med_);
  if (do_fast_el_iso_) sf_fast_el_iso_ = LoadSF<TH2D>(in_fast_el_iso_, hist_fast_el_iso_);
}

void LeptonWeighter::FullSim(pico_tree &pico, float &w_lep, vector<float> &sys_lep){
  pair<double, double> sf(1., 0.);
  for(size_t i = 0; i < pico.out_mu_sig().size(); ++i){
    if(pico.out_mu_sig().at(i)){
      sf = MergeSF(sf, GetMuonScaleFactor(pico, i));
    }
  }
  for(size_t i = 0; i < pico.out_el_sig().size(); ++i){
    if(pico.out_el_sig().at(i)){
      sf = MergeSF(sf, GetElectronScaleFactor(pico, i));
    }
  }
  w_lep = sf.first;
  sys_lep = vector<float>{static_cast<float>(sf.first+sf.second),
                          static_cast<float>(sf.first-sf.second)};
}

void LeptonWeighter::FastSim(pico_tree &pico, float &w_fs_lep, vector<float> &sys_fs_lep){
  pair<double, double> sf(1., 0.);
  for(size_t i = 0; i < pico.out_mu_sig().size(); ++i){
    if(pico.out_mu_sig().at(i)){
      sf = MergeSF(sf, GetMuonScaleFactorFS(pico, i));
    }
  }
  for(size_t i = 0; i < pico.out_el_sig().size(); ++i){
    if(pico.out_el_sig().at(i)){
      sf = MergeSF(sf, GetElectronScaleFactorFS(pico, i));
    }
  }
  w_fs_lep = sf.first;
  sys_fs_lep = vector<float>{static_cast<float>(sf.first+sf.second),
                             static_cast<float>(sf.first-sf.second)};
}

std::pair<double, double> LeptonWeighter::GetMuonScaleFactor(pico_tree &pico, size_t imu){
  //https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#Data_leading_order_FullSim_MC_co
  //ID, iso, tracking SFs applied
  //No stat error, 3% systematic from ID, iso
  double pt = pico.out_mu_pt().at(imu);
  double eta = pico.out_mu_eta().at(imu);
  double abseta = fabs(eta);
  vector<pair<double, double> > sfs;
  if (do_full_mu_med_) {
    sfs.push_back(GetSF(sf_full_mu_med_, pt, abseta, false));
    sfs.push_back(make_pair(1., 0.03));//Systematic uncertainty
  } 
  if (do_full_mu_iso_) {
    sfs.push_back(GetSF(sf_full_mu_iso_, pt, abseta, false));
    sfs.push_back(make_pair(1., 0.03));//Systematic uncertainty
  }
  if (do_full_mu_trk_) {
    sfs.push_back(GetSF(sf_full_mu_trk_, pt, abseta, false));
    sfs.push_back(make_pair(1., 0.03));//Systematic uncertainty
  }
  return accumulate(sfs.cbegin(), sfs.cend(), make_pair(1., 0.), MergeSF);
}

std::pair<double, double> LeptonWeighter::GetElectronScaleFactor(pico_tree &pico, size_t iel){
  //https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#Data_leading_order_FullSim_M_AN1
  //ID, iso, tracking SFs applied
  //ID iso systematics built-in
  //Tracking SFs from https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale
  //3% tracking systematic below 20 GeV
  double pt = pico.out_el_scpt().at(iel);
  double eta = pico.out_el_sceta().at(iel);
  // double abseta = fabs(eta);
  vector<pair<double, double> > sfs;
  //Axes swapped, asymmetric in eta
  if (do_full_el_med_) sfs.push_back(GetSF(sf_full_el_med_, eta, pt));
  if (do_full_el_iso_) sfs.push_back(GetSF(sf_full_el_iso_, eta, pt));
  if (do_full_el_trk_) sfs.push_back(GetSF(sf_full_el_trk_, eta, pt));
  sfs.push_back(make_pair(1., pt<20. || pt >80. ? 0.01 : 0.));//Systematic uncertainty
  return accumulate(sfs.cbegin(), sfs.cend(), make_pair(1., 0.), MergeSF);
}

std::pair<double, double> LeptonWeighter::GetMuonScaleFactorFS(pico_tree &pico, size_t imu){
  //https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#FullSim_FastSim_TTBar_MC_compari
  //ID, iso SFs applied
  //No stat error, 2% systematic from ID, iso
  double pt = pico.out_mu_pt().at(imu);
  double abseta = fabs(pico.out_mu_eta().at(imu));
  vector<pair<double, double> > sfs;
  if (do_fast_mu_med_) {
    sfs.push_back(GetSF(sf_fast_mu_med_, pt, abseta, false));
    sfs.push_back(make_pair(1., 0.02));
  }
  if (do_fast_mu_iso_) {
    sfs.push_back(GetSF(sf_fast_mu_iso_, pt, abseta, false));
    sfs.push_back(make_pair(1., 0.02));
  }
  return accumulate(sfs.cbegin(), sfs.cend(), make_pair(1., 0.), MergeSF);
}

std::pair<double, double> LeptonWeighter::GetElectronScaleFactorFS(pico_tree &pico, size_t iel){
  //https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#FullSim_FastSim_TTBar_MC_com_AN1
  //ID, iso SFs applied
  //No stat error, 2% systematic from ID, iso
  double pt = pico.out_el_scpt().at(iel);
  double abseta = fabs(pico.out_el_sceta().at(iel));
  vector<pair<double, double> > sfs;
  if (do_fast_el_med_) {
    sfs.push_back(GetSF(sf_fast_el_med_, pt, abseta, false));
    sfs.push_back(make_pair(1., 0.02));//Systematic uncertainty
  }
  if (do_fast_el_iso_) {
    sfs.push_back(GetSF(sf_fast_el_iso_, pt, abseta, false));
    sfs.push_back(make_pair(1., 0.02));//Systematic uncertainty
  }
  return accumulate(sfs.cbegin(), sfs.cend(), make_pair(1., 0.), MergeSF);
}

