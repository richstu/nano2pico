#include "egamma_weighter.hpp"

#include <string>
#include <numeric>

#include "TFile.h"
#include "TGraphAsymmErrors.h"

#include "utilities.hpp"

#include "correction.h"

using namespace std;

EgammaWeighter::EgammaWeighter(int year, bool preVFP){
  if (year==2016 && preVFP) {
    in_file_ = "data/zgamma/2016preVFP_UL/electron.json";
    in_file_photon_ = "data/zgamma/2016preVFP_UL/photon.json";
    key_ = "2016preVFP";
  } else if (year==2016) {
    in_file_ = "data/zgamma/2016postVFP_UL/electron.json";
    in_file_photon_ = "data/zgamma/2016postVFP_UL/photon.json";
    key_ = "2016postVFP";
  } else if (year==2017) {
    in_file_ = "data/zgamma/2017_UL/electron.json";
    in_file_photon_ = "data/zgamma/2017_UL/photon.json";
    key_ = "2017";
    cout << "OK" << endl;
  } else if (year==2018) {
    in_file_ = "data/zgamma/2018_UL/electron.json";
    in_file_photon_ = "data/zgamma/2018_UL/photon.json";
    key_ = "2018";
  }
}

void EgammaWeighter::GetSF(pico_tree &pico, float &w_lep){
  double sf_tot = 1.0;
  auto cs = correction::CorrectionSet::from_file(in_file_);
  auto Map = cs->at("UL-Electron-ID-SF");
  for(size_t i = 0; i < pico.out_el_sig().size(); ++i){
    if(pico.out_el_sig().at(i) && pico.out_el_pt().at(i) > 20.0){
      auto sf = Map->evaluate({key_, "sf", "wp90iso", std::abs(pico.out_el_eta().at(i)), pico.out_el_pt().at(i)});
      sf_tot *= sf;
    }
  }
  w_lep = sf_tot;
}

void EgammaWeighter::IDSF(pico_tree &pico, float &w_photon_id){
  double sf_tot = 1.0;
  auto cs = correction::CorrectionSet::from_file(in_file_photon_);
  auto Map = cs->at("UL-Photon-ID-SF");
  for(size_t i = 0; i < pico.out_photon_pt().size(); ++i){
    if(pico.out_photon_sig().at(i) && pico.out_photon_pt().at(i) > 20.0){
      auto sf = Map->evaluate({key_, "sf", "wp90", std::abs(pico.out_photon_eta().at(i)), pico.out_photon_pt().at(i)});
      sf_tot *= sf;
    }
  }
  w_photon_id = sf_tot;
}

void EgammaWeighter::CSEVSF(pico_tree &pico, float &w_photon_csev){
  double sf_tot = 1.0;
  string category = "";
  auto cs = correction::CorrectionSet::from_file(in_file_photon_);
  auto Map = cs->at("UL-Photon-CSEV-SF");
  for(size_t i = 0; i < pico.out_photon_pt().size(); ++i){
    if(pico.out_photon_sig().at(i)){
      if (std::abs(pico.out_photon_eta().at(i)) < 1.5){
	category = "EBInc";
	if (std::abs(pico.out_photon_r9().at(i)) > 0.94) {
	  category = "EBHighR9";
	} else {
	  category = "EBLowR9";
	}
      } else {
	category = "EEInc";
        if (std::abs(pico.out_photon_r9().at(i)) > 0.94) {
          category = "EEHighR9";
        } else {
          category = "EELowR9";
        }
      }
      auto sf = Map->evaluate({key_, "sf", "MVA", category});
      sf_tot *= sf;
    }
  }
  w_photon_csev = sf_tot;
}

