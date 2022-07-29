#include "event_weighter.hpp"

#include <string>
#include <numeric>

#include "TFile.h"
#include "TGraphAsymmErrors.h"

#include "utilities.hpp"

using namespace std;

EventWeighter::EventWeighter(int year, bool preVFP){
  if (year==2016 && preVFP) {
    in_file_electron_ = "data/zgamma/2016preVFP_UL/electron.json";
    in_file_photon_   = "data/zgamma/2016preVFP_UL/photon.json";
    in_file_muon_     = "data/zgamma/2016preVFP_UL/muon_Z.json";
    in_file_pu_       = "data/zgamma/2016preVFP_UL/puWeights.json";
    in_file_btag_     = "data/zgamma/2016preVFP_UL/btagging.json";
    key_              = "2016preVFP";
    puName_           = "Collisions16_UltraLegacy_goldenJSON";
  } else if (year==2016) {
    in_file_electron_ = "data/zgamma/2016postVFP_UL/electron.json";
    in_file_photon_   = "data/zgamma/2016postVFP_UL/photon.json";
    in_file_muon_     = "data/zgamma/2016postVFP_UL/muon_Z.json";
    in_file_pu_       = "data/zgamma/2016postVFP_UL/puWeights.json";
    in_file_btag_     = "data/zgamma/2016postVFP_UL/btagging.json";
    key_              = "2016postVFP";
    puName_           = "Collisions16_UltraLegacy_goldenJSON";
  } else if (year==2017) {
    in_file_electron_ = "data/zgamma/2017_UL/electron.json";
    in_file_photon_   = "data/zgamma/2017_UL/photon.json";
    in_file_muon_     = "data/zgamma/2017_UL/muon_Z.json";
    in_file_pu_       = "data/zgamma/2017_UL/puWeights.json";
    in_file_btag_     = "data/zgamma/2017_UL/btagging.json";
    key_              = "2017";
    puName_           = "Collisions17_UltraLegacy_goldenJSON";
  } else if (year==2018) {
    in_file_electron_ = "data/zgamma/2018_UL/electron.json";
    in_file_photon_   = "data/zgamma/2018_UL/photon.json";
    in_file_muon_     = "data/zgamma/2018_UL/muon_Z.json";
    in_file_pu_       = "data/zgamma/2018_UL/puWeights.json";
    in_file_btag_     = "data/zgamma/2018_UL/btagging.json";
    key_              = "2018";
    puName_           = "Collisions18_UltraLegacy_goldenJSON";
  }
  cs_electron_        = correction::CorrectionSet::from_file(in_file_electron_);
  cs_photon_          = correction::CorrectionSet::from_file(in_file_photon_);
  cs_muon_            = correction::CorrectionSet::from_file(in_file_muon_);
  cs_pileup_          = correction::CorrectionSet::from_file(in_file_pu_);
  cs_btag_            = correction::CorrectionSet::from_file(in_file_btag_);
  map_electron_       = cs_electron_->at("UL-Electron-ID-SF");
  map_photon_id_      = cs_photon_->at("UL-Photon-ID-SF");
  map_photon_csev_    = cs_photon_->at("UL-Photon-CSEV-SF");
  map_muon_looseid_   = cs_muon_->at("NUM_LooseID_DEN_genTracks");
  map_muon_highptid_  = cs_muon_->at("NUM_HighPtID_DEN_genTracks");
  map_muon_iso_       = cs_muon_->at("NUM_LooseRelIso_DEN_LooseID");
  map_btag_           = cs_btag_->at("deepCSV_mujets");
  // DeepJet can be used instead of DeepCSV
  // map_btag_           = cs_btag_->at("deepJet_mujets");
  map_pileup_         = cs_pileup_->at(puName_);
}

// Electron MVA ID Scale Factors
void EventWeighter::ElectronIDSF(pico_tree &pico, float &w_el_id){
  double sf_tot = 1.0;
  for(size_t i = 0; i < pico.out_el_sig().size(); ++i){
    if(pico.out_el_sig().at(i) && pico.out_el_pt().at(i) > 10.0){
      auto sf = map_electron_->evaluate({key_, "sf", "wp90iso", std::abs(pico.out_el_eta().at(i)), pico.out_el_pt().at(i)});
      sf_tot *= sf;
    }
  }
  w_el_id = sf_tot;
}

// Photon MVA ID Scale Factors
void EventWeighter::PhotonIDSF(pico_tree &pico, float &w_photon_id){
  double sf_tot = 1.0;
  for(size_t i = 0; i < pico.out_photon_pt().size(); ++i){
    if(pico.out_photon_sig().at(i) && pico.out_photon_pt().at(i) > 20.0){
      auto sf = map_photon_id_->evaluate({key_, "sf", "wp90", std::abs(pico.out_photon_eta().at(i)), pico.out_photon_pt().at(i)});
      sf_tot *= sf;
    }
  }
  w_photon_id = sf_tot;
}

// Photon CSEV Scale Factors
void EventWeighter::PhotonCSEVSF(pico_tree &pico, float &w_photon_csev){
  double sf_tot = 1.0;
  string category = "";
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
      auto sf = map_photon_csev_->evaluate({key_, "sf", "MVA", category});
      sf_tot *= sf;
    }
  }
  w_photon_csev = sf_tot;
}

// Muon ID Scale Factors
void EventWeighter::MuonIDSF(pico_tree &pico, float &w_muon_id){
  double sf_tot = 1.0;
  for(size_t i = 0; i < pico.out_mu_pt().size(); ++i){
    if(pico.out_mu_pt().at(i) > 15.){
      if(pico.out_mu_id().at(i)){
  	auto sf = map_muon_looseid_->evaluate({key_ + "_UL", std::abs(pico.out_mu_eta().at(i)), pico.out_mu_pt().at(i), "sf"});
  	sf_tot *= sf;
      }
      else if(pico.out_mu_highptid().at(i) && pico.out_mu_pt().at(i) > 200.){
  	auto sf = map_muon_highptid_->evaluate({key_ + "_UL", std::abs(pico.out_mu_eta().at(i)), pico.out_mu_pt().at(i), "sf"});
  	sf_tot *= sf;
      }
    }
  }
  w_muon_id = sf_tot;
}

// Muon Iso Scale Factors
void EventWeighter::MuonIsoSF(pico_tree &pico, float &w_muon_iso){
  double sf_tot = 1.0;
  for(size_t i = 0; i < pico.out_mu_pt().size(); ++i){
    if(pico.out_mu_pt().at(i) > 15.){
      if(pico.out_mu_reliso().at(i) < 0.35){
  	auto sf = map_muon_iso_->evaluate({key_ + "_UL", std::abs(pico.out_mu_eta().at(i)), pico.out_mu_pt().at(i), "sf"});
  	sf_tot *= sf;
      }
    }
  }
  w_muon_iso = sf_tot;
}

// Pileup Scale Factors
void EventWeighter::PileupSF(pico_tree &pico, float &w_pu){
  auto sf = map_pileup_->evaluate({float(pico.out_npu_tru_mean()), "nominal"});
  w_pu = sf;
}

// b-tagging Scale Factors
void EventWeighter::bTaggingSF(pico_tree &pico, float &w_btag){
  double sf_tot = 1.0;
  auto n_jets = pico.out_jet_isgood().size();
  for(size_t i = 0; i < n_jets; ++i){
    int hadronFlavour = abs(pico.out_jet_hflavor().at(i));
    if(pico.out_jet_isgood().at(i) && hadronFlavour==5 && std::abs(pico.out_jet_eta().at(i)) < 2.4){
      sf_tot *= map_btag_->evaluate({"central", "M", hadronFlavour, std::abs(pico.out_jet_eta().at(i)), pico.out_jet_pt().at(i)}); // M stands for Medium WP
    }
  }
  w_btag = sf_tot;
}
