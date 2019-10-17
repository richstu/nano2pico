#include "btag_weighter.hpp"

#include <cmath>

#include "utilities.hpp"

using namespace std;

const vector<BTagEntry::OperatingPoint> BTagWeighter::op_pts_{BTagEntry::OP_LOOSE, BTagEntry::OP_MEDIUM, BTagEntry::OP_TIGHT};
const vector<BTagEntry::JetFlavor> BTagWeighter::flavors_{BTagEntry::FLAV_B, BTagEntry::FLAV_C, BTagEntry::FLAV_UDSG};

namespace{
  template<typename T, typename... Args>
    std::unique_ptr<T> MakeUnique(Args&&... args){
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
  }
}

BTagWeighter::BTagWeighter(int year_, bool isFastsim_, bool doDeepFlav_, const vector<float> &btag_wpts): 
  year(year_),
  isFastsim(isFastsim_),
  doDeepFlav(doDeepFlav_),
  wp_loose(btag_wpts[0]),
  wp_medium(btag_wpts[0]),
  wp_tight(btag_wpts[0]){

  // setup SFs and WPs depending on the year
  btag_efficiencies_deep_ = vector<TH3D>(); btag_efficiencies_deep_.resize(op_pts_.size());
  TString beff_file = "";
  if (year==2016) {
    if (doDeepFlav) {
      calib_deep_full_ = unique_ptr<BTagCalibration>(new BTagCalibration("deep_flavor", "data/DeepJet_2016LegacySF_WP_V1.csv"));
      calib_deep_fast_ = unique_ptr<BTagCalibration>(new BTagCalibration("deep_flavor", "data/DeepFlav_13TEV_16SL_18_3_2019.csv"));
      beff_file = "data/btagEfficiency_DeepFlavor_2016.root";
    } else {
      calib_deep_full_ = unique_ptr<BTagCalibration>(new BTagCalibration("csvv2_deep", "data/DeepCSV_2016LegacySF_WP_V1.csv"));
      calib_deep_fast_ = unique_ptr<BTagCalibration>(new BTagCalibration("csvv2_deep", "data/deepcsv_13TEV_16SL_18_3_2019.csv"));
      beff_file = "data/btagEfficiency_DeepCSV_2016.root";
    }
  } else if (year==2017) {  
    if (doDeepFlav) {
      calib_deep_full_ = unique_ptr<BTagCalibration>(new BTagCalibration("deep_flavor", "data/DeepFlavour_94XSF_WP_V3_B_F.csv"));
      calib_deep_fast_ = unique_ptr<BTagCalibration>(new BTagCalibration("deep_flavor", "data/DeepFlav_13TEV_17SL_18_3_2019.csv"));
      beff_file = "data/btagEfficiency_DeepFlavor_2017.root";
    } else {
      calib_deep_full_ = unique_ptr<BTagCalibration>(new BTagCalibration("csvv2_deep", "data/DeepCSV_94XSF_WP_V4_B_F.csv"));
      calib_deep_fast_ = unique_ptr<BTagCalibration>(new BTagCalibration("csvv2_deep", "data/deepcsv_13TEV_17SL_18_3_2019.csv"));
      beff_file = "data/btagEfficiency_DeepCSV_2017.root";
    }
  } else {
    if (doDeepFlav) {
      calib_deep_full_ = unique_ptr<BTagCalibration>(new BTagCalibration("deep_flavor", "data/DeepJet_102XSF_WP_V1.csv"));
      calib_deep_fast_ = unique_ptr<BTagCalibration>(new BTagCalibration("deep_flavor", "data/DeepFlav_13TEV_18SL_7_5_2019.csv"));
      beff_file = "data/btagEfficiency_DeepFlavor_2018.root";
    } else {
      calib_deep_full_ = unique_ptr<BTagCalibration>(new BTagCalibration("csvv2_deep", "data/DeepCSV_102XSF_WP_V1.csv")); 
      calib_deep_fast_ = unique_ptr<BTagCalibration>(new BTagCalibration("csvv2_deep", "data/deepcsv_13TEV_18SL_7_5_2019.csv"));
      beff_file = "data/btagEfficiency_DeepCSV_2018.root";
    }
  }

  TFile file_deep(beff_file, "read");
  for(size_t i = 0; i < op_pts_.size(); ++i){
    const auto op = op_pts_.at(i);
    
    readers_deep_full_[op] = MakeUnique<BTagCalibrationReader>(op, "central", vector<string>{"up", "down"});
    readers_deep_full_.at(op)->load(*calib_deep_full_, BTagEntry::FLAV_UDSG, "incl");
    readers_deep_full_.at(op)->load(*calib_deep_full_, BTagEntry::FLAV_C, "comb");
    readers_deep_full_.at(op)->load(*calib_deep_full_, BTagEntry::FLAV_B, "comb");

    readers_deep_fast_[op] = MakeUnique<BTagCalibrationReader>(op, "central", vector<string>{"up", "down"});
    readers_deep_fast_.at(op)->load(*calib_deep_fast_, BTagEntry::FLAV_UDSG, "fastsim");
    readers_deep_fast_.at(op)->load(*calib_deep_fast_, BTagEntry::FLAV_C, "fastsim");
    readers_deep_fast_.at(op)->load(*calib_deep_fast_, BTagEntry::FLAV_B, "fastsim");

    string hist_deep;
    switch(op){
    case BTagEntry::OP_LOOSE:
      if(doDeepFlav)
        hist_deep = "btagEfficiency_DeepFlavor_loose";
      else
        hist_deep = "btagEfficiency_DeepCSV_loose";
      break;
    case BTagEntry::OP_MEDIUM:
      if(doDeepFlav)
        hist_deep = "btagEfficiency_DeepFlavor_medium";
      else
        hist_deep = "btagEfficiency_DeepCSV_medium";
      break;
    case BTagEntry::OP_TIGHT:
      if(doDeepFlav)
        hist_deep = "btagEfficiency_DeepFlavor_tight";
      else
        hist_deep = "btagEfficiency_DeepCSV_tight";
      break;
    case BTagEntry::OP_RESHAPING:
      hist_deep = "btagEfficiency_DeepCSV_reshaping";
      break;
    default:
      hist_deep = "btagEfficiency";
      break;
    }
    btag_efficiencies_deep_.at(i) = *static_cast<const TH3D*>(file_deep.Get(hist_deep.c_str()));
  }
}

double BTagWeighter::EventWeight(pico_tree &pico, BTagEntry::OperatingPoint op,
                                 const string &bc_full_syst, const string &udsg_full_syst,
                                 const string &bc_fast_syst, const string &udsg_fast_syst) const{
  double product = 1.;
  auto n_jets = pico.out_jet_islep().size();
  for(size_t i = 0; i < n_jets; ++i){
    if(!pico.out_jet_islep().at(i)){
      product *= JetBTagWeight(pico, i, op, bc_full_syst, udsg_full_syst, bc_fast_syst, udsg_fast_syst);
    }
  }
  return product;
}

double BTagWeighter::EventWeight(pico_tree &pico, BTagEntry::OperatingPoint op,
                                 const string &bc_full_syst, const string &udsg_full_syst) const{
  double product = 1.;
  auto n_jets = pico.out_jet_islep().size();
  for(size_t i = 0; i < n_jets; ++i){
    if(!pico.out_jet_islep().at(i)){
      product *= JetBTagWeight(pico, i, op, bc_full_syst, udsg_full_syst);
    }
  }
  return product;
}

double BTagWeighter::EventWeight(pico_tree &pico, const vector<BTagEntry::OperatingPoint> &ops,
                                 const string &bc_full_syst, const string &udsg_full_syst) const{
  double product = 1.;
  auto n_jets = pico.out_jet_islep().size();
  for(size_t i = 0; i < n_jets; ++i){
    if(!pico.out_jet_islep().at(i)){
      product *= JetBTagWeight(pico, i, ops,bc_full_syst, udsg_full_syst);
    }
  }
  return product;
}

double BTagWeighter::EventWeight(pico_tree &pico, const vector<BTagEntry::OperatingPoint> &ops,
                                 const string &bc_full_syst, const string &udsg_full_syst,
                                 const string &bc_fast_syst, const string &udsg_fast_syst) const{
  double product = 1.;
  auto n_jets = pico.out_jet_islep().size();
  for(size_t i = 0; i < n_jets; ++i){
    if(!pico.out_jet_islep().at(i)){
      product *= JetBTagWeight(pico, i, ops, bc_full_syst, udsg_full_syst, bc_fast_syst, udsg_fast_syst);
    }
  }
  return product;
}

double BTagWeighter::JetBTagWeight(pico_tree &pico, size_t ijet, BTagEntry::OperatingPoint op,
                                   const string &bc_full_syst, const string &udsg_full_syst,
                                   const string &bc_fast_syst, const string &udsg_fast_syst) const{
  return JetBTagWeight(pico, ijet, vector<BTagEntry::OperatingPoint>{op},
                       bc_full_syst, udsg_full_syst, bc_fast_syst, udsg_fast_syst);
}

double BTagWeighter::JetBTagWeight(pico_tree &pico, size_t ijet, BTagEntry::OperatingPoint op,
                                   const string &bc_full_syst, const string &udsg_full_syst) const{
  return JetBTagWeight(pico, ijet, vector<BTagEntry::OperatingPoint>{op},
                       bc_full_syst, udsg_full_syst,"central", "central");
}

double BTagWeighter::JetBTagWeight(pico_tree &pico, size_t ijet, const vector<BTagEntry::OperatingPoint> &ops,
                                   const string &bc_full_syst, const string &udsg_full_syst) const{
  return JetBTagWeight(pico, ijet, ops,
                       bc_full_syst, udsg_full_syst,"central", "central");
}

double BTagWeighter::JetBTagWeight(pico_tree &pico, size_t ijet, const vector<BTagEntry::OperatingPoint> &ops,
                                   const string &bc_full_syst, const string &udsg_full_syst,
                                   const string &bc_fast_syst, const string &udsg_fast_syst) const{
  // procedure from https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#1a_Event_reweighting_using_scale
  int hadronFlavour = abs(pico.out_jet_hflavor().at(ijet));
  BTagEntry::JetFlavor flav;
  string full_syst, fast_syst;
  switch(hadronFlavour){
    case 5: flav = BTagEntry::FLAV_B; break;
    case 4: flav = BTagEntry::FLAV_C; break;
    default: flav = BTagEntry::FLAV_UDSG; break;
  }
  switch(flav){
    case BTagEntry::FLAV_B:
    case BTagEntry::FLAV_C:
      full_syst = bc_full_syst;
      fast_syst = bc_fast_syst;
      break;
    case BTagEntry::FLAV_UDSG:
      full_syst = udsg_full_syst;
      fast_syst = udsg_fast_syst;
      break;
    default:
      ERROR("Did not recognize BTagEntry::JetFlavor "+std::to_string(static_cast<int>(flav)));
  }

  vector<float> opcuts;
  for (auto &iop: ops) {
    if (iop==BTagEntry::OP_LOOSE) opcuts.push_back(wp_loose);
    else if (iop==BTagEntry::OP_MEDIUM) opcuts.push_back(wp_medium); 
    else if (iop==BTagEntry::OP_TIGHT) opcuts.push_back(wp_tight);
  }

  float csv = pico.out_jet_deepcsv().at(ijet);

  int tag = -1;
  for (unsigned iop(0); iop<opcuts.size(); iop++) 
    if (csv>opcuts[iop]) tag = iop;

  const map<BTagEntry::OperatingPoint, unique_ptr<BTagCalibrationReader> > *ireaders_full = &readers_deep_full_;

  const map<BTagEntry::OperatingPoint, unique_ptr<BTagCalibrationReader> > *ireaders_fast = &readers_deep_fast_;

  double jet_pt = pico.out_jet_pt().at(ijet);
  double jet_eta = pico.out_jet_eta().at(ijet);
  double eff1(1), eff2(0), sf1(1), sf2(1), sf1_fs(1), sf2_fs(1);
  if (tag >= 0){
    BTagEntry::OperatingPoint iop = ops[tag];
    eff1 = GetMCTagEfficiency(hadronFlavour, jet_pt, jet_eta, iop);
    sf1 = ireaders_full->at(iop)->eval_auto_bounds(full_syst, flav, jet_eta, jet_pt);
    if (isFastsim) sf1_fs = ireaders_fast->at(iop)->eval_auto_bounds(fast_syst, flav, jet_eta, jet_pt);
  }
  if (tag < int(ops.size())-1) {
    BTagEntry::OperatingPoint iop = ops[tag+1];
    eff2 = GetMCTagEfficiency(hadronFlavour, jet_pt, jet_eta, iop);
    sf2 = ireaders_full->at(iop)->eval_auto_bounds(full_syst, flav, jet_eta, jet_pt);
    if (isFastsim) sf2_fs = ireaders_fast->at(iop)->eval_auto_bounds(fast_syst, flav, jet_eta, jet_pt);
  }

  double eff1_fs(eff1/sf1_fs), eff2_fs(eff2/sf2_fs);
  double result = (sf1*sf1_fs*eff1_fs-sf2*sf2_fs*eff2_fs)/(eff1_fs-eff2_fs);
  if(std::isnan(result) || std::isinf(result)){
    result = 1.;
    DBG("SF is NaN or inf. Setting to 1.!");
  }
  return result;
}

double BTagWeighter::GetMCTagEfficiency(int pdgId, float pT, float eta,
                                        BTagEntry::OperatingPoint op) const{
  size_t rdr_idx = distance(op_pts_.cbegin(), find(op_pts_.cbegin(), op_pts_.cend(), op));
  pdgId = abs(pdgId);
  if(pdgId != 4 && pdgId != 5){
    // in the ghost clustering scheme to determine flavor, there are only b, c and other (id=0) flavors
    pdgId = 0;
  }

  int bin;
  float eff;
  bin = btag_efficiencies_deep_.at(rdr_idx).FindFixBin(fabs(eta), pT, pdgId);
  eff = btag_efficiencies_deep_.at(rdr_idx).GetBinContent(bin);
  
  return eff;
}
