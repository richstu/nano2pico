#include <ctime>

#include <iostream>
#include <iomanip>
#include <bitset>
#include <regex>

#include <getopt.h>

#include "TError.h"

#include "nano_tree.hpp"
#include "pico_tree.hpp"
#include "corrections_tree.hpp"
#include "utilities.hpp"
#include "cross_sections.hpp"

#include "mc_producer.hpp"
#include "el_producer.hpp"
#include "mu_producer.hpp"
#include "dilep_producer.hpp"
#include "tk_producer.hpp"
#include "photon_producer.hpp"
#include "jetmet_producer.hpp"
#include "hig_producer.hpp"
#include "zgamma_producer.hpp"
#include "gammagamma_producer.hpp"
#include "bb_producer.hpp"
#include "bbgammagamma_producer.hpp"
#include "in_json.hpp"

#include "btag_weighter.hpp"
#include "lepton_weighter.hpp"
#include "prefire_weighter.hpp"
#include "photon_weighter.hpp"
#include "event_tools.hpp"
#include "isr_tools.hpp"
#include "event_weighter.hpp"
#include "trigger_weighter.hpp"

using namespace std;

namespace {
  string in_file = "";
  string in_dir = "";
  string out_dir = "";
  int nent_test = -1;
  bool debug = false;
  // requirements for jets to be counted in njet, mofified for Zgamma below
  float min_jet_pt = 30.0;
  float max_jet_eta =  2.4;
}

void WriteDataQualityFilters(nano_tree& nano, pico_tree& pico);
void CopyTriggerDecisions(nano_tree& nano, pico_tree& pico);
void Initialize(corrections_tree& wgt_sums);
void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  // gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches       
  GetOptions(argc, argv);

  if(in_file=="" || in_dir=="" || out_dir == "") {
    cout<<"ERROR: Input file, sum-of-weights and/or output directory not specified. Exit."<<endl;
    exit(1);
  }

  //bool isData = Contains(in_file, "Run201") ? true : false;
  bool isData = Contains(in_file, "Run20") ? true : false; //Changed to allow for Run 3 data
  bool isFastsim = Contains(in_file, "Fast") ? true : false;
  bool isSignal = Contains(in_file, "TChiHH") || Contains(in_file, "T5qqqqZH") ? true : false;
  bool isZgamma = Contains(out_dir, "zgamma");
  bool isHiggsino = Contains(out_dir, "higgsino");
  int year = -1;
  int isAPV = false;
  int is_preUL = true;
  if (regex_search(in_file, std::regex("RunIISummer\\d\\dUL"))) is_preUL = false;
  if (regex_search(in_file, std::regex("UL201\\d"))) is_preUL = false;
  // Find year and isAPV for MC
  if (!isData) { // MC
    if (!is_preUL) { // UL
      if (regex_search(in_file, std::regex("RunIISummer\\d\\dUL16NanoAODAPV"))) isAPV = true;
      if (regex_search(in_file, std::regex("RunIISummer\\d\\dUL16"))) year = 2016;
      else if (regex_search(in_file, std::regex("RunIISummer\\d\\dUL17"))) year = 2017;
      else if (regex_search(in_file, std::regex("RunIISummer\\d\\dUL18"))) year = 2018;
    } else { // Not UL
      if (regex_search(in_file, std::regex("RunIISummer16"))) year = 2016;
      else if (regex_search(in_file, std::regex("RunIIFall17"))) year = 2017;
      else if (regex_search(in_file, std::regex("RunIIAutumn18"))) year = 2018;
      else if (regex_search(in_file, std::regex("Run3Summer22"))) year = 2022;
      else if (regex_search(in_file, std::regex("Run3Summer23"))) year = 2023;
    }
  } else { // Data
    if (Contains(in_file, "HIPM")) isAPV = true;
    if (Contains(in_file, "Run2016")) year = 2016;
    else if (Contains(in_file, "Run2017")) year = 2017;
    else if (Contains(in_file, "Run2018")) year = 2018;
    else if (Contains(in_file, "Run2022")) year = 2022;
    else if (Contains(in_file, "Run2023")) year = 2023;
  }
  if (year < 0) {
    cout<<"ERROR: Add code for new year!"<<endl;
    exit(1);
  }

  bool is2022preEE = false; //Classify data and MC into pre and post EE for 2022
  if(year == 2022){ 
    if(isData){
      if (Contains(in_file, "2022C") || Contains(in_file, "2022D")){
        is2022preEE = true;
      }
    } else {
      if (!Contains(in_file, "Summer2022EE")){
        is2022preEE = true;
      }
    }
  }

  bool is2023preBPix = false;
  if(year == 2023){ 
    if(isData){
      if (Contains(in_file, "2023B") || Contains(in_file, "2023C")){
        is2023preBPix = true;
      }
    } else {
      if (!Contains(in_file, "Summer23BPix")){
        is2023preBPix = true;
      }
    }
  }

  string year_string;
  if (year == 2016 && isAPV)               year_string = "2016APV";
  else if (year == 2016 && !isAPV)         year_string = "2016";
  else if (year == 2017)                   year_string = "2017";
  else if (year == 2018)                   year_string = "2018";
  else if (year == 2022 && is2022preEE)    year_string = "2022";
  else if (year == 2022 && !is2022preEE)   year_string = "2022EE";
  else if (year == 2023 && is2023preBPix)  year_string = "2023";
  else if (year == 2023 && !is2023preBPix) year_string = "2023BPix";
  else {
    cout << "ERROR: unknown year";
    exit(1);
  }

  //if (Contains(in_file, "RunIISummer20")) { 
  //  is_preUL = false;
  //  if (Contains(in_file, "RunIISummer20UL16NanoAODAPV")) isAPV = true;
  //  if (Contains(in_file, "RunIISummer20UL16")) year = 2016;
  //  else if (Contains(in_file, "RunIISummer20UL17")) year = 2017;
  //  else year = 2018;
  //} else if (Contains(in_file, "RunIISummer19")) { 
  //  is_preUL = false;
  //  if (Contains(in_file, "RunIISummer19UL16NanoAODAPV")) isAPV = true;
  //  if (Contains(in_file, "RunIISummer19UL16")) year = 2016;
  //  else if (Contains(in_file, "RunIISummer19UL17")) year = 2017;
  //  else year = 2018;
  //} else if (Contains(in_file, "Run3Summer22")){
  //  is_preUL = false;
  //  if (Contains(in_file, "Run3_2022")){
  //    year = 2022;
  //    cout<<"Using 2018 btag wpts by default currently."<<endl;
  //  } 
  //  else cout<<"Add code for new year!"<<endl;
  //} else {
  //  year = Contains(in_file, "RunIISummer16") ? 2016 : (Contains(in_file, "RunIIFall17") ? 2017 : 2018);
  //}
  //if (isData) {
  //  year = Contains(in_file, "Run2016") ? 2016 : (Contains(in_file, "Run2017") ? 2017 : (Contains(in_file, "Run2018") ? 2018: 2022));
  //}

  vector<vector<int>> VVRunLumi;
  if (isData) {
    switch (year) {
      case 2016:
        if (Contains(in_file, "UL2016")) VVRunLumi = MakeVRunLumi("goldenUL2016");
        else VVRunLumi = MakeVRunLumi("golden2016");
        break;
      case 2017:
        if (Contains(in_file, "UL2017")) VVRunLumi = MakeVRunLumi("goldenUL2017");
        else VVRunLumi = MakeVRunLumi("golden2017");
        break;
      case 2018:
        if (Contains(in_file, "UL2018")) VVRunLumi = MakeVRunLumi("goldenUL2018");
        else VVRunLumi = MakeVRunLumi("golden2018");
        break;
      case 2022:
        if (Contains(in_file, "2022")) VVRunLumi = MakeVRunLumi("golden2022");
        break;
      case 2023:
        if (Contains(in_file, "2023")) VVRunLumi = MakeVRunLumi("golden2023");
        break;
      default:
        cout << "ERROR: no golden cert for given year" << endl;
        exit(1);
    }
  }

  string in_path = in_dir+"/"+in_file;
  string wgt_sums_path = out_dir+"/wgt_sums/wgt_sums_"+in_file;
  string out_path;
  out_path = out_dir+"/raw_pico/raw_pico_"+in_file;

  // Find nanoAOD version
  float nanoaod_version = -1;
  std::smatch nanoad_version_matches;
  bool version_found = std::regex_search(in_file, nanoad_version_matches, std::regex("NanoAOD(?:APVv|v)(\\d+p\\d+|\\d+)"));
  if (version_found) nanoaod_version = std::stof(std::regex_replace(nanoad_version_matches[1].str(), std::regex("p"), "."));
  else {
    bool is_nanoAODv7_found = std::regex_search(in_file, nanoad_version_matches, std::regex("02Apr2020"));
    if (is_nanoAODv7_found) nanoaod_version = 7;
  }
  if (Contains(in_dir, "NanoAODv9UCSB")) nanoaod_version = 9.5;
  if (Contains(in_dir, "NanoAODv12")) nanoaod_version = 12;
  cout<<"Using NanoAOD version: "<<nanoaod_version<<endl;

  time_t begtime, endtime;
  time(&begtime);

  // jet requirements
  if(isZgamma) max_jet_eta  =  4.7;

  // B-tag working points
  // Updated Values May-28-2024 from https://btv-wiki.docs.cern.ch/ScaleFactors/
  // btag_df: WPs for deepJet (DeepFlavourB)
  map<string, vector<float>> btag_df_wpts{
    {"2016APV", vector<float>({0.0508, 0.2598, 0.6502})},
    {"2016", vector<float>({0.0480, 0.2489, 0.6377})},
    {"2017", vector<float>({0.0532, 0.3040, 0.7476})},
    {"2018", vector<float>({0.0490, 0.2783, 0.7100})},
    {"2022", vector<float>({0.0583, 0.3086, 0.7183})},
    {"2022EE", vector<float>({0.0614, 0.3196, 0.73})},
    {"2023", vector<float>({0.0479, 0.2431, 0.6553})},
    {"2023BPix", vector<float>({0.048, 0.2435, 0.6563})}
  };
  // WPs for Run 3 values are for PNet, Run 2 values are for deepCSV (DeepB)
  map<string, vector<float>> btag_wpts{
    {"2016APV", vector<float>({0.2027, 0.6001, 0.8819})},
    {"2016", vector<float>({0.1918, 0.5847, 0.8767})},
    {"2017", vector<float>({0.1355, 0.4506, 0.7738})},
    {"2018", vector<float>({0.1208, 0.4168, 0.7665})},
    {"2022", vector<float>({0.047,  0.245,  0.6734})},
    {"2022EE", vector<float>({0.0499, 0.2605, 0.6915})},  
    {"2023", vector<float>({0.0358, 0.1917, 0.6172})},
    {"2023BPix", vector<float>({0.0359, 0.1919, 0.6133})}
  };

  // Rochester corrections
  string rocco_file = "data/zgamma/2018_UL/RoccoR2018UL.txt";
  if (is_preUL) {
    if (year==2016)
      rocco_file = "data/RoccoR2016.txt";
    else if (year==2017)
      rocco_file = "data/RoccoR2017.txt";
    else if (year==2018)
      rocco_file = "data/RoccoR2018.txt";
    else
      cout<<"WARNING: No rochester corrections for year."<<endl;
  }
  else {
    if (year==2016 && isAPV)
      rocco_file = "data/zgamma/2016preVFP_UL/RoccoR2016aUL.txt";
    else if (year==2016 && !isAPV)
      rocco_file = "data/zgamma/2016postVFP_UL/RoccoR2016bUL.txt";
    else if (year==2017)
      rocco_file = "data/zgamma/2017_UL/RoccoR2017UL.txt";
    else if (year==2018)
      rocco_file = "data/zgamma/2018_UL/RoccoR2018UL.txt";
    else
      cout<<"WARNING: No rochester corrections for year."<<endl;
  }

  //Initialize object producers
  GenParticleProducer mc_producer(year, nanoaod_version);
  ElectronProducer el_producer(year_string, isData, nanoaod_version);
  MuonProducer mu_producer(year, isData, nanoaod_version, rocco_file);
  DileptonProducer dilep_producer(year);
  IsoTrackProducer tk_producer(year);
  PhotonProducer photon_producer(year_string, isData, nanoaod_version);
  JetMetProducer jetmet_producer(year, year_string, nanoaod_version, min_jet_pt, max_jet_eta, 
                                 isData, is_preUL);
  HigVarProducer hig_producer(year);
  ZGammaVarProducer zgamma_producer(year);
  GammaGammaVarProducer gammagamma_producer(year);
  BBVarProducer bb_producer(year);
  BBGammaGammaVarProducer bbgammagamma_producer(year);

  //Initialize scale factor tools
  const string ctr = "central";
  const vector<string> updn = {"up","down"};
  PrefireWeighter prefire_weighter(year, true);
  // Pre-UL scale factors
  const vector<BTagEntry::OperatingPoint> op_all = {BTagEntry::OP_LOOSE, BTagEntry::OP_MEDIUM, BTagEntry::OP_TIGHT};
  BTagWeighter btag_weighter(year, isFastsim, false, btag_wpts[year_string]);
  BTagWeighter btag_df_weighter(year, isFastsim, true, btag_df_wpts[year_string]);
  LeptonWeighter lep_weighter(year, isZgamma);
  LeptonWeighter lep_weighter16gh(year, isZgamma, true);
  PhotonWeighter photon_weighter(year, isZgamma || isHiggsino);
  // UL scale factors
  EventWeighter event_weighter(year_string, btag_df_wpts[year_string]);
  TriggerWeighter trigger_weighter(year_string);
  //cout<<"Is APV: "<<isAPV<<endl;

  // Other tools
  EventTools event_tools(in_path, year, isData, nanoaod_version);
  int event_type = event_tools.GetEventType();
  bool isDY = event_type/100 == 62 ? isZgamma : false;

  ISRTools isr_tools(in_path, year, nanoaod_version, isData);

  // Initialize trees
  nano_tree nano(in_path, nanoaod_version);
  //nano_tree nano(in_path, 9);
  size_t nentries(nent_test>0 ? nent_test : nano.GetEntries());
  cout << "Nano file: " << in_path << endl;
  cout << "Input number of events: " << nentries << endl;
  if (nent_test > nano.GetEntries()) {
    cout << "ERROR: nent: " << nent_test << " is larger than nano.GetEntries(): "<< nano.GetEntries() << endl;
    exit(1);
  }
  // cout << "Running on "<< (isFastsim ? "FastSim" : "FullSim") << endl;
  // cout << "Calculating weights based on " << year << " scale factors." << endl;

  pico_tree pico("", out_path);
  cout << "Writing output to: " << out_path << endl;

  corrections_tree wgt_sums("", wgt_sums_path);
  cout << "Writing sum-of-weights to: " << wgt_sums_path << endl;
  Initialize(wgt_sums);
  wgt_sums.out_nent() = nentries;

  for(size_t entry(0); entry<nentries; ++entry){
    if (debug) cout << "GetEntry: " << entry <<" event = "<<pico.out_event()<< endl;
    nano.GetEntry(entry);
    if (entry%2000==0 || entry == nentries-1) {
      cout<<"Processing event: "<<entry<<endl;
    }

    //skip events that are data but not in the golden json
    if (isData) {
      if(!inJSON(VVRunLumi, nano.run(), nano.luminosityBlock())) continue; 
    }

    bool passed_trig = event_tools.SaveTriggerDecisions(nano, pico, isZgamma);
    if (isData && !passed_trig) {
      continue;
    }
    // event info
    pico.out_event()     = nano.event();
    pico.out_lumiblock() = nano.luminosityBlock();
    pico.out_run()       = nano.run();
    pico.out_type()      = event_type;

    // number of reconstructed primary vertices
    pico.out_npv() = nano.PV_npvs();
    pico.out_npv_good() = nano.PV_npvsGood();
    // number of pileup in mc
    if (!isData) {
      pico.out_npu_tru() = nano.Pileup_nPU();
      pico.out_npu_tru_mean() = nano.Pileup_nTrueInt();
    }

    //pileup energy density
    if (nanoaod_version >= 11 || nanoaod_version == 9.5)
      pico.out_rho() = nano.fixedGridRhoAll();

    // ----------------------------------------------------------------------------------------------
    //            *** Writing physics objects ***
    // N.B. Order in which producers are called matters! E.g. jets are not counted if overlapping 
    // with signal lepton, thus jets must be processed only after leptons have been selected.
    //-----------------------------------------------------------------------------------------------
    if (debug) cout<<"INFO:: Writing leptons, photons and tracks"<<endl;
    vector<int> jet_islep_nano_idx = vector<int>();
    vector<int> jet_isvlep_nano_idx = vector<int>();
    pico.out_nlep() = 0; pico.out_nvlep() = 0; // filled by lepton producers
    vector<int> sig_el_pico_idx = vector<int>();
    vector<int> sig_mu_pico_idx = vector<int>();
    vector<int> photon_el_pico_idx = vector<int>();
    
    vector<int> sig_el_nano_idx = el_producer.WriteElectrons(nano, pico, jet_islep_nano_idx, jet_isvlep_nano_idx, sig_el_pico_idx, photon_el_pico_idx, isZgamma, isFastsim);
    vector<int> sig_mu_nano_idx = mu_producer.WriteMuons(nano, pico, jet_islep_nano_idx, jet_isvlep_nano_idx, sig_mu_pico_idx, isZgamma, isFastsim);
    // save a separate vector with just signal leptons ordered by pt
    struct SignalLepton{ float pt; float eta; float phi; int pdgid;};
    vector<SignalLepton> sig_leps;
    for (auto &iel: sig_el_nano_idx)
      sig_leps.push_back({nano.Electron_pt()[iel], nano.Electron_eta()[iel], 
                          nano.Electron_phi()[iel], nano.Electron_pdgId()[iel]});
    for (auto &imu: sig_mu_nano_idx) 
      sig_leps.push_back({nano.Muon_pt()[imu], nano.Muon_eta()[imu], 
                          nano.Muon_phi()[imu], nano.Muon_pdgId()[imu]});

    auto greaterPt = [](SignalLepton lep1, SignalLepton lep2){ return lep1.pt > lep2.pt;};
    sort(sig_leps.begin(), sig_leps.end(), greaterPt);
    for(auto &ilep : sig_leps) {
      pico.out_lep_pt().push_back(ilep.pt);
      pico.out_lep_eta().push_back(ilep.eta);
      pico.out_lep_phi().push_back(ilep.phi);
      pico.out_lep_pdgid().push_back(ilep.pdgid);
    }
    vector<int> jet_isphoton_nano_idx = vector<int>();
    if(isZgamma || isHiggsino) 
      vector<int> sig_ph_nano_idx = photon_producer.WritePhotons(nano, pico, jet_isphoton_nano_idx,
                                                                 sig_el_nano_idx, sig_mu_nano_idx,
                                                                 photon_el_pico_idx);
    event_tools.WriteStitch(nano, pico);
    tk_producer.WriteIsoTracks(nano, pico, sig_el_nano_idx, sig_mu_nano_idx, isFastsim, is_preUL);
    dilep_producer.WriteDileptons(pico, sig_el_pico_idx, sig_mu_pico_idx);

    if (debug) cout<<"INFO:: Writing gen particles"<<endl;

    if (!isData) mc_producer.WriteGenParticles(nano, pico, isDY);
    isr_tools.WriteISRSystemPt(nano, pico);

    if (debug) cout<<"INFO:: Writing jets, MET and ISR vars"<<endl;

    vector<HiggsConstructionVariables> sys_higvars;
    vector<int> sig_jet_nano_idx = jetmet_producer.WriteJetMet(nano, pico, 
        jet_islep_nano_idx, jet_isvlep_nano_idx, jet_isphoton_nano_idx,
        btag_wpts[year_string], btag_df_wpts[year_string], isFastsim, isSignal, 
        is2022preEE, sys_higvars);
    jetmet_producer.WriteJetSystemPt(nano, pico, sig_jet_nano_idx, btag_wpts[year_string][1], isFastsim); // usually w.r.t. medium WP
    jetmet_producer.WriteFatJets(nano, pico); // jetmet_producer.SetVerbose(nano.nSubJet()>0);
    jetmet_producer.WriteSubJets(nano, pico);
    isr_tools.WriteISRJetMultiplicity(nano, pico);

    // calculate mT only for single lepton events
    pico.out_mt() = -999; 
    if (pico.out_nlep()==1) {
      float MET_pt, MET_phi;
      getMETWithJEC(nano, year, isFastsim, MET_pt, MET_phi, is_preUL);
      if (sig_el_nano_idx.size()>0) {
        pico.out_mt() = GetMT(MET_pt, MET_phi, 
          nano.Electron_pt()[sig_el_nano_idx[0]], nano.Electron_phi()[sig_el_nano_idx[0]]);
      } else {
        pico.out_mt() = GetMT(MET_pt, MET_phi, 
          nano.Muon_pt()[sig_mu_nano_idx[0]], nano.Muon_phi()[sig_mu_nano_idx[0]]);
      }
    } 

    if (pico.out_ntrulep()==1) {
      for (unsigned imc(0); imc<pico.out_mc_id().size(); imc++){
        if (abs(pico.out_mc_id()[imc])==11 || abs(pico.out_mc_id()[imc])==13) {
          // statusflag: 12 = isFirstCopy
          bitset<15> mc_statusFlags(pico.out_mc_statusflag().at(imc));
          if  ((mc_statusFlags[12]==1)) {
            pico.out_mt_tru() = GetMT(pico.out_met_tru(), pico.out_met_tru_phi(), 
                                  pico.out_mc_pt()[imc], pico.out_mc_phi()[imc]);
            break;
          }
        }
      }
    } 
    if (debug) cout<<"INFO:: Writing analysis specific variables"<<endl;
    // might need as input sig_el_nano_idx, sig_mu_nano_idx, sig_ph_nano_idx
    if(isZgamma)
      zgamma_producer.WriteZGammaVars(nano, pico, sig_jet_nano_idx);
  
    if (isHiggsino) gammagamma_producer.WriteGammaGammaVars(pico);
    if (isHiggsino) bb_producer.WriteBBVars(pico, /*doDeepFlav*/false);
    if (isHiggsino) bb_producer.WriteBBVars(pico, /*doDeepFlav*/true);
    if (isHiggsino) bbgammagamma_producer.WriteBBGammaGammaVars(pico);

    //save higgs variables using DeepCSV and DeepFlavor
    hig_producer.WriteHigVars(pico, false, isSignal, sys_higvars, nanoaod_version);
    hig_producer.WriteHigVars(pico, true, isSignal, sys_higvars, nanoaod_version);

    if (debug) cout<<"INFO:: Writing filters and triggers"<<endl;
    // N.B. Jets: pico.out_pass_jets() and pico.out_pass_fsjets() filled in jetmet_producer
    event_tools.WriteDataQualityFilters(nano, pico, sig_jet_nano_idx, min_jet_pt, isFastsim, is_preUL);

    if (isHiggsino) event_tools.WriteTriggerEfficiency(pico);
    if (isZgamma && !isData) {
      std::vector<float> zgamma_trigsfs = trigger_weighter.GetSF(pico);
      pico.out_w_trig() = zgamma_trigsfs[0];
      pico.out_sys_trig().resize(2,0.);
      pico.out_sys_trig()[0] = zgamma_trigsfs[1];
      pico.out_sys_trig()[1] = zgamma_trigsfs[2];
    }

    // ----------------------------------------------------------------------------------------------
    //              *** Calculating weight branches ***
    // ----------------------------------------------------------------------------------------------
    if (debug) cout<<"INFO:: Calculating weights"<<endl;
    float w_lep(1.), w_fs_lep(1.);
    float w_photon(1.);
    vector<float> sys_lep(2,1.), sys_fs_lep(2,1.);
    vector<float> sys_photon(2,1.);

    if (isData) {
      pico.out_w_btag()    = 1.; 
      pico.out_w_btag_df() = 1.; 
      pico.out_w_bhig()    = 1.; 
      pico.out_w_bhig_df() = 1.; 
      pico.out_sys_bchig().resize(2,0); pico.out_sys_udsghig().resize(2,0);
      pico.out_sys_fs_bchig().resize(2,0); pico.out_sys_fs_udsghig().resize(2,0);
      pico.out_w_lep() = 1.;
      pico.out_w_fs_lep() = 1.;
      pico.out_sys_lep().resize(2,0); pico.out_sys_fs_lep().resize(2,0);
      pico.out_w_pu() = 1.;
      pico.out_sys_pu().resize(2, 0);
      pico.out_w_photon() = 1.;
      pico.out_w_trig() = 1.;
      pico.out_sys_photon().resize(2,0);
    } else { // MC
      if ((!is_preUL) || year>=2022) { //UL or run 3
        // ElectronISO SF need to be implemented for non-HToZgamma
        event_weighter.ElectronSF(pico);
        event_weighter.MuonSF(pico);
        event_weighter.PileupSF(pico);
        event_weighter.bTaggingSF(pico);
        event_weighter.PhotonSF(pico);
        pico.out_sys_lep().resize(2,1.); 
        pico.out_sys_photon().resize(2, 1.); 
        pico.out_sys_prefire().resize(2, 1.); 
        pico.out_w_lep()          = pico.out_w_el() * pico.out_w_mu();
        pico.out_sys_lep()[0]     = pico.out_sys_el()[0]*pico.out_sys_mu()[0]; 
        pico.out_sys_lep()[1]     = pico.out_sys_el()[1]*pico.out_sys_mu()[1]; 
        pico.out_sys_fs_bchig().resize(2,1.); 
        pico.out_sys_fs_udsghig().resize(2,1.); 
        pico.out_sys_fs_lep().resize(2,1.);
        pico.out_w_btag()    = 1.; 
        pico.out_w_btag_df() = 1.; 
        pico.out_w_bhig()    = 1.; 
        pico.out_w_fs_lep()  = 1.;
        if (year >= 2022) {
          pico.out_w_prefire()      = 1.0;
          pico.out_sys_prefire()[0] = 1.0;
          pico.out_sys_prefire()[1] = 1.0;
        }
        else {
          pico.out_w_prefire()      = nano.L1PreFiringWeight_Nom();
          pico.out_sys_prefire()[0] = nano.L1PreFiringWeight_Up();
          pico.out_sys_prefire()[1] = nano.L1PreFiringWeight_Dn();
        }
      } else { // Pre-UL run 2
        pico.out_w_btag()    = btag_weighter.EventWeight(pico, BTagEntry::OP_MEDIUM, ctr, ctr);; 
        pico.out_w_btag_df() = btag_df_weighter.EventWeight(pico, BTagEntry::OP_MEDIUM, ctr, ctr); 
        pico.out_w_bhig()    = btag_weighter.EventWeight(pico, op_all, ctr, ctr); 
        pico.out_w_bhig_df() = btag_df_weighter.EventWeight(pico, op_all, ctr, ctr); 
        pico.out_sys_bchig().resize(2,0); pico.out_sys_udsghig().resize(2,0);
        pico.out_sys_fs_bchig().resize(2,0); pico.out_sys_fs_udsghig().resize(2,0);
        for(size_t i = 0; i<2; ++i){ 
          pico.out_sys_bchig()[i]   = btag_weighter.EventWeight(pico, op_all, updn[i], ctr);
          pico.out_sys_udsghig()[i] = btag_weighter.EventWeight(pico, op_all, ctr, updn[i]);
          if (isFastsim) {
            pico.out_sys_fs_bchig()[i]   = btag_weighter.EventWeight(pico, op_all, ctr, ctr, updn[i], ctr);
            pico.out_sys_fs_udsghig()[i] = btag_weighter.EventWeight(pico, op_all, ctr, ctr, ctr, updn[i]);
          }
        }
        lep_weighter.FullSim(pico, w_lep, sys_lep);
        pico.out_w_lep() = w_lep;
        pico.out_sys_lep() = sys_lep;
        if (isFastsim) { 
          lep_weighter.FastSim(pico, w_fs_lep, sys_fs_lep);
          pico.out_w_fs_lep() = w_fs_lep;
          pico.out_sys_fs_lep() = sys_fs_lep;
        }
        photon_weighter.FullSim(pico, w_photon, sys_photon);
        pico.out_w_photon() = w_photon;
        pico.out_sys_photon() = sys_photon;
        if (isZgamma)  { 
          if (year==2016) {
            if(nano.event() % 3516 <= 1887) lep_weighter.FullSim(pico, w_lep, sys_lep);
            else lep_weighter16gh.FullSim(pico, w_lep, sys_lep);
          } else {
            lep_weighter.FullSim(pico, w_lep, sys_lep);
          }
          pico.out_w_lep() = w_lep;
          pico.out_sys_lep() = sys_lep;
        }
        pico.out_w_pu() = 1.; // To be implemented
        pico.out_sys_pu().resize(2, 0.); // Need to be implemented
        // N.B. out_w_prefire should not be renormalized because it models an inefficiency, 
        // i.e. we *should* get less events!
        float w_prefire=1.;
        std::vector<float> sys_prefire(2, 1.);
        prefire_weighter.EventWeight(nano, w_prefire, sys_prefire, isFastsim);
        pico.out_w_prefire() = w_prefire;
        pico.out_sys_prefire() = sys_prefire;
      } // Pre-UL
    } // MC
    if (!isZgamma) pico.out_w_photon() = 1.0;

    // to be calculated in Step 2: merge_corrections
    if (!isData)
      pico.out_w_lumi() = nano.Generator_weight()>0 ? 1:-1;
    else
      pico.out_w_lumi() = 1.;

    //copy LHE scale variation weights and PS weights
    if (!isData) {
      pico.out_sys_murf() = nano.LHEScaleWeight();
      pico.out_sys_ps() = nano.PSWeight();
    }

    isr_tools.WriteISRWeights(pico);


    // do not include w_prefire, or anything that should not be renormalized! Will be set again in Step 3
    if (isZgamma) {
      pico.out_weight() = pico.out_w_lumi() *
                          pico.out_w_lep() * pico.out_w_bhig() * pico.out_w_photon()  *
                          pico.out_w_isr() * pico.out_w_pu() * pico.out_w_trig();
    } else {
      pico.out_weight() = pico.out_w_lumi() *
                          pico.out_w_lep() * pico.out_w_fs_lep() * pico.out_w_bhig() *
                          pico.out_w_isr() * pico.out_w_pu();
    }

    // ----------------------------------------------------------------------------------------------
    //              *** Add up weights to save for renormalization step ***
    // ----------------------------------------------------------------------------------------------
    if (debug) cout<<"INFO:: Writing sum of weights"<<endl;
    if (!isData) {
      wgt_sums.out_weight() += pico.out_weight();
      // taking care of samples with negative weights
      wgt_sums.out_neff() += nano.Generator_weight()>0 ? 1:-1;

      // leptons, keeping track of 0l and 1l totals separately to determine the SF for 0l events
      if(pico.out_nlep()==0){
        wgt_sums.out_nent_zlep() += 1.;
        wgt_sums.out_tot_weight_l0() += pico.out_weight()*(nano.Generator_weight()>0 ? 1:-1); // multiplying by GenWeight to remove the sign...
      }else{
        wgt_sums.out_tot_weight_l1() += pico.out_weight()*(nano.Generator_weight()>0 ? 1:-1);
        wgt_sums.out_w_lep() += w_lep;
        if(isFastsim) wgt_sums.out_w_fs_lep() += w_fs_lep;
        for(size_t i = 0; i<pico.out_sys_lep().size(); ++i){
          wgt_sums.out_sys_lep()[i] += sys_lep[i];
          wgt_sums.out_sys_fs_lep()[i] += sys_fs_lep[i];
        }
      }
      if (pico.out_nel()>0) {
        wgt_sums.out_neff_el() += nano.Generator_weight()>0 ? 1:-1;
        if (pico.out_trig_single_el() || pico.out_trig_double_el()) {
          wgt_sums.out_neff_pass_eltrigs() += nano.Generator_weight()>0 ? 1:-1;
        }
      }
      wgt_sums.out_w_el()      += pico.out_w_el();
      wgt_sums.out_w_mu()      += pico.out_w_mu();
      wgt_sums.out_w_photon()  += pico.out_w_photon();
      wgt_sums.out_w_btag()    += pico.out_w_btag();
      wgt_sums.out_w_btag_df() += pico.out_w_btag_df();
      wgt_sums.out_w_bhig()    += pico.out_w_bhig();
      wgt_sums.out_w_bhig_df() += pico.out_w_bhig_df();
      wgt_sums.out_w_isr()     += pico.out_w_isr();
      wgt_sums.out_w_pu()      += pico.out_w_pu();
      wgt_sums.out_w_trig()    += pico.out_w_trig();

      for(size_t i = 0; i<2; ++i){ 
        wgt_sums.out_sys_el()[i]         += pico.out_sys_el()[i];
        wgt_sums.out_sys_mu()[i]         += pico.out_sys_mu()[i];
        wgt_sums.out_sys_photon()[i]     += pico.out_sys_photon()[i];
        wgt_sums.out_sys_trig()[i]       += pico.out_sys_trig()[i];
        wgt_sums.out_sys_bchig()[i]      += pico.out_sys_bchig()[i];
        wgt_sums.out_sys_udsghig()[i]    += pico.out_sys_udsghig()[i];
        wgt_sums.out_sys_fs_bchig()[i]   += pico.out_sys_fs_bchig()[i];
        wgt_sums.out_sys_fs_udsghig()[i] += pico.out_sys_fs_udsghig()[i];
        wgt_sums.out_sys_isr()[i]        += pico.out_sys_isr()[i];
        wgt_sums.out_sys_pu()[i]         += pico.out_sys_pu()[i];
      }
      for(size_t i = 0; i<pico.out_sys_murf().size(); ++i){ 
        wgt_sums.out_sys_murf()[i] += pico.out_sys_murf()[i];
      }
      for(size_t i = 0; i<pico.out_sys_ps().size(); ++i){ 
        wgt_sums.out_sys_ps()[i] += pico.out_sys_ps()[i];
      }
    }
    
    if (debug) cout<<"INFO:: Filling tree"<<endl;
    pico.Fill();
  } // loop over events

  wgt_sums.Fill();
  wgt_sums.Write();
  pico.Write();

  cout<<endl;
  time(&endtime); 
  cout<<"Time passed: "<<hoursMinSec(difftime(endtime, begtime))<<endl<<endl;  
}

void Initialize(corrections_tree &wgt_sums){
  wgt_sums.out_neff()              = 0;
  wgt_sums.out_nent_zlep()         = 0;
  wgt_sums.out_neff_el()           = 0;
  wgt_sums.out_neff_pass_eltrigs() = 0;
  wgt_sums.out_tot_weight_l0()     = 0.;
  wgt_sums.out_tot_weight_l1()     = 0.;

  wgt_sums.out_weight()      = 0.;
  wgt_sums.out_w_lumi()      = 0.;
  wgt_sums.out_w_el()        = 0.;
  wgt_sums.out_w_mu()        = 0.;
  wgt_sums.out_w_lep()       = 0.;
  wgt_sums.out_w_fs_lep()    = 0.;
  wgt_sums.out_w_photon()    = 0.;
  wgt_sums.out_w_btag()      = 0.;
  wgt_sums.out_w_btag_df()   = 0.;
  wgt_sums.out_w_bhig()      = 0.;
  wgt_sums.out_w_bhig_df()   = 0.;
  wgt_sums.out_w_isr()       = 0.;
  wgt_sums.out_w_pu()        = 0.;
  wgt_sums.out_w_trig()      = 0.;
  wgt_sums.out_w_zvtx_pass() = 0.;
  wgt_sums.out_w_zvtx_fail() = 0.;
  // w_prefire should not be normalized

  wgt_sums.out_sys_el().resize(2,0);
  wgt_sums.out_sys_mu().resize(2,0);
  wgt_sums.out_sys_lep().resize(2,0);
  wgt_sums.out_sys_fs_lep().resize(2,0);
  wgt_sums.out_sys_photon().resize(2,0);
  wgt_sums.out_sys_bchig().resize(2,0);
  wgt_sums.out_sys_udsghig().resize(2,0);
  wgt_sums.out_sys_fs_bchig().resize(2,0);
  wgt_sums.out_sys_fs_udsghig().resize(2,0);
  wgt_sums.out_sys_isr().resize(2,0);
  wgt_sums.out_sys_pu().resize(2,0);
  wgt_sums.out_sys_trig().resize(2,0);
  wgt_sums.out_sys_murf().resize(9,0);
  wgt_sums.out_sys_ps().resize(4,0);
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"in_file", required_argument, 0,'f'}, 
      {"in_dir",  required_argument, 0,'i'},
      {"out_dir", required_argument, 0,'o'},
      {"nent",    required_argument, 0, 0},
      {"debug",    no_argument, 0, 'd'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "f:i:o:d", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'f':
      in_file = optarg;
      break;
    case 'i':
      in_dir = optarg;
      break;
    case 'd':
      debug = true;
      break;
    case 'o':
      out_dir = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "nent"){
        nent_test = atoi(optarg);
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
        exit(1);
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
