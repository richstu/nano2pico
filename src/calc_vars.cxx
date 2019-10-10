#include <ctime>

#include <iostream>

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
#include "ph_producer.hpp"
#include "jet_producer.hpp"
#include "fjet_producer.hpp"
#include "hig_producer.hpp"

#include "btag_weighter.hpp"
#include "lepton_weighter.hpp"

using namespace std;

namespace {
  string nano_file = "";
  string wgt_sums_file = "";
  string pico_file = "";
  bool isData = false;
  bool isFastsim = false;
  int year = 2016;
  int nent_test = -1;
}

void WriteDataQualityFilters(nano_tree& nano, pico_tree& pico);
void CopyTriggerDecisions(nano_tree& nano, pico_tree& pico);

void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  // gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches       
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  if(nano_file=="" || wgt_sums_file == "" || pico_file == "") {
    cout<<"ERROR: Input file, sum-of-weights and/or output directory not specified. Exit."<<endl;
    exit(0);
  }

  // Initialize trees
  nano_tree nano(nano_file);
  size_t nentries(nent_test>0 ? nent_test : nano.GetEntries());
  cout << "Nano file: " << nano_file << endl;
  cout << "Input number of events: " << nentries << endl;
  // cout << "Running on "<< (isFastsim ? "FastSim" : "FullSim") << endl;
  // cout << "Calculating weights based on " << year << " scale factors." << endl;

  pico_tree pico("", pico_file);
  cout << "Writing output to: " << pico_file << endl;

  corrections_tree corr("", wgt_sums_file);
  cout << "Writing sum-of-weights to: " << wgt_sums_file << endl;

  corr.out_nent_zlep() = 0.;
  corr.out_tot_weight_l0() = 0.;
  corr.out_tot_weight_l1() = 0.;
  corr.out_nent() = nentries;

  //Initialize object producers
  GenParticleProducer mc_producer(year);
  ElectronProducer el_producer(year);
  MuonProducer mu_producer(year);
  DileptonProducer dilep_producer(year);
  IsoTrackProducer tk_producer(year);
  PhotonProducer ph_producer(year);
  JetProducer jet_producer(year);
  FatJetProducer fjet_producer(year);
  HigVarProducer hig_producer(year);

  BTagWeighter btw(isFastsim, year);
  LeptonWeighter lw(year);
  const string ctr = "central";
  const string vup = "up";
  const string vdown = "down";
  // const auto op_loose = BTagEntry::OP_LOOSE;
  // const auto op_med = BTagEntry::OP_MEDIUM;
  // const auto op_tight = BTagEntry::OP_TIGHT;
  const vector<BTagEntry::OperatingPoint> op_all = {BTagEntry::OP_LOOSE, BTagEntry::OP_MEDIUM, BTagEntry::OP_TIGHT};

  for(size_t entry(0); entry<nentries; ++entry){
    nano.GetEntry(entry);
    if (entry%100==0 || entry == nentries-1) {
      cout<<"Processing event: "<<entry<<endl;
    }
    // if (nano.event()!=534935544) continue;

    pico.out_event() = nano.event();
    // N.B. Order in which producers are called matters! E.g. jets are not counted if overlapping 
    // with signal lepton, thus jets must be processed only after leptons have been selected.
    mc_producer.WriteGenParticles(nano, pico);

    vector<int> jet_islep_nano_idx = vector<int>();
    vector<int> sig_el_nano_idx = el_producer.WriteElectrons(nano, pico, jet_islep_nano_idx);
    vector<int> sig_mu_nano_idx = mu_producer.WriteMuons(nano, pico, jet_islep_nano_idx);
    dilep_producer.WriteDielectrons(nano, pico, sig_el_nano_idx);
    dilep_producer.WriteDimuons(nano, pico, sig_mu_nano_idx);

    tk_producer.WriteIsoTracks(nano, pico, sig_el_nano_idx, sig_mu_nano_idx);
    ph_producer.WritePhotons(nano, pico);

    jet_producer.WriteJets(nano, pico, jet_islep_nano_idx);
    fjet_producer.WriteFatJets(nano, pico);

    hig_producer.WriteHigVars();

    // N.B. Jets: pico.out_pass_jets() and pico.out_pass_fsjets() filled in jet_producer
    // so this should be called after writeJets to allow computing overall pass variable
    WriteDataQualityFilters(nano, pico);

    CopyTriggerDecisions(nano, pico);

    pico.Fill();
  } // loop over events

  corr.Fill();
  corr.Write();
  pico.Write();

  cout<<endl;
  time(&endtime); 
  cout<<"Time passed: "<<hoursMinSec(difftime(endtime, begtime))<<endl<<endl;  
}

void WriteDataQualityFilters(nano_tree& nano, pico_tree& pico){

  pico.out_pass_hbhe() = nano.Flag_HBHENoiseFilter();
  pico.out_pass_hbheiso() = nano.Flag_HBHENoiseIsoFilter();
  pico.out_pass_goodv() = nano.Flag_goodVertices();
  pico.out_pass_cschalo_tight() = nano.Flag_globalSuperTightHalo2016Filter();
  pico.out_pass_eebadsc() = nano.Flag_eeBadScFilter();
  pico.out_pass_ecaldeadcell() = nano.Flag_EcalDeadCellTriggerPrimitiveFilter();
  if (year==2016) {
    pico.out_pass_badcalib() = true;
    pico.out_pass_badpfmu() = nano.Flag_BadPFMuonSummer16Filter();
    pico.out_pass_badchhad() = nano.Flag_BadChargedCandidateSummer16Filter();
  } else {
    pico.out_pass_badcalib() = nano.Flag_ecalBadCalibFilterV2();
    pico.out_pass_badpfmu() = nano.Flag_BadPFMuonFilter();
    pico.out_pass_badchhad() = nano.Flag_BadChargedCandidateFilter();
  }
  pico.out_pass_mubadtrk() = nano.Flag_muonBadTrackFilter();

  // Combined pass variable
  pico.out_pass() = pico.out_pass_ra2_badmu() && pico.out_pass_badpfmu() && 
                    pico.out_met()/pico.out_met_calo()<5 &&
                    pico.out_pass_goodv() &&
                    pico.out_pass_hbhe() && pico.out_pass_hbheiso() && 
                    pico.out_pass_ecaldeadcell() && pico.out_pass_badcalib();

  if (isFastsim) {
    pico.out_pass() = pico.out_pass() && pico.out_pass_fsjets();
  } else {
    pico.out_pass() = pico.out_pass() && pico.out_pass_jets() && pico.out_pass_cschalo_tight();
    if (isData) 
      pico.out_pass() = pico.out_pass() && pico.out_pass_eebadsc();
  }
  
  return;
}

void CopyTriggerDecisions(nano_tree& nano, pico_tree& pico){
  pico.out_HLT_IsoMu24() = nano.HLT_IsoMu24();
  pico.out_HLT_IsoMu27() = nano.HLT_IsoMu27();
  pico.out_HLT_Mu50() = nano.HLT_Mu50();
  pico.out_HLT_Ele27_WPTight_Gsf() = nano.HLT_Ele27_WPTight_Gsf();
  pico.out_HLT_Ele35_WPTight_Gsf() = nano.HLT_Ele35_WPTight_Gsf();
  pico.out_HLT_Ele115_CaloIdVT_GsfTrkIdT() = nano.HLT_Ele115_CaloIdVT_GsfTrkIdT();

  // MET triggers
  pico.out_HLT_PFMET110_PFMHT110_IDTight() = nano.HLT_PFMET110_PFMHT110_IDTight();
  pico.out_HLT_PFMET120_PFMHT120_IDTight() = nano.HLT_PFMET120_PFMHT120_IDTight();
  pico.out_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight() = nano.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight();
  pico.out_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight() = nano.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight();
  pico.out_HLT_PFMET120_PFMHT120_IDTight_PFHT60() = nano.HLT_PFMET120_PFMHT120_IDTight_PFHT60();
  pico.out_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60() = nano.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60();

  // Jet trigger
  pico.out_HLT_PFJet500() = nano.HLT_PFJet500();

  // ZGamma triggers
  pico.out_HLT_Mu17_Photon30_IsoCaloId() = nano.HLT_Mu17_Photon30_IsoCaloId();
  pico.out_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() = nano.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ();
  pico.out_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL() = nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL();
  pico.out_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ() = nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ();
  return;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"nent", required_argument, 0, 0},
      {"nano_file", required_argument, 0, 'n'},
      {"wgt_sums_file", required_argument, 0, 'w'}, 
      {"pico_file", required_argument, 0, 'p'},
      {"year", required_argument, 0, 'y'},  
      {"isFastsim", no_argument, 0, 0},       
      {"isData", no_argument, 0, 0},       
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "n:w:p:y:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'n':
      nano_file = optarg;
      break;
    case 'w':
      wgt_sums_file = optarg;
      break;
    case 'p':
      pico_file = optarg;
      break;
    case 'y':
      year = atoi(optarg);
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "isFastsim"){
        isFastsim = true;
      } if(optname == "isData"){
        isData = true;
      } else if(optname == "nent"){
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
