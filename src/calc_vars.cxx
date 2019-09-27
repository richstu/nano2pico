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
  bool isFastSim = false;
  int year = 2016;
  const size_t nent_test = 5;
}

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
  // cout << "Running on "<< (isFastSim ? "FastSim" : "FullSim") << endl;
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

  BTagWeighter btw(isFastSim, year);
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

    // N.B. Order in which producers are called matters! E.g. jets are not counted if overlapping 
    // with signal lepton, thus jets must be processed only after leptons have been selected.
    mc_producer.WriteGenParticles(nano, pico);

    el_producer.WriteElectrons(nano, pico);
    mu_producer.WriteMuons(nano, pico);
    dilep_producer.WriteDileptons();
    tk_producer.WriteIsoTracks(nano, pico);
    ph_producer.WritePhotons(nano, pico);

    jet_producer.WriteJets(nano, pico);
    fjet_producer.WriteFatJets(nano, pico);

    hig_producer.WriteHigVars();

    pico.Fill();
  } // loop over events

  corr.Fill();
  corr.Write();
  pico.Write();

  cout<<endl;
  time(&endtime); 
  cout<<"Time passed: "<<hoursMinSec(difftime(endtime, begtime))<<endl<<endl;  
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"nano_file", required_argument, 0, 'n'},
      {"wgt_sums_file", required_argument, 0, 'w'}, 
      {"pico_file", required_argument, 0, 'p'},
      {"year", required_argument, 0, 'y'},  
      {"fastsim", no_argument, 0, 0},       
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
      if(optname == "fastsim"){
        isFastSim = true;
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
