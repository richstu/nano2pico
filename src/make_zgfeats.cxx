#include <ctime>

#include <iostream>
#include <iomanip>

#include <getopt.h>

#include "TError.h"
#include "Math/Vector4D.h"

#include "pico_tree.hpp"
#include "zgfeats_tree.hpp"
#include "utilities.hpp"

using namespace std;

namespace {
  string in_file = "merged_pico_higloose_met150_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_higmc_higloose_nfiles_78.root";
  string in_dir = "/cms29r0/pico/NanoAODv5/zgamma_channelIslandsv2/2016/mc/skim_zcand/";
  string out_dir = "features";
  int nent_test = -1;
  bool debug = false;
}

void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches       
  GetOptions(argc, argv);

  if(in_file=="" || in_dir=="" || out_dir == "") {
    cout<<"ERROR: Input file, sum-of-weights and/or output directory not specified. Exit."<<endl;
    exit(0);
  }

  string in_path = in_dir+"/"+in_file;
  string out_path = out_dir+"/features_"+in_file;

  time_t begtime, endtime;
  time(&begtime);

 // Initialize trees
  pico_tree pico(in_path);
  size_t nentries(nent_test>0 ? nent_test : pico.GetEntries());
  cout << "Pico input file: " << in_path << endl;
  cout << "Input number of events: " << nentries << endl;

  zgfeats_tree features("", out_path);
  cout << "Writing output to: " << out_path << endl;

  for(size_t entry(0); entry<nentries; ++entry){
    if (debug) cout << "GetEntry: " << entry << endl;
    pico.GetEntry(entry);
    if (entry%1000==0 || entry == nentries-1) {
      cout<<"Processing event: "<<entry<<endl;
    }

    if (debug) cout<<"INFO:: Filling global features."<<endl;
    features.out_event()   = pico.event();
    features.out_nlep()    = pico.nlep();
    features.out_njet()    = pico.njet();
    features.out_nphoton() = pico.nphoton();

    //--------------------------------------------------------------
    // Save two highest pT jets (and combination)
    //--------------------------------------------------------------
    if (debug) cout<<"INFO:: Filling Ak4 jet features."<<endl;
    int jetn = 0;
    for (size_t ijet(0); ijet<pico.jet_pt().size(); ijet++) {
      if(pico.jet_isgood()[ijet]) {
        if(jetn == 0) {
          features.out_jet1_pt()  = pico.jet_pt()[ijet];
          features.out_jet1_eta() = pico.jet_eta()[ijet];
          features.out_jet1_phi() = pico.jet_phi()[ijet];
          features.out_jet1_m()   = pico.jet_m()[ijet];
          jetn++;
        }
        else if(jetn == 1) {
          features.out_jet2_pt()  = pico.jet_pt()[ijet];
          features.out_jet2_eta() = pico.jet_eta()[ijet];
          features.out_jet2_phi() = pico.jet_phi()[ijet];
          features.out_jet2_m()   = pico.jet_m()[ijet];
          features.out_dijet_pt()   = pico.dijet_pt();
          features.out_dijet_eta()  = pico.dijet_eta();
          features.out_dijet_phi()  = pico.dijet_phi();
          features.out_dijet_m()    = pico.dijet_m();
          features.out_dijet_dr()   = pico.dijet_dr();
          features.out_dijet_dphi() = pico.dijet_dphi();
          features.out_dijet_deta() = pico.dijet_deta();
          continue;
        }
      }
    }
    //--------------------------------------------------------------
    // Save highest pT photon candidate
    //--------------------------------------------------------------
    if (debug) cout<<"INFO:: Filling photon features."<<endl;
    for(size_t iph(0); iph < pico.photon_pt().size(); iph++)
      if(pico.photon_sig()[iph]) {
        features.out_photon_pt()    = pico.photon_pt()[iph];
        features.out_photon_eta()   = pico.photon_eta()[iph];
        features.out_photon_phi()   = pico.photon_phi()[iph];
        features.out_photon_drmin() = pico.photon_drmin()[iph];
        features.out_photon_drmax() = pico.photon_drmax()[iph];
        continue;
      }
    //--------------------------------------------------------------
    // Find best Z candidate and save attributes of Zcand and leptons
    //--------------------------------------------------------------
    if (debug) cout<<"INFO:: Filling Z candidate features."<<endl;
    double mindm(100000);
    int biz(-1);
    for(int iz(0); iz < pico.nll(); iz++)
      if(abs(pico.ll_m()[iz] - 91.1876) < mindm)
        biz = iz;
    if(biz != -1) {
      features.out_ll_pt()    = pico.ll_pt()[biz];
      features.out_ll_eta()   = pico.ll_eta()[biz];
      features.out_ll_phi()   = pico.ll_phi()[biz];
      features.out_ll_m()     = pico.ll_m()[biz];
      features.out_ll_dr()    = pico.ll_dr()[biz];
      features.out_ll_dphi()  = pico.ll_dphi()[biz];
      features.out_ll_deta()  = pico.ll_deta()[biz];
      features.out_ll_lepid() = pico.ll_lepid()[biz];
      int il1(pico.ll_i1()[biz]), il2(pico.ll_i2()[biz]);
      if(pico.ll_lepid()[biz] == 11) {
        features.out_lep1_pt()  = pico.el_pt()[il1];
        features.out_lep1_eta() = pico.el_eta()[il1];
        features.out_lep1_phi() = pico.el_phi()[il1];
        features.out_lep2_pt()  = pico.el_pt()[il2];
        features.out_lep2_eta() = pico.el_eta()[il2];
        features.out_lep2_phi() = pico.el_phi()[il2];
      }
      else{
        features.out_lep1_pt()  = pico.mu_pt()[il1];
        features.out_lep1_eta() = pico.mu_eta()[il1];
        features.out_lep1_phi() = pico.mu_phi()[il1];
        features.out_lep2_pt()  = pico.mu_pt()[il2];
        features.out_lep2_eta() = pico.mu_eta()[il2];
        features.out_lep2_phi() = pico.mu_phi()[il2];
      }
    }
    //--------------------------------------------------------------
    // Find llg with best Z candidate and save attributes
    //--------------------------------------------------------------
    if (debug) cout<<"INFO:: Filling Higgs candidate features."<<endl;
    for(int izg(0); izg < pico.nllphoton(); izg++)
      if(pico.llphoton_ill()[izg] == biz) {
        features.out_llphoton_pt()     = pico.llphoton_pt()[izg];
        features.out_llphoton_eta()    = pico.llphoton_eta()[izg];
        features.out_llphoton_phi()    = pico.llphoton_phi()[izg];
        features.out_llphoton_m()      = pico.llphoton_m()[izg];
        features.out_llphoton_dr()     = pico.llphoton_dr()[izg];
        features.out_llphoton_dphi()   = pico.llphoton_dphi()[izg];
        features.out_llphoton_deta()   = pico.llphoton_deta()[izg];
        //features.out_llphoton_costhj() = pico.llphoton_costhj()[izg];
        features.out_llphoton_costhj() = pico.llphoton_costheta()[izg];
        continue;
      }
      
    if (debug) cout<<"INFO:: Filling tree"<<endl;
    features.Fill();
  } // loop over events

  features.Write();

  cout<<endl;
  time(&endtime); 
  cout<<"Time passed: "<<hoursMinSec(difftime(endtime, begtime))<<endl<<endl;  
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
