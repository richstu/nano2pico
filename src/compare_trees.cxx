#include <ctime>

#include <iostream>

#include <getopt.h>

#include "TError.h"
#include "TH1D.h"
#include "TFile.h"

#include "baby_tree.hpp"
#include "pico_tree.hpp"
#include "utilities.hpp"

using namespace std;

namespace {
  string baby_file = "";
  string pico_file = "";
  int nent_test = -1;
}

void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  // gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches       
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  pico_tree pico(pico_file);
  baby_tree baby(baby_file);

  map<string, TH1D> histos;
  histos["nel"] = TH1D("nel","nel; dx(pico, baby)",6,-3,3);
  histos["nmu"] = TH1D("nmu","nmu",6,-3,3);
  histos["mu_reliso"] = TH1D("mu_reliso","mu_reliso",40,-20,20);
  histos["mu_miniso"] = TH1D("mu_miniso","mu_miniso",40,-5,5);
  // histos[""] = TH1D("njet","njet; Jets dN(pico, baby)",6,-3,3));
  // histos[""] = TH1D("jet_pt","jet_pt; Jets dPt(pico, baby) [GeV]",40,-2,2));
  // histos[""] = TH1D("jet_eta","jet_eta; Jets deta(pico, baby)",40,-2,2));

  size_t nentries(nent_test>0 ? nent_test : pico.GetEntries());
  for(size_t entry(0); entry<nentries; ++entry){
    pico.GetEntry(entry);
    baby.GetEntry(entry);

    histos["nel"].Fill(pico.nel()-baby.nels());

    histos["nmu"].Fill(pico.nmu()-baby.nmus());
    size_t mu_pt_size = pico.mu_pt().size() < baby.mus_pt().size() ? pico.mu_pt().size() : baby.mus_pt().size(); 
    for (size_t imu(0); imu<mu_pt_size; imu++) {
      histos["mu_reliso"].Fill(pico.mu_reliso()[imu]*pico.mu_pt()[imu]-baby.mus_miniso()[imu]*baby.mus_pt()[imu]);
      histos["mu_miniso"].Fill(pico.mu_miniso()[imu]-baby.mus_miniso()[imu]);
    }

    // size_t jet_pt_size = pico.jet_pt().size() < baby.jets_pt().size() ? pico.jet_pt().size() : baby.jets_pt().size(); 
    // for (size_t ijet(0); ijet<jet_pt_size; ijet++) {
    //   histos[4].Fill(pico.njet()-baby.njets());
    //   histos[5].Fill(pico.jet_pt()[ijet]-baby.jets_pt()[ijet]);
    //   histos[6].Fill(pico.jet_eta()[ijet]-baby.jets_eta()[ijet]);
    // }
  }

  TFile fhist("histos.root","recreate");
  for (auto &h: histos) h.second.Write();
  fhist.Close();

  cout<<endl;
  time(&endtime); 
  cout<<"Time passed: "<<hoursMinSec(difftime(endtime, begtime))<<endl<<endl; 

}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"nent", required_argument, 0, 0},
      {"baby", required_argument, 0, 'b'},
      {"pico", required_argument, 0, 'p'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "b:p:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'b':
      baby_file = optarg;
      break;
    case 'p':
      pico_file = optarg;
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