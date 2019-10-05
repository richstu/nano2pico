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
  histos["nvjet"] = TH1D("nvjet","nvjet; Jets dN(pico, baby)",8,-4,4);
  histos["njet"] = TH1D("njet","njet; Jets dN(pico, baby)",8,-4,4);
  histos["jet_pt"] = TH1D("jet_pt","jet_pt; Jets dPt(pico, baby) [GeV]",100,-1,1);
  histos["jet_eta"] = TH1D("jet_eta","jet_eta; Jets deta(pico, baby)",40,-2,2);
  histos["p_isele"] = TH1D("","",2,0,2);
  histos["b_isele"] = TH1D("","",2,0,2);
  histos["p_ismu"] = TH1D("","",2,0,2);
  histos["b_ismu"] = TH1D("","",2,0,2);

  size_t nentries(nent_test>0 ? nent_test : pico.GetEntries());
  for(size_t entry(0); entry<nentries; ++entry){
    pico.GetEntry(entry);
    baby.GetEntry(entry);

    histos["nel"].Fill(pico.nel()-baby.nels());
    histos["nmu"].Fill(pico.nmu()-baby.nmus());

    bool isele = (pico.nel()-baby.nels())==0 && pico.nmu()==0 && baby.nmus()==0;
    bool ismu = (pico.nmu()-baby.nmus())==0 && pico.nel()==0 && baby.nels()==0;
    
    histos["njet"].Fill(pico.njet()-baby.njets());
    histos["nvjet"].Fill(pico.jet_pt().size()-baby.jets_pt().size());
    size_t jet_pt_size = pico.jet_pt().size() < baby.jets_pt().size() ? pico.jet_pt().size() : baby.jets_pt().size(); 
    for (size_t ijet(0); ijet<jet_pt_size; ijet++) {
      histos["jet_pt"].Fill(pico.jet_pt()[ijet]-baby.jets_pt()[ijet]);
      histos["jet_eta"].Fill(pico.jet_eta()[ijet]-baby.jets_eta()[ijet]);
      if (isele) {
        histos["p_isele"].Fill(pico.jet_islep()[ijet]);
        histos["b_isele"].Fill(baby.jets_islep()[ijet]);
      } else if (ismu) {
        histos["p_ismu"].Fill(pico.jet_islep()[ijet]);
        histos["b_ismu"].Fill(baby.jets_islep()[ijet]);
      }
    }
  }

  cout<<"Difference in number of signal electrons (total = "<<histos["nel"].GetEntries()<<")"
      <<histos["nel"].GetMean()<<"+-"<<histos["nel"].GetRMS()<<endl;
  cout<<"Mean+-rms of the nvjet disctribution is "<<histos["nvjet"].GetMean()<<"+-"<<histos["nvjet"].GetRMS()<<endl<<endl;
  cout<<"Pico isele jets -> "<<histos["p_isele"].GetMean()<<"+-"<<histos["p_isele"].GetRMS()<<endl;
  cout<<"Baby isele jets -> "<<histos["b_isele"].GetMean()<<"+-"<<histos["b_isele"].GetRMS()<<endl<<endl;
  cout<<"Pico ismu jets -> "<<histos["p_ismu"].GetMean()<<"+-"<<histos["p_ismu"].GetRMS()<<endl;
  cout<<"Baby ismu jets -> "<<histos["b_ismu"].GetMean()<<"+-"<<histos["b_ismu"].GetRMS()<<endl;
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