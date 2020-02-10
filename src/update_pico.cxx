#include <ctime>
#include <iostream>

#include <getopt.h>

#include "TTree.h"
#include "TFile.h"
#include "TObject.h"
#include "TBranch.h"

#include "utilities.cpp"

using namespace std;

namespace {
  TString dnnout_file = "";
  TString pico_file = "";
}

void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  TFile fdnn(dnnout_file,"read"); 
  TTree *tree_dnn = static_cast<TTree*>(fdnn.Get("tree")); 

  TFile fpico(pico_file,"update"); 
  TTree *tree_pico = static_cast<TTree*>(fpico.Get("tree")); 
  cout<<"Updating tree in:"<<pico_file<<endl;
  
  float hig_am_dnn_pico; 
  TBranch *bhig_am_dnn_pico = tree_pico->Branch("hig_am_dnn",&hig_am_dnn_pico,"hig_am_dnn/F");

  double hig_am_dnn;
  tree_dnn->SetBranchAddress("hig_am_dnn",&hig_am_dnn); 
  size_t nentries = tree_pico->GetEntries(); 
  for (size_t i=0;i<nentries;i++) { 
    tree_pico->GetEntry(i); 
    tree_dnn->GetEntry(i); 
    hig_am_dnn_pico = hig_am_dnn; 
    bhig_am_dnn_pico->Fill(); 
  } 
  tree_pico->Print(); 
  tree_pico->Write("",TObject::kWriteDelete);  

  cout<<endl;
  time(&endtime); 
  cout<<"Time passed: "<<hoursMinSec(difftime(endtime, begtime))<<endl<<endl; 

}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"dnnout_file", required_argument, 0, 'd'},
      {"pico_file", required_argument, 0, 'p'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "d:p:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'd':
      dnnout_file = optarg;
      break;
    case 'p':
      pico_file = optarg;
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      exit(1);
      break;
    }
  }
}