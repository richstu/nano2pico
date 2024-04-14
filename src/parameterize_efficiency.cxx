#include <iostream>
#include <iomanip>
#include <getopt.h>
#include "TChain.h"
#include "TBranch.h"
#include "TH3D.h"
#include "TFile.h"

namespace { 
  TString infile = "";
  int year = 2016;
}

using namespace std;

void GetOptions(int argc, char *argv[]);
void UpdateProgressBar(float progress);

int main(int argc, char *argv[]) {
    GetOptions(argc, argv);
    TChain c("Events");
    infile = infile +"TTJets_*4.root"; // Don't need full TTJet sample, only use a fraction
    c.Add(infile);
    // Set appropriate working points
    vector<double> DeepCSV_WPs;
    vector<double> DeepFlavor_WPs;
    if(year == 2016) {
      DeepCSV_WPs =    {0.2217, 0.6321, 0.8953};
      DeepFlavor_WPs = {0.0614, 0.3093, 0.7221}; }
    else if(year == 2017) {
      DeepCSV_WPs =    {0.1522, 0.4941, 0.8001};
      DeepFlavor_WPs = {0.0521, 0.3033, 0.7489}; }
    else { 
      DeepCSV_WPs =    {0.1241, 0.4184, 0.7527};
      DeepFlavor_WPs = {0.0494, 0.2770, 0.7264}; }
    double eta_cuts[3] = {0., 1.2, 2.4};
    double pt_cuts[9] = {30., 50., 70., 100., 140., 200., 300., 670., 1.e4};
    double flavor_cuts[4] = {-0.5, 3.5, 4.5, 5.5};
    vector<string> WPs = {"loose", "medium", "tight"};
    // Define TH3Ds for each algorthm and working point
    vector<TH3D> numerators = {}, numerators_df = {};
    for(size_t iwp = 0; iwp < WPs.size(); iwp++) {
      TString name("btagEfficiency_DeepCSV_"+WPs[iwp]);
      TString name_df("btagEfficiency_DeepFlavor_"+WPs[iwp]);
      numerators.push_back(TH3D(name,name, 2, eta_cuts, 8, pt_cuts, 3, flavor_cuts));
      numerators_df.push_back(TH3D(name_df,name_df, 2, eta_cuts, 8, pt_cuts, 3, flavor_cuts));
      }
    TH3D denominator = TH3D("btagEfficiency", "btagEfficiency", 2, eta_cuts, 8, pt_cuts, 3, flavor_cuts);
    // Assign branch variables and turn off unused branches
    int nJet;
    array<float,100> Jet_pt, Jet_eta, Jet_btagDeepB, Jet_btagDeepFlavB;
    array<int,100> Jet_hadronFlavour, Jet_nElectrons, Jet_nMuons;
    TBranch *b_Jet_pt = nullptr, *b_Jet_eta = nullptr, *b_Jet_btagDeepB = nullptr, *b_Jet_btagDeepFlavB = nullptr;
    TBranch *b_Jet_hadronFlavour = nullptr, *b_Jet_nElectrons = nullptr, *b_Jet_nMuons = nullptr;
    c.SetBranchStatus("*",0);
    c.SetBranchStatus("nJet",             1); c.SetBranchAddress("nJet",             &nJet);
    c.SetBranchStatus("Jet_pt",           1); c.SetBranchAddress("Jet_pt",           &Jet_pt,            &b_Jet_pt);
    c.SetBranchStatus("Jet_eta",          1); c.SetBranchAddress("Jet_eta",          &Jet_eta,           &b_Jet_eta);
    c.SetBranchStatus("Jet_btagDeepB",    1); c.SetBranchAddress("Jet_btagDeepB",    &Jet_btagDeepB,     &b_Jet_btagDeepB);
    c.SetBranchStatus("Jet_btagDeepFlavB",1); c.SetBranchAddress("Jet_btagDeepFlavB",&Jet_btagDeepFlavB, &b_Jet_btagDeepFlavB);
    c.SetBranchStatus("Jet_hadronFlavour",1); c.SetBranchAddress("Jet_hadronFlavour",&Jet_hadronFlavour, &b_Jet_hadronFlavour);
    c.SetBranchStatus("Jet_nElectrons",   1); c.SetBranchAddress("Jet_nElectrons",   &Jet_nElectrons,    &b_Jet_nElectrons);
    c.SetBranchStatus("Jet_nMuons",       1); c.SetBranchAddress("Jet_nMuons",       &Jet_nMuons,        &b_Jet_nMuons);
    int num_entries = c.GetEntries();
    float flavor, pt, eta, deepcsv, deepflav;
    for(int entry = 0; entry < c.GetEntries(); entry++) { // Loop over events
        c.GetEntry(entry);
        if(entry % (num_entries/100) == 0) 
            UpdateProgressBar(float(entry)/num_entries);
        for(int ijet = 0; ijet < nJet; ijet++) {
            if (Jet_nElectrons[ijet] != 0 || Jet_nMuons[ijet] != 0 || Jet_pt[ijet] < 30f || abs(Jet_eta[ijet]) > 2.4f) 
              continue;
            flavor = abs(Jet_hadronFlavour[ijet]);
            pt = Jet_pt[ijet];
            eta = abs(Jet_eta[ijet]);
            deepcsv  = Jet_btagDeepB[ijet];
            deepflav = Jet_btagDeepFlavB[ijet];
            denominator.Fill(eta, pt, flavor);
            for(int iwp = 0; iwp < 3; iwp++) {
                if(deepcsv > DeepCSV_WPs[iwp])
                  numerators[iwp].Fill(eta, pt, flavor);
                if(deepflav > DeepFlavor_WPs[iwp])
                  numerators_df[iwp].Fill(eta, pt, flavor);
            }
        }
    }
    cout << endl << "Done!" << endl;
    // Divide for efficiency 
    for(int iwp = 0; iwp < 3; iwp++) {
        numerators[iwp].Sumw2(); // Needed for errors to scale correctly
        numerators_df[iwp].Sumw2();
        numerators[iwp] = numerators[iwp]/denominator;
        numerators_df[iwp] = numerators_df[iwp]/denominator;
    }
    // Write efficiencies to outfile
    TString outname = "data/btagEfficiency_DeepCSV_"+to_string(year)+".root";            
    TFile* out_file = new TFile(outname, "recreate");
    for(int iwp = 0; iwp < 3; iwp++) 
        numerators[iwp].Write();
    out_file->Close();

    outname = "data/btagEfficiency_DeepFlavor_"+to_string(year)+".root";            
    TFile* out_file_df = new TFile(outname, "recreate");
    for(int iwp = 0; iwp < 3; iwp++) 
        numerators_df[iwp].Write();
    out_file_df->Close();
    }

void UpdateProgressBar(float progress) { // Because why not?
  cout << " Running: [";
  progress = int(progress*100);
  for(int p=0; p<100;++p) 
    if(p<progress) cout << "=";
    else if(p==progress) cout << ">";
    else cout << " ";
  cout << "] " << int(progress) << " %\r";
  cout.flush();
  }

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"year", required_argument, 0, 'y'},
      {"infile", required_argument, 0, 'i'}       
    };
    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "y:i:", long_options, &option_index);
    if(opt == -1) 
      break;
    string optname;
    switch(opt){
    case 'y':
      year = atoi(optarg);
      break;
    case 'i':
      infile = optarg;
      break;
    case 0:
      printf("Need to specify MC path and year!\n");
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}

