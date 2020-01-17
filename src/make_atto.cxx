#include <ctime>

#include <iostream>
#include <iomanip>

#include <getopt.h>

#include "TError.h"
#include "Math/Vector4D.h"

#include "nano_tree.hpp"
#include "atto_tree.hpp"
#include "utilities.hpp"

using namespace std;

namespace {
  string in_file = "TChiHH_HToBB_HToBB_3D_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv5__PUMoriond17_run2_nanoAOD_94X2016_102X_mcRun2_asymptotic_v7_file501.root";
  string in_dir = "../data/";
  string out_dir = "out";
  int nent_test = -1;
  bool debug = false;
  // requirements for jets to be counted in njet, mofified for Zgamma below
  float min_jet_pt = 30.0;
  float max_jet_eta =  2.4;
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
  string out_path = out_dir+"/raw_atto_"+in_file;

  int year = Contains(in_file, "RunIISummer16") ? 2016 : (Contains(in_file, "RunIIFall17") ? 2017 : 2018);

  time_t begtime, endtime;
  time(&begtime);

  // B-tag working points
  map<int, vector<float>> btag_wpts{
    {2016, vector<float>({0.2217, 0.6321, 0.8953})},
    {2017, vector<float>({0.1522, 0.4941, 0.8001})},
    {2018, vector<float>({0.1241, 0.4184, 0.7527})}
  };

 // Initialize trees
  nano_tree nano(in_path);
  size_t nentries(nent_test>0 ? nent_test : nano.GetEntries());
  cout << "Nano file: " << in_path << endl;
  cout << "Input number of events: " << nentries << endl;

  atto_tree atto("", out_path);
  cout << "Writing output to: " << out_path << endl;

  for(size_t entry(0); entry<nentries; ++entry){
    if (debug) cout << "GetEntry: " << entry << endl;
    nano.GetEntry(entry);
    if (entry%1000==0 || entry == nentries-1) {
      cout<<"Processing event: "<<entry<<endl;
    }

    atto.out_event() = nano.event();
    atto.out_met()   = nano.MET_pt();

    //--------------------------------------------------------------
    //         Save jets ordered by deepCSV value    
    //--------------------------------------------------------------
    vector<pair<int, float>>  ordered_idx;
    for (int ijet(0); ijet<nano.nJet(); ijet++) {
      if (nano.Jet_pt()[ijet] <= min_jet_pt || fabs(nano.Jet_eta()[ijet]) > max_jet_eta) continue;
      ordered_idx.push_back(make_pair(ijet, nano.Jet_btagDeepB()[ijet]));
    }

    // enough jets to make two higgses?
    if (ordered_idx.size()<4) continue;

    sort(ordered_idx.begin(), ordered_idx.end(), 
          [](const pair<int, float> &a, const pair<int, float> &b) -> bool {
            return a.second > b.second;
          });


    int max_idx = ordered_idx.size()<=5 ? ordered_idx.size() : 5;
    vector<ROOT::Math::PtEtaPhiMVector> jets_lv;
    for (int idx(0); idx<max_idx; idx++) {
      int ijet = ordered_idx[idx].first;
      ROOT::Math::PtEtaPhiMVector lv(nano.Jet_pt()[ijet], nano.Jet_eta()[ijet], 
                                            nano.Jet_phi()[ijet], nano.Jet_mass()[ijet]);
      jets_lv.push_back(lv);
    }

    atto.out_njet() = 0; atto.out_nbl() = 0; atto.out_nbm() = 0; atto.out_nbt() = 0;
    for(int idx(0); idx<max_idx; ++idx){
      int ijet = ordered_idx[idx].first;

      atto.out_jet_pt().push_back(nano.Jet_pt()[ijet]);
      atto.out_jet_eta().push_back(nano.Jet_eta()[ijet]);
      atto.out_jet_phi().push_back(nano.Jet_phi()[ijet]);
      atto.out_jet_m().push_back(nano.Jet_mass()[ijet]);
      atto.out_jet_deepcsv().push_back(nano.Jet_btagDeepB()[ijet]);
      atto.out_jet_qgl().push_back(nano.Jet_qgl()[ijet]);

      atto.out_njet()++;
      if (nano.Jet_btagDeepB()[ijet] > btag_wpts[year][0]) atto.out_nbl()++; 
      if (nano.Jet_btagDeepB()[ijet] > btag_wpts[year][1]) atto.out_nbm()++; 
      if (nano.Jet_btagDeepB()[ijet] > btag_wpts[year][2]) atto.out_nbt()++;

      for(int jdx(idx+1); jdx<max_idx; ++jdx){
        int jjet = ordered_idx[jdx].first;
        ROOT::Math::PtEtaPhiMVector jj_lv = jets_lv[idx]+jets_lv[jdx];
        atto.out_mjj().push_back(jj_lv.M());
        atto.out_drjj().push_back(dR(nano.Jet_eta()[ijet], nano.Jet_eta()[jjet],nano.Jet_phi()[ijet], nano.Jet_phi()[jjet]));
      }
    } // end jet loop

    // if (max_idx==4) {
    //   atto.out_jet_pt().push_back(-9);
    //   atto.out_jet_eta().push_back(-9);
    //   atto.out_jet_phi().push_back(-9);
    //   atto.out_jet_m().push_back(-9);
    //   atto.out_jet_deepcsv().push_back(-9);
    //   atto.out_jet_qgl().push_back(-9);
    //   for (int i(0); i<4; i++) atto.out_mjj().push_back(-9);
    //   for (int i(0); i<4; i++) atto.out_drjj().push_back(-9);
    // }

    //make possible combinations with top 4 jets as we do in cut-based for reference
    vector<vector<unsigned>>  hcombs = {{0,1,2,3},{0,2,1,3},{0,3,1,2}};
    float dm(9999.), am(-9999.);
    for (unsigned ic(0); ic<hcombs.size(); ic++){
      ROOT::Math::PtEtaPhiMVector h1 = jets_lv[hcombs[ic][0]] + jets_lv[hcombs[ic][1]]; 
      ROOT::Math::PtEtaPhiMVector h2 = jets_lv[hcombs[ic][2]] + jets_lv[hcombs[ic][3]]; 

      float idm = fabs(h1.M()-h2.M());
      float iam = (h1.M()+h2.M())/2.;
      if (idm < dm) {
        dm = idm;
        atto.out_hig_dm() = idm;
        atto.out_hig_am() = iam;
      }
    }     

    atto.out_nbacc() = 0;
    for (int imc(0); imc<nano.nGenPart(); imc++) {
      if (nano.GenPart_pdgId()[imc]==1000025) atto.out_mchi() = nano.GenPart_mass()[imc];
      else if (nano.GenPart_pdgId()[imc]==1000022) atto.out_mlsp() = nano.GenPart_mass()[imc];
      else if (nano.GenPart_pdgId()[imc]==25) {
        atto.out_mhiggs() = nano.GenPart_mass()[imc];
        if (atto.out_h1_pt()<0) {
          atto.out_h1_pt() = nano.GenPart_pt()[imc];
        } else {
          if (nano.GenPart_pt()[imc] > atto.out_h1_pt()){
            atto.out_h2_pt() = atto.out_h1_pt();
            atto.out_h1_pt() = nano.GenPart_pt()[imc];
          } else {
            atto.out_h2_pt() = nano.GenPart_pt()[imc];
          }
        }
      } else if (abs(nano.GenPart_pdgId()[imc])==5 && 
                 nano.GenPart_pdgId()[nano.GenPart_genPartIdxMother()[imc]]==25 &&
                 fabs(nano.GenPart_eta()[imc]) <= max_jet_eta &&
                 nano.GenPart_pt()[imc] > min_jet_pt) {
        atto.out_nbacc()++;
      }
    }

    if (debug) cout<<"INFO:: Filling tree"<<endl;
    atto.Fill();
  } // loop over events

  atto.Write();

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
