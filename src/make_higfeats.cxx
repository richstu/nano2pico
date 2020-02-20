#include <ctime>

#include <iostream>
#include <iomanip>

#include <getopt.h>

#include "TError.h"
#include "Math/Vector4D.h"

#include "pico_tree.hpp"
#include "higfeats_tree.hpp"
#include "utilities.hpp"

using namespace std;

namespace {
  string in_file = "merged_pico_higloose_met150_SMS-TChiHH_mChi-400_mLSP-0_higmc_higloose_nfiles_1.root";
  string in_dir = "/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/SMS-TChiHH_2D/merged_higmc_higloose/";
  string out_dir = "out";
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
  string out_path = out_dir+"/higfeats_"+in_file;

  time_t begtime, endtime;
  time(&begtime);

 // Initialize trees
  pico_tree pico(in_path);
  size_t nentries(nent_test>0 ? nent_test : pico.GetEntries());
  cout << "Pico input file: " << in_path << endl;
  cout << "Input number of events: " << nentries << endl;

  higfeats_tree higfeats("", out_path);
  cout << "Writing output to: " << out_path << endl;

  for(size_t entry(0); entry<nentries; ++entry){
    if (debug) cout << "GetEntry: " << entry << endl;
    pico.GetEntry(entry);
    if (entry%1000==0 || entry == nentries-1) {
      cout<<"Processing event: "<<entry<<endl;
    }

    if (debug) cout<<"INFO:: Filling event number."<<endl;
    higfeats.out_event() = pico.event();

    //--------------------------------------------------------------
    //         Save AK4 jets ordered by deepCSV value    
    //--------------------------------------------------------------
    if (debug) cout<<"INFO:: Filling Ak4 jet higfeats."<<endl;
    vector<pair<unsigned, float>>  ordered_idx;
    for (unsigned ijet(0); ijet<pico.jet_pt().size(); ijet++) {
      if (fabs(pico.jet_eta()[ijet]) > 2.4) continue; // because this cut is not applied in pico!
      ordered_idx.push_back(make_pair(ijet, pico.jet_deepcsv()[ijet]));
    }

    higfeats.out_njet() = ordered_idx.size();

    sort(ordered_idx.begin(), ordered_idx.end(), 
          [](const pair<unsigned, float> &a, const pair<unsigned, float> &b) -> bool {
            return a.second > b.second;
          });

    vector<ROOT::Math::PtEtaPhiMVector> jets_lv;
    unsigned max_idx = ordered_idx.size()<=5 ? ordered_idx.size() : 5;
    for (unsigned idx(0); idx<max_idx; idx++) {
      unsigned ijet = ordered_idx[idx].first;
      higfeats.out_jet_brank_pt().push_back(pico.jet_pt()[ijet]);
      higfeats.out_jet_brank_eta().push_back(pico.jet_eta()[ijet]);
      higfeats.out_jet_brank_phi().push_back(SignedDeltaPhi(pico.jet_phi()[0], pico.jet_phi()[ijet]));
      // don't write negative masses since unphysical and it interferes with taking the log
      // only 5 events have small negative mass in training sample
      higfeats.out_jet_brank_m().push_back(pico.jet_m()[ijet]<=0 ? 1e-5 : pico.jet_m()[ijet]);
      higfeats.out_jet_brank_deepcsv().push_back(pico.jet_deepcsv()[ijet]);
      higfeats.out_jet_brank_qgl().push_back(pico.jet_qgl()[ijet]);

      ROOT::Math::PtEtaPhiMVector lv(pico.jet_pt()[ijet], pico.jet_eta()[ijet], 
                                            pico.jet_phi()[ijet], pico.jet_m()[ijet]);
      jets_lv.push_back(lv);
    }

    if (debug) cout<<"INFO:: Filling Ak4 jet pair higfeats."<<endl;
    for(unsigned idx(0); idx<max_idx; ++idx){
      unsigned ijet = ordered_idx[idx].first;
      for(unsigned jdx(idx+1); jdx<max_idx; ++jdx){
        unsigned jjet = ordered_idx[jdx].first;
        ROOT::Math::PtEtaPhiMVector jj_lv = jets_lv[idx]+jets_lv[jdx];
        higfeats.out_mjj().push_back(jj_lv.M());
        higfeats.out_drjj().push_back(dR(pico.jet_eta()[ijet], pico.jet_eta()[jjet],pico.jet_phi()[ijet], pico.jet_phi()[jjet]));
      }
    } // end jet loop

    //--------------------------------------------------------------
    //         Save AK8 jets ordered by deepDoubleB value    
    //--------------------------------------------------------------
    if (debug) cout<<"INFO:: Filling Ak8 jet higfeats."<<endl;
    vector<pair<unsigned, float>>  ordered_fjet_idx;
    for (unsigned ifjet(0); ifjet<pico.fjet_pt().size(); ifjet++) {
      if (pico.fjet_pt()[ifjet] <= 300. || fabs(pico.fjet_eta()[ifjet]) > 2.4) continue;
      ordered_fjet_idx.push_back(make_pair(ifjet, pico.fjet_deep_md_hbb_btv()[ifjet]));
    }

    sort(ordered_fjet_idx.begin(), ordered_fjet_idx.end(), 
          [](const pair<unsigned, float> &a, const pair<unsigned, float> &b) -> bool {
            return a.second > b.second;
          });

    unsigned max_fjet_idx = ordered_fjet_idx.size()<=2 ? ordered_fjet_idx.size() : 2;
    for(unsigned idx(0); idx<max_fjet_idx; ++idx){
      unsigned ifjet = ordered_fjet_idx[idx].first;

      higfeats.out_fjet_brank_pt().push_back(pico.fjet_pt()[ifjet]);
      higfeats.out_fjet_brank_eta().push_back(pico.fjet_eta()[ifjet]);
      higfeats.out_fjet_brank_phi().push_back(SignedDeltaPhi(pico.fjet_phi()[0], pico.fjet_phi()[ifjet]));
      // don't write negative masses since unphysical and it interferes with taking the log
      // only 5 events have small negative mass in training sample
      higfeats.out_fjet_brank_m().push_back(pico.fjet_m()[ifjet]<=0 ? 1e-5 : pico.fjet_m()[ifjet]);
      higfeats.out_fjet_brank_ddb().push_back(pico.fjet_deep_md_hbb_btv()[ifjet]);
    }

    //     Other branches of interest
    //------------------------------
    higfeats.out_mprod() = pico.mprod();
    higfeats.out_mlsp() = pico.mlsp();

    higfeats.out_nbacc() = 0;
    for (unsigned imc(0); imc<pico.mc_pt().size(); imc++) {
      if (pico.mc_id()[imc]==25) {
        // this will get triggered multiple times but that's ok, should be all the same
        higfeats.out_mhiggs() = pico.mc_mass()[imc];
      } else if (abs(pico.mc_id()[imc])==5 && 
                  pico.mc_id()[pico.mc_momidx()[imc]]==25 &&
                  fabs(pico.mc_eta()[imc]) <= 2.4 && 
                  pico.mc_pt()[imc] > 30) {
        higfeats.out_nbacc()++;
      }
    }

    if (debug) cout<<"INFO:: Filling tree"<<endl;
    higfeats.Fill();
  } // loop over events

  higfeats.Write();

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
