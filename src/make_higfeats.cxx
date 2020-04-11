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
  string in_file_path = ""; //will be used insted of in_file and in_dir if specified
  string in_file = "merged_pico_higloose_met150_SMS-TChiHH_mChi-400_mLSP-0_higmc_higloose_nfiles_1.root";
  string in_dir = "/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/SMS-TChiHH_2D/merged_higmc_higloose/";
  string out_dir = "out";
  int nent_test = -1;
  bool debug = false;
}

int comboEncoder(pair<int, int> h1_jets, pair<int, int> h2_jets);
void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches       
  GetOptions(argc, argv);

  if((in_file_path=="" && (in_file=="" || in_dir=="")) || out_dir == "") {
    cout<<"ERROR: Input file and/or output directory not specified. Exit."<<endl;
    exit(0);
  }

  string out_file_path = "";
  if (in_file_path!="") {
    size_t found = in_file_path.find_last_of("/");
    out_file_path = out_dir+"/higfeats_"+in_file_path.substr(found+1);
  } else {
    in_file_path = in_dir+"/"+in_file;
    out_file_path = out_dir+"/higfeats_"+in_file;
  }

  time_t begtime, endtime;
  time(&begtime);

 // Initialize trees
  pico_tree pico(in_file_path);
  size_t nentries(nent_test>0 ? nent_test : pico.GetEntries());
  cout << "Pico input file: " << in_file_path << endl;
  cout << "Input number of events: " << nentries << endl;

  higfeats_tree higfeats("", out_file_path);
  cout << "Writing output to: " << out_file_path << endl;

  for(size_t entry(0); entry<nentries; ++entry){
    if (debug) cout << "GetEntry: " << entry << endl;
    pico.GetEntry(entry);
    if (entry%1000==0 || entry == nentries-1) {
      cout<<"Processing event: "<<entry<<endl;
    }

    if (pico.njet()<4) continue;

    if (debug) cout<<"INFO:: Filling event number."<<endl;
    higfeats.out_event() = pico.event();
    higfeats.out_nbt() = pico.nbt();
    higfeats.out_nbm() = pico.nbm();
    higfeats.out_nbl() = pico.nbl();

    higfeats.out_hig_cand_am() = -999.;
    if (pico.hig_cand_am().size()>0)
      higfeats.out_hig_cand_am() = pico.hig_cand_am()[0];

    //--------------------------------------------------------------
    //     Save AK4 jets ordered by deepCSV value    
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
    //    Save AK8 jets ordered by deepDoubleB value    
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

    //--------------------------------------------------------------
    //     Other branches of interest
    //--------------------------------------------------------------
    higfeats.out_mprod() = pico.mprod();
    higfeats.out_mlsp() = pico.mlsp();

    higfeats.out_nbacc() = 0;
    for (unsigned imc(0); imc<pico.mc_pt().size(); imc++) {
      if (pico.mc_id()[imc]==25) {
        // this will get triggered multiple times but that's ok, should be all the same
        higfeats.out_mhiggs() = pico.mc_mass()[imc];
      } else if (abs(pico.mc_id()[imc])==5 && pico.mc_mom()[imc]==25 &&
                  fabs(pico.mc_eta()[imc]) <= 2.4 && pico.mc_pt()[imc] > 30) {
        higfeats.out_nbacc()++;
      }
    }

    //--------------------------------------------------------------
    //     Match Higgs decay products to jets
    //--------------------------------------------------------------
    pair<int, int> h1_jets(make_pair(-1,-1)), h2_jets(make_pair(-1,-1));
    int h1_idx(-1), h2_idx(-1);
    set<int> selected_jets; 

    //printing jets
    // for (unsigned idx(0); idx<max_idx; idx++) {
    //   unsigned ijet = ordered_idx[idx].first;
    //   cout<<"Jet "<<ijet<<" ("<<pico.jet_pt()[ijet]<<","<<pico.jet_eta()[ijet]
    //                      <<","<<pico.jet_phi()[ijet]<<","<<pico.jet_m()[ijet]<<")"<<endl;
    // }

    for (unsigned imc(0); imc<pico.mc_pt().size(); imc++) {
      if (abs(pico.mc_id()[imc])==5 && pico.mc_mom()[imc]==25) {
        // cout<<"b "<<imc<<" ("<<pico.mc_pt()[imc]<<","<<pico.mc_eta()[imc]
        //                <<","<<pico.mc_phi()[imc]<<","<<pico.mc_mass()[imc]<<")"<<endl;
        // cout<<"found b (idx "<<imc<<","<<pico.mc_status()[imc]<<") with mom idx "<<pico.mc_momidx()[imc];
        float dr_min(9999);
        unsigned closest_jet(-1);
        // loop through original jets, since the output phi is transformed
        for (unsigned idx(0); idx<max_idx; idx++) {
          unsigned ijet = ordered_idx[idx].first;
          if (selected_jets.find(idx) != selected_jets.end()) continue;
          float dr_ = dR(pico.mc_eta()[imc], pico.jet_eta()[ijet], 
                         pico.mc_phi()[imc], pico.jet_phi()[ijet]);
          if (dr_ < dr_min) {
            dr_min = dr_;
            closest_jet = idx;
          }
        }
        if (dr_min<0.4) {
          // cout<<" matching with jet "<<closest_jet<<endl;
          selected_jets.insert(closest_jet);
          if (h1_idx==-1 || h1_idx==pico.mc_momidx()[imc]) {
            h1_idx = pico.mc_momidx()[imc];
            if (h1_jets.first==-1) h1_jets.first = closest_jet;
            else if (h1_jets.second==-1) h1_jets.second = closest_jet;
            else cout<<"ERROR: Found more than two decay products for the Higgs"<<endl;
          } else if (h2_idx==-1 || h2_idx==pico.mc_momidx()[imc]) {
            h2_idx = pico.mc_momidx()[imc];
            if (h2_jets.first==-1) h2_jets.first = closest_jet;
            else if (h2_jets.second==-1) h2_jets.second = closest_jet;
            else cout<<"ERROR: Found more than two decay products for the Higgs"<<endl;
          } else {
            cout<<"ERROR: Found more than two Higgs particles"<<endl;
          }
        } 
      }
    }
     
    higfeats.out_combo_code() = (selected_jets.size()==4) ? comboEncoder(h1_jets, h2_jets) : -1;
    // cout<<"INFO: Combination code is: "<<higfeats.out_combo_code()<<endl;

    if (debug) cout<<"INFO:: Filling tree"<<endl;
    higfeats.Fill();
  } // loop over events

  higfeats.Write();

  cout<<endl;
  time(&endtime); 
  cout<<"Time passed: "<<hoursMinSec(difftime(endtime, begtime))<<endl<<endl;  
}

int comboEncoder(pair<int, int> h1_jets, pair<int, int> h2_jets) {
  // order pairs to enable encoding
  array<int,4> a;
  if (h1_jets.first > h1_jets.second) h1_jets = std::make_pair(h1_jets.second, h1_jets.first);
  if (h2_jets.first > h2_jets.second) h2_jets = std::make_pair(h2_jets.second, h2_jets.first);
  if (h2_jets.first > h1_jets.first) {
    a = {h1_jets.first, h1_jets.second, h2_jets.first, h2_jets.second};
  } else {
    a = {h2_jets.first, h2_jets.second, h1_jets.first, h1_jets.second};
  } 

  // 5 combinations of 4 jets - each has 3 possible pairings
  // combo 0,1,2,3
  if      (a[0]==0 && a[1]==1 && a[2]==2 && a[3]==3) return 0; //pair option 0
  else if (a[0]==0 && a[1]==2 && a[2]==1 && a[3]==3) return 1; //pair option 1
  else if (a[0]==0 && a[1]==3 && a[2]==1 && a[3]==2) return 2; //pair option 2
  // combo 0,1,2,4
  else if (a[0]==0 && a[1]==1 && a[2]==2 && a[3]==4) return 3;
  else if (a[0]==0 && a[1]==2 && a[2]==1 && a[3]==4) return 4;
  else if (a[0]==0 && a[1]==4 && a[2]==1 && a[3]==2) return 5;
  // combo 0,1,3,4
  else if (a[0]==0 && a[1]==1 && a[2]==3 && a[3]==4) return 6;
  else if (a[0]==0 && a[1]==3 && a[2]==1 && a[3]==4) return 7;
  else if (a[0]==0 && a[1]==4 && a[2]==1 && a[3]==3) return 8;
  // combo 0,2,3,4
  else if (a[0]==0 && a[1]==2 && a[2]==3 && a[3]==4) return 9;
  else if (a[0]==0 && a[1]==3 && a[2]==2 && a[3]==4) return 10;
  else if (a[0]==0 && a[1]==4 && a[2]==2 && a[3]==3) return 11;
  // combo 1,2,3,4
  else if (a[0]==1 && a[1]==2 && a[2]==3 && a[3]==4) return 12;
  else if (a[0]==1 && a[1]==3 && a[2]==2 && a[3]==4) return 13;
  else if (a[0]==1 && a[1]==4 && a[2]==2 && a[3]==3) return 14;
  else {
    std::cout<<"ERROR: Could not find matching combination for ("<<a[0]<<","<<a[1]<<") ("<<a[2]<<","<<a[3]<<")"<<endl;
    return -1;
  }
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"in_file_path", required_argument, 0,'p'}, //alternative way to specify input, more convenient sometime
      {"in_file", required_argument, 0,'f'}, 
      {"in_dir",  required_argument, 0,'i'},
      {"out_dir", required_argument, 0,'o'},
      {"nent",    required_argument, 0, 0},
      {"debug",    no_argument, 0, 'd'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "p:f:i:o:d", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'p':
      in_file_path = optarg;
      break;
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
