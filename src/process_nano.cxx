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
#include "photon_producer.hpp"
#include "jet_producer.hpp"
#include "fjet_producer.hpp"
#include "hig_producer.hpp"
#include "zgamma_producer.hpp"

#include "btag_weighter.hpp"
#include "lepton_weighter.hpp"
#include "event_tools.hpp"

using namespace std;

namespace {
  string nano_file = "";
  string wgt_sums_file = "";
  string pico_file = "";
  bool isData = false;
  bool isFastsim = false;
  bool doSystematics = false;
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

  // B-tag working points
  map<int, vector<float>> btag_wpts{
    {2016, vector<float>({0.2217, 0.6321, 0.8953})},
    {2017, vector<float>({0.1522, 0.4941, 0.8001})},
    {2018, vector<float>({0.1241, 0.4184, 0.7527})}
  };
  map<int, vector<float>> btag_df_wpts{
    {2016, vector<float>({0.0614, 0.3093, 0.7221})},
    {2017, vector<float>({0.0521, 0.3033, 0.7489})},
    {2018, vector<float>({0.0494, 0.2770, 0.7264})}
  };

  //Initialize object producers
  GenParticleProducer mc_producer(year);
  ElectronProducer el_producer(year);
  MuonProducer mu_producer(year);
  DileptonProducer dilep_producer(year);
  IsoTrackProducer tk_producer(year);
  PhotonProducer photon_producer(year);
  JetProducer jet_producer(year);
  FatJetProducer fjet_producer(year);
  HigVarProducer hig_producer(year);
  ZGammaVarProducer zgamma_producer(year);

  //Initialize scale factor tools
  const string ctr = "central";
  const vector<string> updn = {"up","down"};
  const vector<BTagEntry::OperatingPoint> op_all = {BTagEntry::OP_LOOSE, BTagEntry::OP_MEDIUM, BTagEntry::OP_TIGHT};
  BTagWeighter btag_weighter(year, isFastsim, false, btag_wpts[year]);
  BTagWeighter btag_df_weighter(year, isFastsim, true, btag_df_wpts[year]);
  LeptonWeighter lep_weighter(year);

  // Event tools
  EventTools event_tools(nano_file, year);
  int event_type = event_tools.GetEventType();

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

  corr.out_neff() = 0;
  corr.out_nent_zlep() = 0;
  corr.out_tot_weight_l0() = 0.;
  corr.out_tot_weight_l1() = 0.;
  corr.out_nent() = nentries;

  for(size_t entry(0); entry<nentries; ++entry){
    nano.GetEntry(entry);
    if (entry%1000==0 || entry == nentries-1) {
      cout<<"Processing event: "<<entry<<endl;
    }

    // event info
    pico.out_event() = nano.event();
    pico.out_lumiblock() = nano.luminosityBlock();
    pico.out_run() = nano.run();
    pico.out_type() = event_type;
    pico.out_stitch() = isData ? true: event_tools.GetStitch(nano);
    // number of reconstructed primary vertices
    pico.out_npv() = nano.PV_npvs();

    // ----------------------------------------------------------------------------------------------
    //            *** Writing physics objects ***
    // N.B. Order in which producers are called matters! E.g. jets are not counted if overlapping 
    // with signal lepton, thus jets must be processed only after leptons have been selected.
    //-----------------------------------------------------------------------------------------------
    mc_producer.WriteGenParticles(nano, pico);
    mc_producer.WriteGenInfo(nano, pico);

    vector<int> jet_islep_nano_idx = vector<int>();
    vector<int> sig_el_nano_idx = el_producer.WriteElectrons(nano, pico, jet_islep_nano_idx);
    vector<int> sig_mu_nano_idx = mu_producer.WriteMuons(nano, pico, jet_islep_nano_idx);
    pico.out_nlep() = sig_el_nano_idx.size() + sig_mu_nano_idx.size();

    vector<int> jet_isphoton_nano_idx = vector<int>();
    vector<int> sig_ph_nano_idx = photon_producer.WritePhotons(nano, pico, jet_isphoton_nano_idx);

    tk_producer.WriteIsoTracks(nano, pico, sig_el_nano_idx, sig_mu_nano_idx);

    dilep_producer.WriteDielectrons(nano, pico, sig_el_nano_idx);
    dilep_producer.WriteDimuons(nano, pico, sig_mu_nano_idx);

    vector<int> sig_jet_nano_idx = jet_producer.WriteJets(nano, pico, jet_islep_nano_idx, 
                                                          btag_wpts[year], btag_df_wpts[year]);
    jet_producer.WriteJetSys(nano, pico, sig_jet_nano_idx, btag_wpts[year][1]); // usually w.r.t. medium WP
    fjet_producer.WriteFatJets(nano, pico);

    // Copy MET directly from NanoAOD
    pico.out_met() = nano.MET_pt();
    pico.out_met_phi() = nano.MET_phi();
    pico.out_met_calo() = nano.CaloMET_pt();
    pico.out_met_tru() = nano.GenMET_pt();
    pico.out_met_tru_phi() = nano.GenMET_phi();
 
    // calculate mT only for single lepton events
    pico.out_mt() = -999; 
    if (pico.out_nlep()==1) {
      if (sig_el_nano_idx.size()>0) {
        pico.out_mt() = GetMT(nano.MET_pt(), nano.MET_phi(), 
          nano.Electron_pt()[sig_el_nano_idx[0]], nano.Electron_phi()[sig_el_nano_idx[0]]);
      } else {
        pico.out_mt() = GetMT(nano.MET_pt(), nano.MET_phi(), 
          nano.Muon_pt()[sig_mu_nano_idx[0]], nano.Muon_phi()[sig_mu_nano_idx[0]]);
      }
    } 
    if (pico.out_ntrulep()==1) {
      for (unsigned imc(0); imc<pico.out_mc_id().size(); imc++){
        if (abs(pico.out_mc_id()[imc])==11 || abs(pico.out_mc_id()[imc])==13) {
          pico.out_mt_tru() = GetMT(pico.out_met_tru(), pico.out_met_tru_phi(), 
                                pico.out_mc_pt()[imc], pico.out_mc_phi()[imc]);
          break;
        }
      }
    } 

    // might need as input sig_el_nano_idx, sig_mu_nano_idx, sig_ph_nano_idx
    zgamma_producer.WriteZGammaVars();

    //save higgs variables using DeepCSV and DeepFlavor
    hig_producer.WriteHigVars(pico, /*DeepFlavor*/ false);
    hig_producer.WriteHigVars(pico, true);
    hig_producer.WriteDPhiVars();

    // N.B. Jets: pico.out_pass_jets() and pico.out_pass_fsjets() filled in jet_producer
    event_tools.WriteDataQualityFilters(nano, pico, sig_jet_nano_idx, isData, isFastsim);
    
    event_tools.CopyTriggerDecisions(nano, pico);

    // ----------------------------------------------------------------------------------------------
    //              *** Calculating weight branches ***
    // ----------------------------------------------------------------------------------------------

    float w_lep(1.), w_fs_lep(1.);
    vector<float> sys_lep(2,1.), sys_fs_lep(2,1.);
    lep_weighter.FullSim(pico, w_lep, sys_lep);
    if(isFastsim) lep_weighter.FastSim(pico, w_fs_lep, sys_fs_lep);
    pico.out_w_lep() = w_lep;
    pico.out_w_fs_lep() = w_fs_lep;
    if (doSystematics) { // don't worry about systematics just yet
      pico.out_sys_lep() = sys_lep; 
      pico.out_sys_fs_lep() = sys_fs_lep;
    }

    pico.out_w_btag() = btag_weighter.EventWeight(pico, BTagEntry::OP_MEDIUM, ctr, ctr);
    pico.out_w_btag_df() = btag_df_weighter.EventWeight(pico, BTagEntry::OP_MEDIUM, ctr, ctr);
    pico.out_w_bhig() = btag_weighter.EventWeight(pico, op_all, ctr, ctr);
    pico.out_w_bhig_df() = btag_df_weighter.EventWeight(pico, op_all, ctr, ctr);
    if (doSystematics) { // don't worry about systematics just yet
      for(size_t i = 0; i<2; ++i){ 
        pico.out_sys_bctag()[i] = btag_weighter.EventWeight(pico, BTagEntry::OP_MEDIUM, updn[i], ctr);
        pico.out_sys_bchig()[i] = btag_weighter.EventWeight(pico, op_all, updn[i], ctr);
        pico.out_sys_udsgtag()[i] = btag_weighter.EventWeight(pico, BTagEntry::OP_MEDIUM, ctr, updn[i]);
        pico.out_sys_udsghig()[i] = btag_weighter.EventWeight(pico, op_all, ctr, updn[i]);
      }
    }

    // will be assigned in subsequent step? 
    int neff_dataset(100e3); // @todo, if want to run in 1 step, need to get neff for full sample in advance
    double xsec(0.); 
    if (Contains(nano_file, "TChiHH")){
      double exsec(0.);
      int mchi = 100;//GetGluinoMass(nano_file); @todo
      xsec::higgsinoCrossSection(mchi, xsec, exsec);
    }else{
      xsec = xsec::crossSection(nano_file, (year==2016));  
    }
    const float lumi = 1000.; // calc. weight for 1fb-1 of total lumi
    pico.out_w_lumi() = xsec*lumi/neff_dataset; 

    // @todo, copy weights from babymaker
    pico.out_w_pu() = 1.;

    // pico.out_w_isr() = ;

    // N.B. out_w_prefire should not be renormalized because it models an inefficiency, 
    // i.e. we SHOULD get less events!
    // This is on the @todo list, to be implemented based on: 
    // https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/common/PrefireCorr.py
    pico.out_w_prefire() = 1.;

    pico.out_weight() = nano.Generator_weight() * pico.out_w_lumi() * w_lep *
                        (isFastsim ? w_fs_lep : 1.) * pico.out_w_bhig() * pico.out_w_pu();
    if (year==2016 || isFastsim)
      pico.out_weight() *= pico.out_w_isr();

    // ----------------------------------------------------------------------------------------------
    //              *** Add up weights to save for renormalization step ***
    // ----------------------------------------------------------------------------------------------

    // taking care of samples with negative weights
    corr.out_neff() += nano.Generator_weight()>0 ? 1:-1;

    // leptons, keeping track of 0l and 1l totals separately to determine the SF for 0l events
    if(pico.out_nlep()==0){
      corr.out_nent_zlep() += 1.;
      corr.out_tot_weight_l0() += pico.out_weight();
    }else{
      corr.out_tot_weight_l1() += pico.out_weight();
      corr.out_w_lep() += w_lep;
      if(isFastsim) corr.out_w_fs_lep() += w_fs_lep;
      for(size_t i = 0; i<pico.out_sys_lep().size(); ++i){
        corr.out_sys_lep().at(i) += sys_lep.at(i);
        if(isFastsim) corr.out_sys_fs_lep().at(i) += sys_fs_lep.at(i);
      }
    }

    // b-tag weights
    corr.out_w_btag() += pico.out_w_btag();
    corr.out_w_btag_df() += pico.out_w_btag_df();
    corr.out_w_bhig() += pico.out_w_bhig();
    corr.out_w_bhig_df() += pico.out_w_bhig_df();
    if (doSystematics) { // don't worry about systematics just yet
      for(size_t i = 0; i<2; ++i){ 
        corr.out_sys_bctag()[i] += pico.out_sys_bctag()[i];
        corr.out_sys_bchig()[i] += pico.out_sys_bchig()[i];
        corr.out_sys_udsgtag()[i] += pico.out_sys_udsgtag()[i];
        corr.out_sys_udsghig()[i] += pico.out_sys_udsghig()[i];
      }
    }

    // general event vars
    corr.out_weight() += pico.out_weight();
    corr.out_w_isr() += pico.out_w_isr();
    corr.out_w_pu() += pico.out_w_pu();

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
