#include <ctime>

#include <iostream>
#include <iomanip>

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
#include "met_producer.hpp"
#include "hig_producer.hpp"
#include "zgamma_producer.hpp"
#include "in_json.hpp"

#include "btag_weighter.hpp"
#include "lepton_weighter.hpp"
#include "prefire_weighter.hpp"
#include "photon_weighter.hpp"
#include "event_tools.hpp"
#include "isr_tools.hpp"

using namespace std;

namespace {
  string in_file = "";
  string in_dir = "";
  string out_dir = "";
  int nent_test = -1;
  bool debug = false;
  // requirements for jets to be counted in njet, mofified for Zgamma below
  float min_jet_pt = 30.0;
  float max_jet_eta =  2.4;
}

void WriteDataQualityFilters(nano_tree& nano, pico_tree& pico);
void CopyTriggerDecisions(nano_tree& nano, pico_tree& pico);
void Initialize(corrections_tree& wgt_sums);
void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  // gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches       
  GetOptions(argc, argv);

  if(in_file=="" || in_dir=="" || out_dir == "") {
    cout<<"ERROR: Input file, sum-of-weights and/or output directory not specified. Exit."<<endl;
    exit(0);
  }

  bool isData = Contains(in_file, "Run201") ? true : false;
  bool isFastsim = Contains(in_file, "Fast") ? true : false;
  int year = Contains(in_file, "RunIISummer16") ? 2016 : (Contains(in_file, "RunIIFall17") ? 2017 : 2018);
  if (isData) {
    year = Contains(in_file, "Run2016") ? 2016 : (Contains(in_file, "Run2017") ? 2017 : 2018);
  }

  vector<vector<int>> VVRunLumi;
  if (isData) {
    switch (year) {
      case 2016:
        VVRunLumi = MakeVRunLumi("golden2016");
        break;
      case 2017:
        VVRunLumi = MakeVRunLumi("golden2017");
        break;
      case 2018:
        VVRunLumi = MakeVRunLumi("golden2018");
        break;
      default:
        cout << "ERROR: no golden cert for given year" << endl;
        return 0;
    }
  }

  string in_path = in_dir+"/"+in_file;
  string wgt_sums_path = out_dir+"/wgt_sums/wgt_sums_"+in_file;
  string out_path;
  out_path = out_dir+"/raw_pico/raw_pico_"+in_file;

  bool isZgamma = Contains(out_dir, "zgamma");

  time_t begtime, endtime;
  time(&begtime);

  // jet requirements
  if(isZgamma) max_jet_eta  =  4.7;

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
  ElectronProducer el_producer(year, isData);
  MuonProducer mu_producer(year, isData);
  DileptonProducer dilep_producer(year);
  IsoTrackProducer tk_producer(year);
  PhotonProducer photon_producer(year);
  JetProducer jet_producer(year, min_jet_pt, max_jet_eta, isData);
  MetProducer met_producer(year, isData);
  HigVarProducer hig_producer(year);
  ZGammaVarProducer zgamma_producer(year);

  //Initialize scale factor tools
  const string ctr = "central";
  const vector<string> updn = {"up","down"};
  const vector<BTagEntry::OperatingPoint> op_all = {BTagEntry::OP_LOOSE, BTagEntry::OP_MEDIUM, BTagEntry::OP_TIGHT};
  BTagWeighter btag_weighter(year, false, false, btag_wpts[year]); // not applying the FastSim scale factors for now since they seem to have NaN's...
  BTagWeighter btag_df_weighter(year, false, true, btag_df_wpts[year]);
  // BTagWeighter btag_weighter(year, isFastsim, false, btag_wpts[year]);
  // BTagWeighter btag_df_weighter(year, isFastsim, true, btag_df_wpts[year]);
  LeptonWeighter lep_weighter(year, isZgamma);
  LeptonWeighter lep_weighter16gh(year, isZgamma, true);
  PrefireWeighter prefire_weighter(year, true);
  PhotonWeighter photon_weighter(year, isZgamma);

  // Other tools
  EventTools event_tools(in_path, year);
  int event_type = event_tools.GetEventType();

  ISRTools isr_tools(in_path, year);

  // Initialize trees
  nano_tree nano(in_path);
  size_t nentries(nent_test>0 ? nent_test : nano.GetEntries());
  cout << "Nano file: " << in_path << endl;
  cout << "Input number of events: " << nentries << endl;
  // cout << "Running on "<< (isFastsim ? "FastSim" : "FullSim") << endl;
  // cout << "Calculating weights based on " << year << " scale factors." << endl;

  pico_tree pico("", out_path);
  cout << "Writing output to: " << out_path << endl;

  corrections_tree wgt_sums("", wgt_sums_path);
  cout << "Writing sum-of-weights to: " << wgt_sums_path << endl;
  Initialize(wgt_sums);
  wgt_sums.out_nent() = nentries;

  for(size_t entry(0); entry<nentries; ++entry){
    if (debug) cout << "GetEntry: " << entry <<" event = "<< endl;
    nano.GetEntry(entry);
    if (entry%1000==0 || entry == nentries-1) {
      cout<<"Processing event: "<<entry<<endl;
    }

    //skip events that are data but not in the golden json
    if (isData) {
      if(!inJSON(VVRunLumi, nano.run(), nano.luminosityBlock())) continue; 
    }

    bool passed_trig = event_tools.SaveTriggerDecisions(nano, pico, isZgamma);
    if (isData && !passed_trig) continue;

    // event info
    pico.out_event()     = nano.event();
    pico.out_lumiblock() = nano.luminosityBlock();
    pico.out_run()       = nano.run();
    pico.out_type()      = event_type;
    if (isData) pico.out_stitch() = true;
    else event_tools.WriteStitch(nano, pico);
    // number of reconstructed primary vertices
    pico.out_npv() = nano.PV_npvs();
    pico.out_npv_good() = nano.PV_npvsGood();
    // number of pileup in mc
    if (!isData) {
      pico.out_npu_tru() = nano.Pileup_nPU();
      pico.out_npu_tru_mean() = nano.Pileup_nTrueInt();
    }

    // ----------------------------------------------------------------------------------------------
    //            *** Writing physics objects ***
    // N.B. Order in which producers are called matters! E.g. jets are not counted if overlapping 
    // with signal lepton, thus jets must be processed only after leptons have been selected.
    //-----------------------------------------------------------------------------------------------
    if (debug) cout<<"INFO:: Writing leptons, photons and tracks"<<endl;
    vector<int> jet_islep_nano_idx = vector<int>();
    vector<int> jet_isvlep_nano_idx = vector<int>();
    pico.out_nlep() = 0; pico.out_nvlep() = 0; // filled by lepton producers
    vector<int> sig_el_pico_idx = vector<int>();
    vector<int> sig_mu_pico_idx = vector<int>();
    vector<int> sig_el_nano_idx = el_producer.WriteElectrons(nano, pico, jet_islep_nano_idx, jet_isvlep_nano_idx, sig_el_pico_idx, isZgamma);
    vector<int> sig_mu_nano_idx = mu_producer.WriteMuons(nano, pico, jet_islep_nano_idx, jet_isvlep_nano_idx, sig_mu_pico_idx, isZgamma);

    // save a separate vector with just signal leptons ordered by pt
    struct SignalLepton{ float pt; float eta; float phi; int pdgid;};
    vector<SignalLepton> sig_leps;
    for (auto &iel: sig_el_nano_idx)
      sig_leps.push_back({nano.Electron_pt()[iel], nano.Electron_eta()[iel], 
                          nano.Electron_phi()[iel], nano.Electron_pdgId()[iel]});
    for (auto &imu: sig_mu_nano_idx) 
      sig_leps.push_back({nano.Muon_pt()[imu], nano.Muon_eta()[imu], 
                          nano.Muon_phi()[imu], nano.Muon_pdgId()[imu]});

    auto greaterPt = [](SignalLepton lep1, SignalLepton lep2){ return lep1.pt > lep2.pt;};
    sort(sig_leps.begin(), sig_leps.end(), greaterPt);
    for(auto &ilep : sig_leps) {
      pico.out_lep_pt().push_back(ilep.pt);
      pico.out_lep_eta().push_back(ilep.eta);
      pico.out_lep_phi().push_back(ilep.phi);
      pico.out_lep_pdgid().push_back(ilep.pdgid);
    }

    vector<int> jet_isphoton_nano_idx = vector<int>();
    if(isZgamma) 
      vector<int> sig_ph_nano_idx = photon_producer.WritePhotons(nano, pico, jet_isphoton_nano_idx,
                                                                 sig_el_nano_idx, sig_mu_nano_idx);

    tk_producer.WriteIsoTracks(nano, pico, sig_el_nano_idx, sig_mu_nano_idx);

    dilep_producer.WriteDileptons(pico, sig_el_pico_idx, sig_mu_pico_idx);

    if (debug) cout<<"INFO:: Writing gen particles"<<endl;

    pico.out_stitch_dy() = true;
    if (!isData)
	    mc_producer.WriteGenParticles(nano, pico);
    isr_tools.WriteISRSystemPt(nano, pico);

    if (debug) cout<<"INFO:: Writing jets, MET and ISR vars"<<endl;
    vector<int> sig_jet_nano_idx = jet_producer.WriteJets(nano, pico, jet_islep_nano_idx, jet_isvlep_nano_idx, jet_isphoton_nano_idx,
                                                          btag_wpts[year], btag_df_wpts[year]);
    jet_producer.WriteJetSystemPt(nano, pico, sig_jet_nano_idx, btag_wpts[year][1]); // usually w.r.t. medium WP
    if(!isZgamma){
      jet_producer.WriteFatJets(nano, pico); // jet_producer.SetVerbose(nano.nSubJet()>0);
      jet_producer.WriteSubJets(nano, pico);
    }
    isr_tools.WriteISRJetMultiplicity(nano, pico);

    met_producer.WriteMet(nano, pico);
    // Copy MET and ME ISR directly from NanoAOD
    //pico.out_met()         = nano.MET_pt();
    //pico.out_met_phi()     = nano.MET_phi();
    //pico.out_met_calo()    = nano.CaloMET_pt();
    //pico.out_met_tru()     = nano.GenMET_pt();
    //pico.out_met_tru_phi() = nano.GenMET_phi();
    //pico.out_ht_isr_me()   = nano.LHE_HTIncoming();
 
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

    if (debug) cout<<"INFO:: Writing analysis specific variables"<<endl;
    // might need as input sig_el_nano_idx, sig_mu_nano_idx, sig_ph_nano_idx
    if(isZgamma)
      zgamma_producer.WriteZGammaVars(nano, pico, sig_jet_nano_idx);

    //save higgs variables using DeepCSV and DeepFlavor
    hig_producer.WriteHigVars(pico, false);
    hig_producer.WriteHigVars(pico, true);

    if (debug) cout<<"INFO:: Writing filters and triggers"<<endl;
    // N.B. Jets: pico.out_pass_jets() and pico.out_pass_fsjets() filled in jet_producer
    event_tools.WriteDataQualityFilters(nano, pico, sig_jet_nano_idx, min_jet_pt, isData, isFastsim);

    event_tools.WriteTriggerEfficiency(pico);

    // ----------------------------------------------------------------------------------------------
    //              *** Calculating weight branches ***
    // ----------------------------------------------------------------------------------------------
    if (debug) cout<<"INFO:: Calculating weights"<<endl;
    float w_lep(1.), w_fs_lep(1.);
    float w_photon(1.);
    vector<float> sys_lep(2,1.), sys_fs_lep(2,1.);
    vector<float> sys_photon(2,1.);
    if(isZgamma) {
      photon_weighter.FullSim(pico, w_photon, sys_photon);
      if(nano.event() % 3516 <= 1887) lep_weighter.FullSim(pico, w_lep, sys_lep);
      else lep_weighter16gh.FullSim(pico, w_lep, sys_lep);
    }
    else if (!isData) {
      lep_weighter.FullSim(pico, w_lep, sys_lep);
    }
    if(isFastsim) lep_weighter.FastSim(pico, w_fs_lep, sys_fs_lep);
    pico.out_w_lep()      = w_lep;
    pico.out_w_fs_lep()   = w_fs_lep;
    pico.out_sys_lep()    = sys_lep; 
    pico.out_sys_fs_lep() = sys_fs_lep;
    pico.out_w_photon()   = w_photon;
    pico.out_sys_photon() = sys_photon; 
    if(isZgamma || isData) {
      pico.out_w_btag()    = 1.; 
      pico.out_w_btag_df() = 1.; 
      pico.out_w_bhig()    = 1.; 
      pico.out_w_bhig_df() = 1.; 
      pico.out_sys_bchig().resize(2,0); pico.out_sys_udsghig().resize(2,0);
      pico.out_sys_fs_bchig().resize(2,0); pico.out_sys_fs_udsghig().resize(2,0);
    }
    else {
      pico.out_w_btag()    = btag_weighter.EventWeight(pico, BTagEntry::OP_MEDIUM, ctr, ctr);
      pico.out_w_btag_df() = btag_df_weighter.EventWeight(pico, BTagEntry::OP_MEDIUM, ctr, ctr);
      pico.out_w_bhig()    = btag_weighter.EventWeight(pico, op_all, ctr, ctr);
      pico.out_w_bhig_df() = btag_df_weighter.EventWeight(pico, op_all, ctr, ctr);
      pico.out_sys_bchig().resize(2,0); pico.out_sys_udsghig().resize(2,0);
      pico.out_sys_fs_bchig().resize(2,0); pico.out_sys_fs_udsghig().resize(2,0);
      for(size_t i = 0; i<2; ++i){ 
        pico.out_sys_bchig()[i]   = btag_weighter.EventWeight(pico, op_all, updn[i], ctr);
        pico.out_sys_udsghig()[i] = btag_weighter.EventWeight(pico, op_all, ctr, updn[i]);
        if (isFastsim) {
          pico.out_sys_fs_bchig()[i]   = btag_weighter.EventWeight(pico, op_all, ctr, ctr, updn[i], ctr);
          pico.out_sys_fs_udsghig()[i] = btag_weighter.EventWeight(pico, op_all, ctr, ctr, ctr, updn[i]);
        }
      }
    }

    // to be calculated in Step 2: merge_corrections
    if (!isData)
      pico.out_w_lumi() = nano.Generator_weight()>0 ? 1:-1;
    else
      pico.out_w_lumi() = 1.;

    // @todo, copy weights from babymaker
    pico.out_w_pu() = 1.;
    pico.out_sys_pu().resize(2, 1.);

    isr_tools.WriteISRWeights(pico);

    // N.B. out_w_prefire should not be renormalized because it models an inefficiency, 
    // i.e. we SHOULD get less events!
    float w_prefire=1.;
    std::vector<float> sys_prefire(2, 1.);
    prefire_weighter.EventWeight(nano, w_prefire, sys_prefire);
    pico.out_w_prefire() = w_prefire;
    pico.out_sys_prefire() = sys_prefire;

    // do not include w_prefire, or anything that should not be renormalized! Will be set again in Step 3
    pico.out_weight() = pico.out_w_lumi() *
                        w_lep * w_fs_lep * pico.out_w_bhig() *
                        w_photon  *
                        pico.out_w_isr() * pico.out_w_pu();

    // ----------------------------------------------------------------------------------------------
    //              *** Add up weights to save for renormalization step ***
    // ----------------------------------------------------------------------------------------------
    if (debug) cout<<"INFO:: Writing sum of weights"<<endl;
    if (!isData) {
      wgt_sums.out_weight() += pico.out_weight();
      // taking care of samples with negative weights
      wgt_sums.out_neff() += nano.Generator_weight()>0 ? 1:-1;

      // leptons, keeping track of 0l and 1l totals separately to determine the SF for 0l events
      if(pico.out_nlep()==0){
        wgt_sums.out_nent_zlep() += 1.;
        wgt_sums.out_tot_weight_l0() += pico.out_weight()*(nano.Generator_weight()>0 ? 1:-1); // multiplying by GenWeight to remove the sign...
      }else{
        wgt_sums.out_tot_weight_l1() += pico.out_weight()*(nano.Generator_weight()>0 ? 1:-1);
        wgt_sums.out_w_lep() += w_lep;
        if(isFastsim) wgt_sums.out_w_fs_lep() += w_fs_lep;
        for(size_t i = 0; i<pico.out_sys_lep().size(); ++i){
          wgt_sums.out_sys_lep()[i] += sys_lep[i];
          wgt_sums.out_sys_fs_lep()[i] += sys_fs_lep[i];
        }
      }
      if(pico.out_nphoton()>0) {
        wgt_sums.out_w_photon() += w_photon;
        for(size_t i = 0; i < pico.out_sys_photon().size(); ++i) {
          wgt_sums.out_sys_photon()[i] += sys_photon[i];
        }
      }
      wgt_sums.out_w_btag()    += pico.out_w_btag();
      wgt_sums.out_w_btag_df() += pico.out_w_btag_df();
      wgt_sums.out_w_bhig()    += pico.out_w_bhig();
      wgt_sums.out_w_bhig_df() += pico.out_w_bhig_df();
      wgt_sums.out_w_isr()     += pico.out_w_isr();
      wgt_sums.out_w_pu()      += pico.out_w_pu();

      for(size_t i = 0; i<2; ++i){ 
        wgt_sums.out_sys_bchig()[i]      += pico.out_sys_bchig()[i];
        wgt_sums.out_sys_udsghig()[i]    += pico.out_sys_udsghig()[i];
        wgt_sums.out_sys_fs_bchig()[i]   += pico.out_sys_fs_bchig()[i];
        wgt_sums.out_sys_fs_udsghig()[i] += pico.out_sys_fs_udsghig()[i];
        wgt_sums.out_sys_isr()[i]        += pico.out_sys_isr()[i];
        wgt_sums.out_sys_pu()[i]         += pico.out_sys_pu()[i];
      }
    }

    if (debug) cout<<"INFO:: Filling tree"<<endl;
    pico.Fill();
  } // loop over events

  wgt_sums.Fill();
  wgt_sums.Write();
  pico.Write();

  cout<<endl;
  time(&endtime); 
  cout<<"Time passed: "<<hoursMinSec(difftime(endtime, begtime))<<endl<<endl;  
}

void Initialize(corrections_tree &wgt_sums){
  wgt_sums.out_neff()          = 0;
  wgt_sums.out_nent_zlep()     = 0;
  wgt_sums.out_tot_weight_l0() = 0.;
  wgt_sums.out_tot_weight_l1() = 0.;

  wgt_sums.out_weight()     = 0.;
  wgt_sums.out_w_lumi()     = 0.;
  wgt_sums.out_w_lep()      = 0.;
  wgt_sums.out_w_fs_lep()   = 0.;
  wgt_sums.out_w_photon()   = 0.;
  wgt_sums.out_w_btag()     = 0.;
  wgt_sums.out_w_btag_df()  = 0.;
  wgt_sums.out_w_bhig()     = 0.;
  wgt_sums.out_w_bhig_df()  = 0.;
  wgt_sums.out_w_isr()      = 0.;
  wgt_sums.out_w_pu()       = 0.;
  // w_prefire should not be normalized!!

  wgt_sums.out_sys_lep().resize(2,0);
  wgt_sums.out_sys_fs_lep().resize(2,0);
  wgt_sums.out_sys_photon().resize(2,0);
  wgt_sums.out_sys_bchig().resize(2,0);
  wgt_sums.out_sys_udsghig().resize(2,0);
  wgt_sums.out_sys_fs_bchig().resize(2,0);
  wgt_sums.out_sys_fs_udsghig().resize(2,0);
  wgt_sums.out_sys_isr().resize(2,0);
  wgt_sums.out_sys_pu().resize(2,0);
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
