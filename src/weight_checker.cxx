/**
 * Script to validate weights and make plots
 */

#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "pico_tree.hpp"
#include "event_weighter.hpp"
#include "trigger_weighter.hpp"

#include "TCanvas.h"

using std::cin;
using std::cout;
using std::endl;
using std::isnan;
using std::isinf;
using std::map;
using std::string;
using std::to_string;
using std::vector;

const float sf_big_threshold = 12.0;

unsigned max(unsigned a, unsigned b) {
  if (a>b) return a;
  return b;
}

bool sf_is_bad(float value) {
  return isnan(value) || isinf(value) || (fabs(value)>sf_big_threshold);
}

void incr_vector(vector<unsigned>& vec, unsigned max) {
  for (unsigned ivec = 0; ivec < vec.size(); ivec++) {
    if (vec[ivec]<(max-1)) {
      vec[ivec] += 1;
      return;
    }
    else {
      vec[ivec] = 0;
      //carry by continuing loop
    }
  }
  //set vector to something that will terminate loop by stop_vector
  vec.push_back(1);
}

//divides both content and error of each bin in hist1 by hist2
void average_hist(TH1D &hist1, TH1D &hist2) {
  for (int ibin = 0; ibin < hist1.GetNbinsX()+2 ; ibin++) {
    double sum_w = hist1.GetBinContent(ibin);
    double sqrt_sum_w2 = hist1.GetBinError(ibin);
    double n = hist2.GetBinContent(ibin);
    hist1.SetBinError(ibin, sqrt(sqrt_sum_w2*sqrt_sum_w2/n-sum_w*sum_w/n/n));
    hist1.SetBinContent(ibin, sum_w/n);
  }
}

bool stop_vector(vector<unsigned>& vec, unsigned vec_len) {
  //immediately terminate for empty vectors
  //otherwise terminate when incr_vector reaches end
  if (vec.size() > vec_len) return false;
  return true;
}

void validate_weights() {

  //deep flavor
  map<int, vector<float>> btag_wps{
    {2016, vector<float>({0.0614, 0.3093, 0.7221})},
    {2017, vector<float>({0.0521, 0.3033, 0.7489})},
    {2018, vector<float>({0.0494, 0.2770, 0.7264})},
    {2022, vector<float>({0.0494, 0.2770, 0.7264})},
    {2023, vector<float>({0.0494, 0.2770, 0.7264})}
  };

  vector<EventWeighter> weighters;
  weighters.push_back(EventWeighter("2016APV", btag_wps[2016]));
  weighters.push_back(EventWeighter("2016", btag_wps[2016]));
  weighters.push_back(EventWeighter("2017", btag_wps[2017]));
  weighters.push_back(EventWeighter("2018", btag_wps[2018]));
  weighters.push_back(EventWeighter("2022", btag_wps[2022]));
  weighters.push_back(EventWeighter("2022EE", btag_wps[2022]));
  weighters.push_back(EventWeighter("2023", btag_wps[2023]));
  weighters.push_back(EventWeighter("2023BPix", btag_wps[2023]));

  vector<TriggerWeighter> trigger_weighters;
  trigger_weighters.push_back(TriggerWeighter("2016APV"));
  trigger_weighters.push_back(TriggerWeighter("2016"));
  trigger_weighters.push_back(TriggerWeighter("2017"));
  trigger_weighters.push_back(TriggerWeighter("2018"));
  trigger_weighters.push_back(TriggerWeighter("2022"));
  trigger_weighters.push_back(TriggerWeighter("2022EE"));
  trigger_weighters.push_back(TriggerWeighter("2023"));
  trigger_weighters.push_back(TriggerWeighter("2023BPix"));

  vector<int> years = {2016,2016,2017,2018,2022,2022,2023,2023};
  vector<string> years_string = {"2016APV","2016","2017","2018","2022","2022EE","2023","2023BPix"};

  vector<float> el_pt_bins = {7.0,15.0,20.0,35.0,50.0,100.0,500.0};
  vector<float> el_eta_bins = {-2.5,-2.0,-1.566,-1.444,-0.8,0.0,0.8,1.444,1.566,2.0,2.5};
  vector<float> el_phi_bins = {-3.1416,-1.2,-0.8};
  vector<float> mu_pt_bins = {5.0,6.0,7.0,8.0,10.0,12.0,15.0,20.0,25.0,30.0,35.0,40.0,50.0,60.0,80.0,120.0,200.0};
  vector<float> mu_eta_bins = {0.0,0.9,1.2,2.1,2.4};
  vector<float> jet_pt_bins = {30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,500.0,600.0,800.0,1000.0};
  vector<float> jet_eta_bins = {-2.4,-1.92,-1.44,-0.96,-0.48,0.0,0.48,0.96,1.44,1.92,2.4};
  vector<float> jet_flav_bins = {1,4,5};
  //vector<float> ph_pt_bins = {};
  //vector<float> ph_eta_bins = {};
  vector<float> photon_eta_bins = {-2.4,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.4};
  vector<float> photon_phi_bins = {-3.1416,-1.2,-0.8,3.1416};
  vector<float> photon_pt_bins = {15.0,20.0,35.0,50.0,100.0,200.0,500.0};
  vector<float> photon_r9_bins = {0.0,0.94,2.0};
  vector<bool> photon_ebee_bins = {true, false};
  vector<bool> photon_eveto_bins = {true, false};
  vector<vector<bool>> trig_decision_bins;
  vector<bool> truefalse = {true, false};
  for (bool tf1 : truefalse) { 
    for (bool tf2 : truefalse) {
      for (bool tf3 : truefalse) {
        for (bool tf4 : truefalse) {
          trig_decision_bins.push_back({tf1,tf2,tf3,tf4});
        }
      }
    }
  }
  vector<float> trig_el_pt_bins = {7.0,11.0,12.0,13.0,14.0,16.0,20.0,21.0,22.0,
                                   23.0,24.0,25.0,26.0,27.0,28.0,29.0,30.0,
                                   31.0,32.0,33.0,34.0,35.0,38.0,40.0,45.0,
                                   50.0,80.0,100.0,120.0,200.0,500.0};
  vector<float> trig_el_eta_bins = {-2.5,-2.0,-1.566,-1.4442,-0.8,0.0,0.8,
                                    1.4442,1.566,2.0,2.5};
  vector<float> trig_el_phi_bins = {-3.1416,-1.2,-0.8};
  vector<float> trig_mu_pt_bins = {5.0,7.75,8.0,8.1,8.25,8.5,10.0,15.0,16.75,
                                   17.0,17.1,17.25,18.0,20.0,23.0,23.75,24.0,
                                   24.25,24.5,25.0,26.0,26.75,27.0,27.25,27.5,
                                   29.0,30.0,32.0,40.0,60.0,100.0,120.0,200.0,
                                   500.0};
  vector<float> trig_mu_eta_bins = {0.0,0.9,1.2,2.1,2.4};

  pico_tree pico("","temp.root");

  bool check_electron_weights = true;
  bool check_muon_weights = true;
  bool check_photon_weights = true;
  bool check_trigger_weights = true;
  bool check_btag_weights = true;
  bool verbose = false;
  unsigned trig_nlep_max = 3;
  unsigned trig_nlep_min = 2;
  bool do_all_trig = true; //consider events not passing baseline
  bool auto_continue = true;

  for (unsigned iyear = 0; iyear < weighters.size(); iyear++) {

    vector<float> jet_df_bins;
    jet_df_bins.push_back(0.0);
    for (float wp : btag_wps[years[iyear]]) {
      jet_df_bins.push_back(wp);
    }
    jet_df_bins.push_back(1.0);

    cout << endl;
    cout << "---------------------------------------------------------------\n";
    cout << years_string[iyear] << endl;
    cout << "---------------------------------------------------------------\n";

    if (check_electron_weights) {
      cout << endl;
      cout << "Electron weights" << endl;
      for (unsigned ipt = 0; ipt < (el_pt_bins.size()-1); ipt++) {
        for (unsigned ieta = 0; ieta < (el_eta_bins.size()-1); ieta++) {
          for (bool el_sig : {true, false}) {
            pico.out_mc_id().clear();
            pico.out_mc_id().push_back(11);
            pico.out_mc_statusflag().clear();
            pico.out_mc_statusflag().push_back(0x2000);
            pico.out_mc_pt().clear();
            pico.out_mc_pt().push_back((el_pt_bins[ipt]+el_pt_bins[ipt+1])/2.0);
            pico.out_mc_eta().clear();
            pico.out_mc_eta().push_back((el_eta_bins[ieta]+el_eta_bins[ieta+1])/2.0);
            pico.out_mc_phi().clear();
            pico.out_mc_phi().push_back(0.0);
            pico.out_el_sig().clear();
            pico.out_el_sig().push_back(el_sig);
            pico.out_el_pt().clear();
            pico.out_el_pt().push_back((el_pt_bins[ipt]+el_pt_bins[ipt+1])/2.0);
            pico.out_el_eta().clear();
            pico.out_el_eta().push_back((el_eta_bins[ieta]+el_eta_bins[ieta+1])/2.0);
            pico.out_el_phi().clear();
            pico.out_el_phi().push_back(0.0);
            weighters[iyear].ElectronSF(pico);
            bool found_bad = sf_is_bad(pico.out_w_el()) || 
                             sf_is_bad(pico.out_sys_el()[0]) || 
                             sf_is_bad(pico.out_sys_el()[1]);
            if (verbose || found_bad) {
              cout << "pt: " << el_pt_bins[ipt] << "--" << el_pt_bins[ipt+1];
              cout << ", eta: " << el_eta_bins[ieta] << "--" << el_eta_bins[ieta+1];
              cout << ", sig: " << el_sig;
              cout << ", sf = " << pico.out_w_el() << ", up = " << pico.out_sys_el()[0];
              cout << ", dn = " << pico.out_sys_el()[1] << endl;
            }
            if (found_bad) {
              cout << "!!! Found bad SF" << endl;
              string temp;
              cin >> temp;
              if (temp != "c")
                return;
            }
          }
        }
      }
    }

    if (check_muon_weights) {
      cout << endl;
      cout << "Muon weights" << endl;
      for (unsigned ipt = 0; ipt < (mu_pt_bins.size()-1); ipt++) {
        for (unsigned ieta = 0; ieta < (mu_eta_bins.size()-1); ieta++) {
          for (bool mu_sig : {true, false}) {
            pico.out_mc_id().clear();
            pico.out_mc_id().push_back(13);
            pico.out_mc_statusflag().clear();
            pico.out_mc_statusflag().push_back(0x2000);
            pico.out_mc_pt().clear();
            pico.out_mc_pt().push_back((mu_pt_bins[ipt]+mu_pt_bins[ipt+1])/2.0);
            pico.out_mc_eta().clear();
            pico.out_mc_eta().push_back((mu_eta_bins[ieta]+mu_eta_bins[ieta+1])/2.0);
            pico.out_mc_phi().clear();
            pico.out_mc_phi().push_back(0.0);
            pico.out_mu_sig().clear();
            pico.out_mu_sig().push_back(mu_sig);
            pico.out_mu_pt().clear();
            pico.out_mu_pt().push_back((mu_pt_bins[ipt]+mu_pt_bins[ipt+1])/2.0);
            pico.out_mu_eta().clear();
            pico.out_mu_eta().push_back((mu_eta_bins[ieta]+mu_eta_bins[ieta+1])/2.0);
            pico.out_mu_phi().clear();
            pico.out_mu_phi().push_back(0.0);
            weighters[iyear].MuonSF(pico);
            bool found_bad = sf_is_bad(pico.out_w_mu()) || 
                             sf_is_bad(pico.out_sys_mu()[0]) || 
                             sf_is_bad(pico.out_sys_mu()[1]);
            if (verbose || found_bad) {
              cout << "pt: " << mu_pt_bins[ipt] << "--" << mu_pt_bins[ipt+1];
              cout << ", eta: " << mu_eta_bins[ieta] << "--" << mu_eta_bins[ieta+1];
              cout << ", sig: " << mu_sig;
              cout << ", sf = " << pico.out_w_mu() << ", up = " << pico.out_sys_mu()[0];
              cout << ", dn = " << pico.out_sys_mu()[1] << endl;
            }
            if (found_bad) {
              cout << "!!! Found bad SF" << endl;
              string temp;
              cin >> temp;
              if (temp != "c")
                return;
            }
          }
        }
      }
    }

    if (check_photon_weights) {
      cout << endl;
      cout << "Photon weights" << endl;
      for (unsigned ipt = 0; ipt < (photon_pt_bins.size()-1); ipt++) {
        for (unsigned ieta = 0; ieta < (photon_eta_bins.size()-1); ieta++) {
          for (unsigned iphi = 0; iphi < (photon_phi_bins.size()-1); iphi++) {
            for (unsigned ir9 = 0; ir9 < (photon_r9_bins.size()-1); ir9++) {
              for (bool eveto : photon_eveto_bins) {
                for (bool issig : {true, false}) {
                  float mean_eta = (photon_eta_bins[ieta]+photon_eta_bins[ieta+1])/2.0;
                  pico.out_photon_pt().clear();
                  pico.out_photon_pt().push_back((photon_pt_bins[ipt]+photon_pt_bins[ipt+1])/2.0);
                  pico.out_photon_phi().clear();
                  pico.out_photon_phi().push_back((photon_phi_bins[iphi]+photon_phi_bins[iphi+1])/2.0);
                  pico.out_photon_eta().clear();
                  pico.out_photon_eta().push_back(mean_eta);
                  pico.out_photon_idmva().clear();
                  pico.out_photon_idmva().push_back(0.9);
                  pico.out_photon_sig().clear();
                  pico.out_photon_sig().push_back(issig);
                  pico.out_photon_isScEtaEB().clear();
                  pico.out_photon_isScEtaEB().push_back(fabs(mean_eta)<1.5);
                  pico.out_photon_isScEtaEE().clear();
                  pico.out_photon_isScEtaEE().push_back(fabs(mean_eta)>1.5);
                  pico.out_photon_drmin().clear();
                  pico.out_photon_drmin().push_back(1.0);
                  pico.out_photon_elveto().clear();
                  pico.out_photon_elveto().push_back(eveto);
                  pico.out_photon_pflavor().clear();
                  pico.out_photon_pflavor().push_back(1);
                  pico.out_photon_r9().clear();
                  pico.out_photon_r9().push_back((photon_r9_bins[ir9]+photon_r9_bins[ir9+1])/2.0);
                  weighters[iyear].PhotonSF(pico);
                  bool found_bad = sf_is_bad(pico.out_w_photon()) ||
                                   sf_is_bad(pico.out_sys_photon()[0]) ||
                                   sf_is_bad(pico.out_sys_photon()[1]);
                  if (verbose || found_bad) {
                    cout << "  r9: " << photon_r9_bins[ir9] << "--" << photon_r9_bins[ir9+1];
                    cout << "  pt: " << photon_pt_bins[ipt] << "--" << photon_pt_bins[ipt+1];
                    cout << "  eta: " << photon_eta_bins[ieta] << "--" << photon_eta_bins[ieta+1];
                    cout << "  veto: " << eveto;
                    cout << ", sf = " << pico.out_w_photon() << ", up = " << pico.out_sys_photon()[0];
                    cout << ", dn = " << pico.out_sys_photon()[1] << endl;
                  }
                  if (found_bad) {
                    cout << "!!! Found bad SF" << endl;
                    string temp;
                    cin >> temp;
                    if (temp != "c")
                      return;
                  }
                }
              }
            }
          }
        }
      }
    }

    ////weighters[iyear].PhotonIDSF(pico,w_photon_id);
    //float w_photon_csev;
    //vector<float> sys_photon_csev;
    //weighters[iyear].PhotonCSEVSF(pico,w_photon_csev,sys_photon_csev);

    if (check_trigger_weights) {
      cout << endl;
      cout << "Trigger weights" << endl;
      //consider up to 4 leptons for now
      for (unsigned nel = 0; nel < trig_nlep_max; nel++) {
        for (unsigned nmu = max(trig_nlep_min-nel,0); nmu < trig_nlep_max-nel; nmu++) {
          for (vector<unsigned> iptel(nel,0); 
               stop_vector(iptel,nel);
               incr_vector(iptel,trig_el_pt_bins.size()-1)) {
            for (vector<unsigned> ietael(nel,0); 
                 stop_vector(ietael,nel);
                 incr_vector(ietael,trig_el_eta_bins.size()-1)) {
              for (vector<unsigned> iptmu(nmu,0); 
                   stop_vector(iptmu,nmu);
                   incr_vector(iptmu,trig_mu_pt_bins.size()-1)) {
                for (vector<unsigned> ietamu(nmu,0); 
                     stop_vector(ietamu,nmu);
                     incr_vector(ietamu,trig_mu_eta_bins.size()-1)) {
                  for (vector<bool>& trig_decision : trig_decision_bins) {
                    //assign variables for this bin
                    pico.out_trig_single_el() = trig_decision[0];
                    pico.out_trig_single_mu() = trig_decision[1];
                    pico.out_trig_double_el() = trig_decision[2];
                    pico.out_trig_double_mu() = trig_decision[3];
                    pico.out_el_sig().clear();
                    pico.out_el_pt().clear();
                    pico.out_el_eta().clear();
                    pico.out_mu_sig().clear();
                    pico.out_mu_pt().clear();
                    pico.out_mu_eta().clear();
                    float lead_el_pt = 0;
                    float subl_el_pt = 0;
                    float lead_mu_pt = 0;
                    float subl_mu_pt = 0;
                    for (unsigned iel = 0; iel < nel; iel++) {
                      pico.out_el_sig().push_back(true);
                      pico.out_el_pt().push_back(
                          (trig_el_pt_bins[iptel[iel]]+trig_el_pt_bins[iptel[iel]+1])/2.0);
                      pico.out_el_eta().push_back(
                          (trig_el_eta_bins[ietael[iel]]+trig_el_eta_bins[ietael[iel]+1])/2.0);
                      if (pico.out_el_pt().back() > lead_el_pt) {
                        subl_el_pt = lead_el_pt;
                        lead_el_pt = pico.out_el_pt().back();
                      }
                      else if (pico.out_el_pt().back() > subl_el_pt) {
                        subl_el_pt = pico.out_el_pt().back();
                      }
                    }
                    for (unsigned imu = 0; imu < nmu; imu++) {
                      pico.out_mu_sig().push_back(true);
                      pico.out_mu_pt().push_back(
                          (trig_mu_pt_bins[iptmu[imu]]+trig_mu_pt_bins[iptmu[imu]+1])/2.0);
                      pico.out_mu_eta().push_back(
                          (trig_mu_eta_bins[ietamu[imu]]+trig_mu_eta_bins[ietamu[imu]+1])/2.0);
                      if (pico.out_mu_pt().back() > lead_mu_pt) {
                        subl_mu_pt = lead_mu_pt;
                        lead_mu_pt = pico.out_mu_pt().back();
                      }
                      else if (pico.out_mu_pt().back() > subl_mu_pt) {
                        subl_mu_pt = pico.out_mu_pt().back();
                      }
                    }
                    //get trigger weight and check sanity
                    trigger_weighters[iyear].GetSF(pico);
                    vector<float> sfs = {pico.out_w_trig(), 
                        pico.out_sys_trig_el()[0], pico.out_sys_trig_el()[1], 
                        pico.out_sys_trig_mu()[0], pico.out_sys_trig_mu()[1]};
                    bool found_bad = (sf_is_bad(sfs[0]) || sf_is_bad(sfs[1]) 
                        || sf_is_bad(sfs[2]) || sf_is_bad(sfs[3]) 
                        || sf_is_bad(sfs[4]));
                    bool trig_in_sr = false;
                    if (pico.out_trig_single_el() && lead_el_pt > 30) trig_in_sr = true;
                    if (pico.out_trig_double_el() && lead_el_pt > 25 && subl_el_pt > 15) trig_in_sr = true;
                    if (pico.out_trig_single_mu() && lead_mu_pt > 25) trig_in_sr = true;
                    if (pico.out_trig_double_mu() && lead_mu_pt > 20 && subl_el_pt > 10) trig_in_sr = true;
                    if (verbose || (found_bad && (trig_in_sr || do_all_trig))) {
                      cout << "el_pt: ";
                      for (unsigned iel = 0; iel<nel; iel++)
                        cout << pico.out_el_pt()[iel] << ", ";
                      cout << endl;
                      cout << "el_eta: ";
                      for (unsigned iel = 0; iel<nel; iel++)
                        cout << pico.out_el_eta()[iel] << ", ";
                      cout << endl;
                      cout << "mu_pt: ";
                      for (unsigned imu = 0; imu<nmu; imu++)
                        cout << pico.out_mu_pt()[imu] << ", ";
                      cout << endl;
                      cout << "mu_eta: ";
                      for (unsigned imu = 0; imu<nmu; imu++)
                        cout << pico.out_mu_eta()[imu] << ", ";
                      cout << endl;
                      cout << "trig_single_el: " << trig_decision[0] << endl;
                      cout << "trig_single_mu: " << trig_decision[1] << endl;
                      cout << "trig_double_el: " << trig_decision[2] << endl;
                      cout << "trig_double_mu: " << trig_decision[3] << endl;
                      cout << "sf = " << sfs[0] << ", elup = " << sfs[1] 
                           << ", eldn = " << sfs[2] << ", muup = " << sfs[3] 
                           << ", mudn = " << sfs[4] << endl;
                    }
                    if (found_bad && (trig_in_sr || do_all_trig) && !auto_continue) {
                      cout << "!!! Found bad SF" << endl;
                      string temp;
                      cin >> temp;
                      if (temp != "c")
                        return;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    //PU doesn't really need checking since straight from file

    if (check_btag_weights) {
      cout << endl;
      cout << "b-tag weights" << endl;
      for (unsigned ipt = 0; ipt < (jet_pt_bins.size()-1); ipt++) {
        for (unsigned ieta = 0; ieta < (jet_eta_bins.size()-1); ieta++) {
          for (unsigned iflav = 0; iflav < jet_flav_bins.size(); iflav++) {
            for (unsigned idf = 0; idf < (jet_df_bins.size()-1); idf++) {
              pico.out_jet_isgood().clear();
              pico.out_jet_isgood().push_back(true);
              pico.out_jet_hflavor().clear();
              pico.out_jet_hflavor().push_back(jet_flav_bins[iflav]);
              pico.out_jet_deepflav().clear();
              pico.out_jet_deepflav().push_back(jet_flav_bins[iflav]);
              pico.out_jet_pt().clear();
              pico.out_jet_pt().push_back((jet_pt_bins[ipt]+jet_pt_bins[ipt+1])/2.0);
              pico.out_jet_eta().clear();
              pico.out_jet_eta().push_back((jet_eta_bins[ieta]+jet_eta_bins[ieta+1])/2.0);
              weighters[iyear].bTaggingSF(pico);
              bool found_bad = sf_is_bad(pico.out_w_bhig_df()) || 
                               sf_is_bad(pico.out_sys_bchig()[0]) || 
                               sf_is_bad(pico.out_sys_bchig()[1]) || 
                               sf_is_bad(pico.out_sys_udsghig()[0]) || 
                               sf_is_bad(pico.out_sys_udsghig()[1]);
              if (verbose || found_bad) {
                cout << "pt: " << jet_pt_bins[ipt] << "--" << jet_pt_bins[ipt+1];
                cout << ", eta: " << jet_eta_bins[ieta] << "--" << jet_eta_bins[ieta+1];
                cout << ", flavor: " << jet_flav_bins[iflav];
                cout << ", df: " << jet_df_bins[idf] << "--" << jet_df_bins[idf+1];
                cout << ", sf = " << pico.out_w_bhig_df();
                cout << ", up = " << pico.out_sys_bchig()[0];
                cout << ", dn = " << pico.out_sys_bchig()[1];
                cout << ", up = " << pico.out_sys_udsghig()[0];
                cout << ", dn = " << pico.out_sys_udsghig()[1] << endl;
              }
              if (found_bad) {
                cout << "!!! Found bad SF" << endl;
                return;
              }
            }
          }
        }
      }
    }

  } // loop over years
}

//makes average SF plots for validation
void make_sf_plots(const char* filename, int year_idx, string output_prefix) {
  //deep flavor
  map<int, vector<float>> btag_wps{
    {2016, vector<float>({0.0614, 0.3093, 0.7221})},
    {2017, vector<float>({0.0521, 0.3033, 0.7489})},
    {2018, vector<float>({0.0494, 0.2770, 0.7264})},
    {2022, vector<float>({0.0494, 0.2770, 0.7264})},
    {2023, vector<float>({0.0494, 0.2770, 0.7264})}
  };

  vector<EventWeighter> weighters;
  weighters.push_back(EventWeighter("2016APV", btag_wps[2016]));
  weighters.push_back(EventWeighter("2016", btag_wps[2016]));
  weighters.push_back(EventWeighter("2017", btag_wps[2017]));
  weighters.push_back(EventWeighter("2018", btag_wps[2018]));
  weighters.push_back(EventWeighter("2022", btag_wps[2022]));
  weighters.push_back(EventWeighter("2022EE", btag_wps[2022]));
  weighters.push_back(EventWeighter("2023", btag_wps[2023]));
  weighters.push_back(EventWeighter("2023BPix", btag_wps[2023]));

  vector<TriggerWeighter> trigger_weighters;
  trigger_weighters.push_back(TriggerWeighter("2016APV"));
  trigger_weighters.push_back(TriggerWeighter("2016"));
  trigger_weighters.push_back(TriggerWeighter("2017"));
  trigger_weighters.push_back(TriggerWeighter("2018"));
  trigger_weighters.push_back(TriggerWeighter("2022"));
  trigger_weighters.push_back(TriggerWeighter("2022EE"));
  trigger_weighters.push_back(TriggerWeighter("2023"));
  trigger_weighters.push_back(TriggerWeighter("2023BPix"));

  vector<int> years = {2016,2016,2017,2018,2022,2022,2023,2023};
  vector<string> years_string = {"2016APV","2016","2017","2018","2022","2022EE","2023","2023BPix"};
  //bins
  vector<double> el_pt_bins = {7.0,15.0,20.0,35.0,50.0,100.0,500.0};
  vector<float> el_eta_bins = {-2.5,-2.0,-1.566,-1.444,-0.8,0.0,0.8,1.444,1.566,2.0,2.5};
  vector<double> mu_pt_bins = {5.0,6.0,7.0,8.0,10.0,12.0,15.0,20.0,25.0,30.0,35.0,40.0,50.0,60.0,80.0,120.0,200.0};
  vector<float> mu_eta_bins = {0.0,0.9,1.2,2.1,2.4};
  vector<double> ph_pt_bins = {15.0,20.0,35.0,50.0,100.0,200.0,500.0};
  vector<float> ph_eta_bins = {-2.4,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.4};

  vector<TH1D> el_sf_hists;
  vector<TH1D> mu_sf_hists;
  vector<TH1D> ph_sf_hists;
  vector<TH1D> el_norm_hists;
  vector<TH1D> mu_norm_hists;
  vector<TH1D> ph_norm_hists;
  for (unsigned ieta = 0; ieta < el_eta_bins.size()-1; ieta++) {
    el_sf_hists.push_back(TH1D(("el_sf_hist"+to_string(ieta)).c_str(), 
        (to_string(el_eta_bins[ieta])+"<#eta<"+to_string(el_eta_bins[ieta+1])
        +"; p_{T} [GeV]; Average SF").c_str(), el_pt_bins.size()-1, 
        &el_pt_bins.at(0)));
    el_norm_hists.push_back(TH1D(("el_norm_hist"+to_string(ieta)).c_str(), "",
        el_pt_bins.size()-1, &el_pt_bins.at(0)));
  }
  for (unsigned ieta = 0; ieta < mu_eta_bins.size()-1; ieta++) {
    mu_sf_hists.push_back(TH1D(("mu_sf_hist"+to_string(ieta)).c_str(), 
        (to_string(mu_eta_bins[ieta])+"<#eta<"+to_string(mu_eta_bins[ieta+1])
        +"; p_{T} [GeV]; Average SF").c_str(), mu_pt_bins.size()-1, 
        &mu_pt_bins.at(0)));
    mu_norm_hists.push_back(TH1D(("mu_norm_hist"+to_string(ieta)).c_str(), "",
        mu_pt_bins.size()-1, &mu_pt_bins.at(0)));
  }
  for (unsigned ieta = 0; ieta < ph_eta_bins.size()-1; ieta++) {
    ph_sf_hists.push_back(TH1D(("ph_sf_hist"+to_string(ieta)).c_str(), 
        (to_string(ph_eta_bins[ieta])+"<#eta<"+to_string(ph_eta_bins[ieta+1])
        +"; p_{T} [GeV]; Average SF").c_str(), ph_pt_bins.size()-1, 
        &ph_pt_bins.at(0)));
    ph_norm_hists.push_back(TH1D(("ph_norm_hist"+to_string(ieta)).c_str(), "",
        ph_pt_bins.size()-1, &ph_pt_bins.at(0)));
  }

  vector<int>   *mc_id(nullptr);
  vector<int>   *mc_statusflag(nullptr);
  vector<float> *mc_pt(nullptr);
  vector<float> *mc_eta(nullptr);
  vector<float> *mc_phi(nullptr);

  vector<bool>  *el_sig(nullptr);
  vector<float> *el_pt(nullptr);
  vector<float> *el_eta(nullptr);
  vector<float> *el_phi(nullptr);

  vector<bool>  *mu_sig(nullptr);
  vector<float> *mu_pt(nullptr);
  vector<float> *mu_eta(nullptr);
  vector<float> *mu_phi(nullptr);

  vector<bool>  *photon_sig(nullptr);
  vector<float> *photon_pt(nullptr);
  vector<float> *photon_eta(nullptr);
  vector<float> *photon_phi(nullptr);
  vector<float> *photon_idmva(nullptr);
  vector<bool>  *photon_isScEtaEB(nullptr);
  vector<bool>  *photon_isScEtaEE(nullptr);
  vector<float> *photon_drmin(nullptr);
  vector<bool>  *photon_elveto(nullptr);
  vector<int>   *photon_pflavor(nullptr);
  vector<float> *photon_r9(nullptr);

  bool trig_single_el;
  bool trig_single_mu;
  bool trig_double_el;
  bool trig_double_mu;

  int nel;
  int nmu;
  int nphoton;

  pico_tree pico("","temp.root");

  TChain tree("tree");
  tree.Add(filename);

  tree.SetBranchAddress("mc_id", &mc_id);
  tree.SetBranchAddress("mc_statusflag", &mc_statusflag);
  tree.SetBranchAddress("mc_pt", &mc_pt);
  tree.SetBranchAddress("mc_eta", &mc_eta);
  tree.SetBranchAddress("mc_phi", &mc_phi);

  tree.SetBranchAddress("el_sig", &el_sig);
  tree.SetBranchAddress("el_pt", &el_pt);
  tree.SetBranchAddress("el_eta", &el_eta);
  tree.SetBranchAddress("el_phi", &el_phi);

  tree.SetBranchAddress("mu_sig", &mu_sig);
  tree.SetBranchAddress("mu_pt", &mu_pt);
  tree.SetBranchAddress("mu_eta", &mu_eta);
  tree.SetBranchAddress("mu_phi", &mu_phi);

  tree.SetBranchAddress("photon_sig", &photon_sig);
  tree.SetBranchAddress("photon_pt", &photon_pt);
  tree.SetBranchAddress("photon_eta", &photon_eta);
  tree.SetBranchAddress("photon_phi", &photon_phi);
  tree.SetBranchAddress("photon_idmva", &photon_idmva);
  tree.SetBranchAddress("photon_isScEtaEB", &photon_isScEtaEB);
  tree.SetBranchAddress("photon_isScEtaEE", &photon_isScEtaEE);
  tree.SetBranchAddress("photon_drmin", &photon_drmin);
  tree.SetBranchAddress("photon_elveto", &photon_elveto);
  tree.SetBranchAddress("photon_pflavor", &photon_pflavor);
  tree.SetBranchAddress("photon_r9", &photon_r9);

  tree.SetBranchAddress("trig_single_el", &trig_single_el);
  tree.SetBranchAddress("trig_single_mu", &trig_single_mu);
  tree.SetBranchAddress("trig_double_el", &trig_double_el);
  tree.SetBranchAddress("trig_double_mu", &trig_double_mu);
  tree.SetBranchAddress("nel", &nel);
  tree.SetBranchAddress("nmu", &nmu);
  tree.SetBranchAddress("nphoton", &nphoton);

  //event loop
  for (long ievent = 0; ievent < tree.GetEntries(); ievent++) {
    tree.GetEvent(ievent);

    pico.out_mc_id() = *mc_id;
    pico.out_mc_statusflag() = *mc_statusflag;
    pico.out_mc_pt() = *mc_pt;
    pico.out_mc_eta() = *mc_eta;
    pico.out_mc_phi() = *mc_phi;

    pico.out_el_sig() = *el_sig;
    pico.out_el_pt() = *el_pt;
    pico.out_el_eta() = *el_eta;
    pico.out_el_phi() = *el_phi;

    pico.out_mu_sig() = *mu_sig;
    pico.out_mu_pt() = *mu_pt;
    pico.out_mu_eta() = *mu_eta;
    pico.out_mu_phi() = *mu_phi;

    pico.out_photon_sig() = *photon_sig;
    pico.out_photon_pt() = *photon_pt;
    pico.out_photon_eta() = *photon_eta;
    pico.out_photon_phi() = *photon_phi;
    pico.out_photon_idmva() = *photon_idmva;
    pico.out_photon_isScEtaEB() = *photon_isScEtaEB;
    pico.out_photon_isScEtaEE() = *photon_isScEtaEE;
    pico.out_photon_drmin() = *photon_drmin;
    pico.out_photon_elveto() = *photon_elveto;
    pico.out_photon_pflavor() = *photon_pflavor;
    pico.out_photon_r9() = *photon_r9;

    pico.out_trig_single_el() = trig_single_el;
    pico.out_trig_single_mu() = trig_single_mu;
    pico.out_trig_double_el() = trig_double_el;
    pico.out_trig_double_mu() = trig_double_mu;

    weighters[year_idx].ElectronSF(pico);
    weighters[year_idx].MuonSF(pico);
    weighters[year_idx].PhotonSF(pico);
    trigger_weighters[year_idx].GetSF(pico);

    if (nel == 2 && nmu == 0) {
      for (unsigned iel = 0; iel < el_sig->size(); iel++) {
        if (el_sig->at(iel)) {
          float eta = el_eta->at(iel);
          float pt = el_pt->at(iel);
          for (unsigned ieta = 0; ieta < el_eta_bins.size()-1; ieta++) {
            if (eta >= el_eta_bins[ieta] && eta < el_eta_bins[ieta+1]) {
              float avr_weight = pow(pico.out_w_el(), 0.5)
                                 * pow(pico.out_w_trig(), 0.5);
              el_sf_hists[ieta].Fill(pt, avr_weight);
              el_norm_hists[ieta].Fill(pt, 1.0);
            }
          }
        }
      }
    }

    if (nel == 0 && nmu == 2) {
      for (unsigned imu = 0; imu < mu_sig->size(); imu++) {
        if (mu_sig->at(imu)) {
          float eta = mu_eta->at(imu);
          float pt = mu_pt->at(imu);
          for (unsigned ieta = 0; ieta < mu_eta_bins.size()-1; ieta++) {
            if (eta >= mu_eta_bins[ieta] && eta < mu_eta_bins[ieta+1]) {
              float avr_weight = pow(pico.out_w_mu(), 0.5)
                                 * pow(pico.out_w_trig(), 0.5);
              mu_sf_hists[ieta].Fill(pt, avr_weight);
              mu_norm_hists[ieta].Fill(pt, 1.0);
            }
          }
        }
      }
    }

    if (nphoton == 1) {
      for (unsigned iph = 0; iph < photon_sig->size(); iph++) {
        if (photon_sig->at(iph)) {
          float eta = photon_eta->at(iph);
          float pt = photon_pt->at(iph);
          for (unsigned ieta = 0; ieta < ph_eta_bins.size()-1; ieta++) {
            if (eta >= ph_eta_bins[ieta] && eta < ph_eta_bins[ieta+1]) {
              ph_sf_hists[ieta].Fill(pt, pico.out_w_photon());
              ph_norm_hists[ieta].Fill(pt, 1.0);
            }
          }
        }
      }
    }

  }

  TCanvas c;
  for (unsigned ieta = 0; ieta < el_eta_bins.size()-1; ieta++) {
    average_hist(el_sf_hists[ieta], el_norm_hists[ieta]);
    el_sf_hists[ieta].Draw();
    c.SaveAs((output_prefix + "_elhist"+to_string(ieta)+".root").c_str());
  }
  for (unsigned ieta = 0; ieta < mu_eta_bins.size()-1; ieta++) {
    average_hist(mu_sf_hists[ieta], mu_norm_hists[ieta]);
    mu_sf_hists[ieta].Draw();
    c.SaveAs((output_prefix + "_muhist"+to_string(ieta)+".root").c_str());
  }
  for (unsigned ieta = 0; ieta < ph_eta_bins.size()-1; ieta++) {
    average_hist(ph_sf_hists[ieta], ph_norm_hists[ieta]);
    ph_sf_hists[ieta].Draw();
    c.SaveAs((output_prefix + "_phhist"+to_string(ieta)+".root").c_str());
  }

}

int main() {
  //validate_weights();
  make_sf_plots("/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_redwood_v0_systsignal/2016APV/mc/unskimmed/pico_GluGluHToZG_ZToLL_M-125_TuneCP5*", 0, "sfs2016APV");
  make_sf_plots("/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_redwood_v0_systsignal/2016/mc/unskimmed/pico_GluGluHToZG_ZToLL_M-125_TuneCP5*", 1, "sfs2016");
  make_sf_plots("/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_redwood_v0_systsignal/2017/mc/unskimmed/pico_GluGluHToZG_ZToLL_M-125_TuneCP5*", 2, "sfs2017");
  make_sf_plots("/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_redwood_v0_systsignal/2018/mc/unskimmed/pico_GluGluHToZG_ZToLL_M-125_TuneCP5*", 3, "sfs2018");
  make_sf_plots("/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0_systsignal/2022/mc/unskimmed/pico_GluGluHtoZG_Zto2L_M-125_TuneCP5*", 4, "sfs2022");
  make_sf_plots("/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0_systsignal/2022EE/mc/unskimmed/pico_GluGluHtoZG_Zto2L_M-125_TuneCP5*", 5, "sfs2022EE");
  make_sf_plots("/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0_systsignal/2023/mc/unskimmed/pico_GluGluHtoZG_Zto2L_M-125_TuneCP5*", 6, "sfs2023");
  make_sf_plots("/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0_systsignal/2023BPix/mc/unskimmed/pico_GluGluHtoZG_Zto2L_M-125_TuneCP5*", 7, "sfs2023BPix");
  return 0;
}
