/**
 * Script to validate weights
 */

#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "pico_tree.hpp"
#include "event_weighter.hpp"
#include "trigger_weighter.hpp"

using std::cin;
using std::cout;
using std::endl;
using std::isnan;
using std::isinf;
using std::map;
using std::string;
using std::vector;

unsigned max(unsigned a, unsigned b) {
  if (a>b) return a;
  return b;
}

bool sf_is_bad(float value) {
  return isnan(value) || isinf(value) || (fabs(value)>25.0);
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

bool stop_vector(vector<unsigned>& vec, unsigned vec_len) {
  //immediately terminate for empty vectors
  //otherwise terminate when incr_vector reaches end
  if (vec.size() > vec_len) return false;
  return true;
}

int main() {

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

  vector<TriggerWeighter> trigger_weighters;
  trigger_weighters.push_back(TriggerWeighter(2016, true));
  trigger_weighters.push_back(TriggerWeighter(2016, false));
  trigger_weighters.push_back(TriggerWeighter(2017, false));
  trigger_weighters.push_back(TriggerWeighter(2018, false));

  vector<int> years = {2016,2016,2017,2018};
  vector<string> years_string = {"2016APV","2016","2017","2018"};

  vector<float> el_pt_bins = {7.0,10.0,20.0,35.0,50.0,100.0,200.0,500.0};
  vector<float> el_eta_bins = {-2.5,-2.0,-1.566,-1.444,-0.8,0.0,0.8,1.444,1.566,2.0,2.5};
  vector<float> mu_pt_bins = {5.0,6.0,8.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,60.0,120.0,200.0};
  vector<float> mu_eta_bins = {0.0,0.9,1.2,2.1,2.4};
  vector<float> jet_pt_bins = {30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,500.0,600.0,800.0,1000.0};
  vector<float> jet_eta_bins = {-2.4,-1.92,-1.44,-0.96,-0.48,0.0,0.48,0.96,1.44,1.92,2.4};
  vector<float> jet_flav_bins = {1,4,5};
  //vector<float> ph_pt_bins = {};
  //vector<float> ph_eta_bins = {};
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
  vector<float> trig_el_eta_bins = {0.0,0.8,1.4442,1.566,2.5};
  vector<float> trig_mu_pt_bins = {5.0,7.75,8.0,8.1,8.25,8.5,10.0,15.0,16.75,
                                   17.0,17.1,17.25,18.0,20.0,23.0,23.75,24.0,
                                   24.25,24.5,25.0,26.0,26.75,27.0,27.25,27.5,
                                   29.0,30.0,32.0,40.0,60.0,100.0,120.0,200.0,
                                   500.0};
  vector<float> trig_mu_eta_bins = {0.0,0.9,1.2,2.1,2.4};

  pico_tree pico("","temp.root");

  bool check_electron_weights = true;
  bool check_muon_weights = true;
  bool check_photon_csev_weights = true;
  bool check_trigger_weights = true;
  bool check_btag_weights = true;
  bool verbose = false;
  unsigned trig_nlep_max = 3;
  unsigned trig_nlep_min = 2;
  bool do_all_trig = false; //consider events not passing baseline

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
                return 1;
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
                return 1;
            }
          }
        }
      }
    }

    if (check_photon_csev_weights) {
      cout << endl;
      cout << "Photon CSEV weights" << endl;
      for (unsigned ir9 = 0; ir9 < (photon_r9_bins.size()-1); ir9++) {
        for (bool photon_iseb : photon_ebee_bins) {
          for (bool eveto : photon_eveto_bins) {
            pico.out_photon_pt().clear();
            pico.out_photon_pt().push_back(20.0);
            pico.out_photon_idmva().clear();
            pico.out_photon_idmva().push_back(0.9);
            pico.out_photon_isScEtaEB().clear();
            pico.out_photon_isScEtaEB().push_back(photon_iseb);
            pico.out_photon_isScEtaEE().clear();
            pico.out_photon_isScEtaEE().push_back(!photon_iseb);
            pico.out_photon_drmin().clear();
            pico.out_photon_drmin().push_back(1.0);
            pico.out_photon_elveto().clear();
            pico.out_photon_elveto().push_back(eveto);
            pico.out_photon_pflavor().clear();
            pico.out_photon_pflavor().push_back(1);
            pico.out_photon_r9().clear();
            pico.out_photon_r9().push_back((photon_r9_bins[ir9]+photon_r9_bins[ir9+1])/2.0);
            float w_photon_csev;
            vector<float> sys_photon_csev;
            sys_photon_csev.resize(2,1.0);
            weighters[iyear].PhotonCSEVSF(pico, w_photon_csev, sys_photon_csev);
            bool found_bad = sf_is_bad(w_photon_csev) ||
                             sf_is_bad(sys_photon_csev[0]) ||
                             sf_is_bad(sys_photon_csev[1]);
            if (verbose || found_bad) {
              cout << "r9: " << photon_r9_bins[ir9] << "--" << photon_r9_bins[ir9+1];
              cout << ", EB: " << photon_iseb << ", veto: " << eveto;
              cout << ", sf = " << w_photon_csev << ", up = " << sys_photon_csev[0];
              cout << ", dn = " << sys_photon_csev[1] << endl;
            }
            if (found_bad) {
              cout << "!!! Found bad SF" << endl;
              string temp;
              cin >> temp;
              if (temp != "c")
                return 1;
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
                    vector<float> sfs = trigger_weighters[iyear].GetSF(pico);
                    bool found_bad = sf_is_bad(sfs[0]) || sf_is_bad(sfs[1]) || sf_is_bad(sfs[2]);
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
                      cout << "sf = " << sfs[0] << ", up = " << sfs[1] << ", dn = " << sfs[2] << endl;
                    }
                    if (found_bad && (trig_in_sr || do_all_trig)) {
                      cout << "!!! Found bad SF" << endl;
                      string temp;
                      cin >> temp;
                      if (temp != "c")
                        return 1;
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
                return 1;
              }
            }
          }
        }
      }
    }

  } // loop over years

  return 0;
}
