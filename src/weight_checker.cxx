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

using std::cout;
using std::endl;
using std::isnan;
using std::isinf;
using std::map;
using std::string;
using std::vector;

bool sf_is_bad(float value) {
  return isnan(value) || isinf(value) || (fabs(value)>20.0);
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
  weighters.push_back(EventWeighter(2016, true, btag_wps[2016]));
  weighters.push_back(EventWeighter(2016, false, btag_wps[2016]));
  weighters.push_back(EventWeighter(2017, false, btag_wps[2017]));
  weighters.push_back(EventWeighter(2018, false, btag_wps[2018]));

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

  pico_tree pico("","temp.root");

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
          cout << "pt: " << el_pt_bins[ipt] << "--" << el_pt_bins[ipt+1];
          cout << ", eta: " << el_eta_bins[ieta] << "--" << el_eta_bins[ieta+1];
          cout << ", sig: " << el_sig;
          cout << ", sf = " << pico.out_w_el() << ", up = " << pico.out_sys_el()[0];
          cout << ", dn = " << pico.out_sys_el()[1] << endl;
          if (sf_is_bad(pico.out_w_el()) || 
              sf_is_bad(pico.out_sys_el()[0]) || 
              sf_is_bad(pico.out_sys_el()[1]))
            return 1;
        }
      }
    }

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
          cout << "pt: " << mu_pt_bins[ipt] << "--" << mu_pt_bins[ipt+1];
          cout << ", eta: " << mu_eta_bins[ieta] << "--" << mu_eta_bins[ieta+1];
          cout << ", sig: " << mu_sig;
          cout << ", sf = " << pico.out_w_mu() << ", up = " << pico.out_sys_mu()[0];
          cout << ", dn = " << pico.out_sys_mu()[1] << endl;
          if (sf_is_bad(pico.out_w_mu()) || 
              sf_is_bad(pico.out_sys_mu()[0]) || 
              sf_is_bad(pico.out_sys_mu()[1]))
            return 1;
        }
      }
    }

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
          cout << "r9: " << photon_r9_bins[ir9] << "--" << photon_r9_bins[ir9+1];
          cout << ", EB: " << photon_iseb << ", veto: " << eveto;
          cout << ", sf = " << w_photon_csev << ", up = " << sys_photon_csev[0];
          cout << ", dn = " << sys_photon_csev[1] << endl;
        }
      }
    }

    ////weighters[iyear].PhotonIDSF(pico,w_photon_id);
    //float w_photon_csev;
    //vector<float> sys_photon_csev;
    //weighters[iyear].PhotonCSEVSF(pico,w_photon_csev,sys_photon_csev);

    //skip for now
    //cout << endl;
    //cout << "Trigger weights" << endl;
    //trigger_weighter[iyear].GetSF(pico);

    //PU doesn't really need checking
    //cout << endl;
    //cout << "PU weights" << endl;
    //weighters[iyear].PileupSF(pico);

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
            cout << "pt: " << jet_pt_bins[ipt] << "--" << jet_pt_bins[ipt+1];
            cout << ", eta: " << jet_eta_bins[ieta] << "--" << jet_eta_bins[ieta+1];
            cout << ", flavor: " << jet_flav_bins[iflav];
            cout << ", df: " << jet_df_bins[idf] << "--" << jet_df_bins[idf+1];
            cout << ", sf = " << pico.out_w_bhig_df();
            cout << ", up = " << pico.out_sys_bchig()[0];
            cout << ", dn = " << pico.out_sys_bchig()[1];
            cout << ", up = " << pico.out_sys_udsghig()[0];
            cout << ", dn = " << pico.out_sys_udsghig()[1] << endl;
            if (sf_is_bad(pico.out_w_bhig_df()) || 
                sf_is_bad(pico.out_sys_bchig()[0]) || 
                sf_is_bad(pico.out_sys_bchig()[1]) || 
                sf_is_bad(pico.out_sys_udsghig()[0]) || 
                sf_is_bad(pico.out_sys_udsghig()[1]))
              return 1;
          }
        }
      }
    }


  } // loop over years

  return 0;
}
