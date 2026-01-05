#include <iostream>
#include <ctime>
#include <getopt.h>

#include "pico_tree.hpp"
#include "corrections_tree.hpp"
#include "utilities.hpp"
#include "cross_sections.hpp"

#include "TError.h"

using namespace std;

namespace {
  string in_file = "";
  string in_dir = "";
  string corr_file = "";
}

void GetOptions(int argc, char *argv[]);

vector<float> Caster(vector<double> dvec){
  size_t vecsize = dvec.size();
  vector<float> fvec;
  for(unsigned int i(0);i<vecsize;i++){
    fvec.push_back(static_cast<float>(dvec[i]));
  }
  return fvec;
}

int main(int argc, char *argv[]){
  // gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches       
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  string in_file_path = in_dir+"/"+in_file;
  cout<<"Running on input file: "<<in_file_path<<endl;

  string out_file = CopyReplaceAll(in_dir, "/raw_pico/","/unskimmed/") + CopyReplaceAll(in_file, "raw_","");
  corr_file = CopyReplaceAll(in_dir, "/raw_pico/","/corrections/")+corr_file;

  bool is_zgamma = Contains(in_dir, "zgamma");

  cout<<"Corrections file: "<<corr_file<<endl;
  corrections_tree corr(corr_file);
  if (corr.GetEntries()!=1) {
    cout<<"ERROR:: Corrections file has "<<corr.GetEntries()<<", expected 1."<<endl;
    exit(1);
  }
  corr.GetEntry(0);

  pico_tree pico(in_file_path, out_file);
  for(long entry(0); entry<pico.GetEntries(); entry++){

    pico.GetEntry(entry);
    if (entry%100000==0) {
      cout<<"Processing event: "<<entry<<endl;
    }

    pico.out_sys_lep().resize(2); pico.out_sys_fs_lep().resize(2);
    if(pico.nlep()==0) { // load from calculated correction
      pico.out_w_lep()         = static_cast<float>(corr.w_lep());
      pico.out_w_fs_lep()      = static_cast<float>(corr.w_fs_lep());
      pico.out_sys_lep()       = Caster(corr.sys_lep());
      pico.out_sys_fs_lep()    = Caster(corr.sys_fs_lep());
    } else { //load from original tree
      pico.out_w_lep()         = pico.w_lep();
      pico.out_w_fs_lep()      = pico.w_fs_lep();
      pico.out_sys_lep()       = pico.sys_lep();
      pico.out_sys_fs_lep()    = pico.sys_fs_lep();
    }
    pico.out_w_phshape()   = pico.w_phshape()*static_cast<float>(corr.w_phshape());
    
    pico.out_w_trig()     = pico.w_trig();
    pico.out_w_pu()       = pico.w_pu()*static_cast<float>(corr.w_pu());
    pico.out_w_nnlo()     = pico.w_nnlo()*static_cast<float>(corr.w_nnlo());

    if (is_zgamma) {
      if (((pico.type() >= 6000 && pico.type() < 7000) ||
           (pico.type() >= 17000 && pico.type() < 18000))
          && pico.nllphoton() >= 1) {
        pico.out_w_isr() = pico.w_isr()*static_cast<float>(corr.w_isr());
      }
      else {
        pico.out_w_isr() = pico.w_isr();
      }
    }
    else {
      pico.out_w_isr()      = pico.w_isr()*static_cast<float>(corr.w_isr());
    }

    float btag_weight = pico.w_bhig();
    if (is_zgamma) {
      pico.out_w_lep() = pico.w_el()*pico.w_mu();
      pico.out_w_fs_lep() = 1.0;
      pico.out_sys_fs_lep()[0] = 1.0;
      pico.out_sys_fs_lep()[1] = 1.0;
      btag_weight = pico.w_btag_df();
    }

    pico.out_sys_trig().resize(2);
    pico.out_sys_trig_el().resize(2);
    pico.out_sys_trig_mu().resize(2);
    if (pico.nel()>0) {
      if (pico.trig_single_el() || pico.trig_double_el()) {
        pico.out_w_trig() = pico.w_trig()*static_cast<float>(corr.w_zvtx_pass());
        pico.out_sys_trig()[0] = pico.sys_trig()[0]*static_cast<float>(corr.w_zvtx_pass());
        pico.out_sys_trig()[1] = pico.sys_trig()[1]*static_cast<float>(corr.w_zvtx_pass());
        pico.out_sys_trig_el()[0] = pico.sys_trig_el()[0]*static_cast<float>(corr.w_zvtx_pass());
        pico.out_sys_trig_el()[1] = pico.sys_trig_el()[1]*static_cast<float>(corr.w_zvtx_pass());
        pico.out_sys_trig_mu()[0] = pico.sys_trig_mu()[0]*static_cast<float>(corr.w_zvtx_pass());
        pico.out_sys_trig_mu()[1] = pico.sys_trig_mu()[1]*static_cast<float>(corr.w_zvtx_pass());
      }
      else {
        pico.out_w_trig() = pico.w_trig()*static_cast<float>(corr.w_zvtx_fail());
        pico.out_sys_trig()[0] = pico.sys_trig()[0]*static_cast<float>(corr.w_zvtx_fail());
        pico.out_sys_trig()[1] = pico.sys_trig()[1]*static_cast<float>(corr.w_zvtx_fail());
        pico.out_sys_trig_el()[0] = pico.sys_trig_el()[0]*static_cast<float>(corr.w_zvtx_fail());
        pico.out_sys_trig_el()[1] = pico.sys_trig_el()[1]*static_cast<float>(corr.w_zvtx_fail());
        pico.out_sys_trig_mu()[0] = pico.sys_trig_mu()[0]*static_cast<float>(corr.w_zvtx_fail());
        pico.out_sys_trig_mu()[1] = pico.sys_trig_mu()[1]*static_cast<float>(corr.w_zvtx_fail());
      }
    }
    else {
      pico.out_sys_trig()[0] = pico.sys_trig()[0];
      pico.out_sys_trig()[1] = pico.sys_trig()[1];
      pico.out_sys_trig_el()[0] = pico.sys_trig_el()[0];
      pico.out_sys_trig_el()[1] = pico.sys_trig_el()[1];
      pico.out_sys_trig_mu()[0] = pico.sys_trig_mu()[0];
      pico.out_sys_trig_mu()[1] = pico.sys_trig_mu()[1];
    }
    pico.out_w_lumi() = pico.w_lumi()>0 ? 1. : -1.; //get the generator weight sign
    pico.out_w_lumi() *= static_cast<float>(corr.w_lumi());

    pico.out_weight() = static_cast<float>(corr.weight()) * pico.out_w_lumi() *
                     pico.out_w_lep() * pico.out_w_fs_lep() * //post-corr values in order for 0l to be correct
                     btag_weight * pico.w_jetpuid() * pico.out_w_trig() * 
                     pico.out_w_isr() * pico.out_w_pu() * pico.w_prefire() * 
                     pico.w_photon() *  pico.out_w_phshape() * pico.w_fakephoton() * 
                     pico.out_w_nnlo();

    pico.out_sys_isr().resize(2); pico.out_sys_pu().resize(2);
    //pico.out_sys_photon().resize(2); pico.out_sys_photon_csev().resize(2);
    for (unsigned i(0); i<2; i++) {        
      pico.out_sys_pu()[i]             = pico.sys_pu()[i]*static_cast<float>(corr.sys_pu()[i]);
      pico.out_sys_isr()[i]            = pico.sys_isr()[i]*static_cast<float>(corr.sys_isr()[i]);
      //pico.out_sys_photon()[i]         = pico.sys_photon()[i]*static_cast<float>(corr.sys_photon()[i]);
      //pico.out_sys_photon_csev()[i]    = pico.sys_photon_csev()[i]*static_cast<float>(corr.sys_photon()[i]);

    } 
    if (is_zgamma) {
      for (unsigned i(0); i<2; i++) {        
        pico.out_sys_lep()[i] = pico.sys_el()[i]*pico.sys_mu()[i];
      }
    }
    pico.out_sys_murf().resize(9);
    for (unsigned i(0); i<pico.out_sys_murf().size(); i++) {        
      if (pico.sys_murf().size() != 0) {
        pico.out_sys_murf()[i]      = pico.sys_murf()[i]*static_cast<float>(corr.sys_murf()[i]);
      } else {
        pico.out_sys_murf()[i]      = 1.0;
      }
    }
    pico.out_sys_ps().resize(4);
    for (unsigned i(0); i<pico.out_sys_ps().size(); i++) {        
      if (pico.sys_ps().size() != 0) {
        pico.out_sys_ps()[i]      = pico.sys_ps()[i]*static_cast<float>(corr.sys_ps()[i]);
      } else {
        pico.out_sys_ps()[i]      = 1.0;
      }
    }
    //pico.out_sys_pdf().resize(102);
    //for (unsigned i(0); i<pico.out_sys_pdf().size(); i++) {        
    //  if (pico.sys_pdf().size() != 0) {
    //    pico.out_sys_pdf()[i]      = pico.sys_pdf()[i]*static_cast<float>(corr.sys_pdf()[i]);
    //  } else {
    //    pico.out_sys_pdf()[i]      = 1.0;
    //  }
    //}


    // for (unsigned i(0); i<pico.w_pdf().size(); i++) 
    //   pico.out_w_pdf()[i] = pico.w_pdf()[i]*corr.w_pdf()[i];    
    
    pico.Fill();

  } // loop over events
  
  pico.Write();

  cout<<endl;
  time(&endtime); 
  cout<<"Time passed: "<<hoursMinSec(difftime(endtime, begtime))<<endl<<endl;  
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"in_file", required_argument, 0, 'f'},  
      {"in_dir", required_argument, 0, 'i'},  
      {"corr_file", required_argument, 0, 'c'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "f:i:c:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'f':
      in_file = optarg;
      break;
    case 'i':
      in_dir = optarg;
      break;
    case 'c':
      corr_file = optarg;
      break;
    case 0:
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
