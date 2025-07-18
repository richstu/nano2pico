#include <vector>
#include <string>
#include <iostream>
#include <getopt.h>

#include "corrections_tree.hpp"
#include "cross_sections.hpp"
#include "utilities.hpp"

using namespace std;

template<typename T, typename U>
void CopySize(const vector<T> &wgt_sums, vector<U> &corr){
  corr = vector<U>(wgt_sums.size(), static_cast<U>(0.));
}

template<typename T, typename U>
void VecAdd(const vector<T> &wgt_sums, vector<U> &corr){
  for(size_t i = 0; i < wgt_sums.size(); ++i){
    corr.at(i) += wgt_sums.at(i);
  }
}

template<typename T>
void Normalize(T &x, double nent){
  x = x ? nent/x : 1.;
}

template<typename T>
void Normalize(vector<T> &v, double nent){
  for(auto &x: v) x = x ? nent/x : 1.;
}

void Initialize(corrections_tree &wgt_sums, corrections_tree &corr);
void AddEntry(corrections_tree &wgt_sums, corrections_tree &corr);
int GetHiggsinoMass(const string &path);
int GetGluinoMass(const string &path);
int GetSUSYMass(const string & patth, const string & tag);
void FixLumi(corrections_tree &corr, const string &corr_path, int year);
void FixISR(corrections_tree &corr, const string &corr_path, int year, bool is_zgamma);
void Normalize(corrections_tree &corr);
void Fix0L(corrections_tree &corr);
void FixZvtx(corrections_tree &corr, int year);

void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){

  GetOptions(argc, argv);

  if(argc-optind+1 < 3){
    cout << "Too few arguments! Usage: " << argv[0]
         << " output_file input_file [more_input_files...]" << endl;
    return 1;
  }
  
  string output_path = argv[optind];
  vector<string> input_paths(argv+optind+1, argv+argc);

  int year = 2018;
  if      (Contains(input_paths[0], "RunIISummer20UL16")) year = 2016;
  else if (Contains(input_paths[0], "RunIISummer20UL17")) year = 2017;
  else if (Contains(input_paths[0], "RunIISummer20UL18")) year = 2018;
  else if (Contains(input_paths[0], "RunIISummer16")) year = 2016;
  else if (Contains(input_paths[0], "RunIIFall17")) year = 2017;
  else if (Contains(input_paths[0], "RunIIAutumn18")) year = 2017;
  else if (Contains(input_paths[0], "Run3Summer22")) year = 2022;
  else if (Contains(input_paths[0], "Run3Summer23")) year = 2023;

  bool is_zgamma = Contains(output_path, "zgamma");
  //if (Contains(input_paths[0], "RunIISummer20")) {
  //  year = Contains(input_paths[0], "RunIISummer20UL16") ? 2016 : (Contains(input_paths[0], "RunIISummer20UL17") ? 2017 : 2018);
  //} else {
  //  year = Contains(input_paths[0], "RunIISummer16") ? 2016 : (Contains(input_paths[0], "RunIIFall17") ? 2017 : 2018);
  //}
  cout << "Running with settings for year = "<<year<<"."<<endl;

  corrections_tree corr("", output_path.c_str());
  corrections_tree wgt_sums(input_paths.front().c_str());
  for(size_t i = 1; i < input_paths.size(); ++i){
    wgt_sums.intree_->Add(input_paths.at(i).c_str());
  }

  size_t num_entries = wgt_sums.GetEntries();
  if(num_entries <= 0){
    cout << "No entries in input files!" << endl;
    return 1;
  } 
  wgt_sums.GetEntry(0);
  Initialize(wgt_sums, corr);

  for(size_t i = 0; i < num_entries; ++i){
    wgt_sums.GetEntry(i);
    AddEntry(wgt_sums, corr);
  }

  FixLumi(corr, output_path, year);
  FixISR(corr, output_path, year, is_zgamma);
  if (!is_zgamma)
    Fix0L(corr);
  FixZvtx(corr,year);

  Normalize(corr);

  corr.Fill();
  corr.Write();
  cout << "Wrote output to " << output_path << endl;
}

void Initialize(corrections_tree &wgt_sums, corrections_tree &corr){
  corr.out_weight() = 0.;
  corr.out_w_lumi() = 0.;
  corr.out_w_el() = 0.;
  corr.out_w_mu() = 0.;
  corr.out_w_lep() = 0.;
  corr.out_w_fs_lep() = 0.;
  corr.out_w_btag() = 0.;
  corr.out_w_btag_df() = 0.;
  corr.out_w_bhig() = 0.;
  corr.out_w_bhig_df() = 0.;
  corr.out_w_isr() = 0.;
  corr.out_w_pu() = 0.;
  corr.out_w_photon() = 0.;
  corr.out_w_phshape() = 0.;
  corr.out_w_trig() = 0.;
  corr.out_w_nnlo() = 0.;
  corr.out_w_zvtx_pass() = 0.;
  corr.out_w_zvtx_fail() = 0.;
  // w_prefire should not be normalized!!

  corr.out_neff() = 0;
  corr.out_neff_el() = 0;
  corr.out_neff_pass_eltrigs() = 0;
  corr.out_nent() = 0;
  corr.out_nent_zlep() = 0;
  corr.out_nent_isr() = 0;
  corr.out_tot_weight_l0() = 0.;
  corr.out_tot_weight_l1() = 0.;

  CopySize(wgt_sums.sys_el(),                 corr.out_sys_el());
  CopySize(wgt_sums.sys_mu(),                 corr.out_sys_mu());
  CopySize(wgt_sums.sys_lep(),                corr.out_sys_lep());
  CopySize(wgt_sums.sys_photon(),             corr.out_sys_photon());
  CopySize(wgt_sums.sys_photon_csev(),        corr.out_sys_photon_csev());
  CopySize(wgt_sums.sys_fs_lep(),             corr.out_sys_fs_lep());
  CopySize(wgt_sums.sys_bchig(),              corr.out_sys_bchig());
  CopySize(wgt_sums.sys_udsghig(),            corr.out_sys_udsghig());
  CopySize(wgt_sums.sys_fs_bchig(),           corr.out_sys_fs_bchig());
  CopySize(wgt_sums.sys_fs_udsghig(),         corr.out_sys_fs_udsghig());
  CopySize(wgt_sums.sys_trig(),               corr.out_sys_trig());
  CopySize(wgt_sums.sys_isr(),                corr.out_sys_isr());
  CopySize(wgt_sums.sys_pu(),                 corr.out_sys_pu());
  CopySize(wgt_sums.sys_murf(),               corr.out_sys_murf());
  CopySize(wgt_sums.sys_ps(),                 corr.out_sys_ps());
  // CopySize(wgt_sums.w_pdf(),                  corr.out_w_pdf());
  // CopySize(wgt_sums.sys_pdf(),                corr.out_sys_pdf());
}


void AddEntry(corrections_tree &wgt_sums, corrections_tree &corr){
  corr.out_neff() += wgt_sums.neff();
  corr.out_nent() += wgt_sums.nent();
  corr.out_nent_zlep() += wgt_sums.nent_zlep();
  corr.out_nent_isr() += wgt_sums.nent_isr();
  corr.out_neff_pass_eltrigs() += wgt_sums.neff_pass_eltrigs();
  corr.out_neff_el() += wgt_sums.neff_el();
  corr.out_tot_weight_l0() += wgt_sums.tot_weight_l0();
  corr.out_tot_weight_l1() += wgt_sums.tot_weight_l1();

  corr.out_weight()            += wgt_sums.weight();
  corr.out_w_el()              += wgt_sums.w_el();
  corr.out_w_mu()              += wgt_sums.w_mu();
  corr.out_w_lep()             += wgt_sums.w_lep();
  corr.out_w_fs_lep()          += wgt_sums.w_fs_lep();
  corr.out_w_bhig()            += wgt_sums.w_bhig();
  corr.out_w_btag()            += wgt_sums.w_btag();
  corr.out_w_bhig_df()         += wgt_sums.w_bhig_df();
  corr.out_w_btag_df()         += wgt_sums.w_btag_df();
  corr.out_w_trig()            += wgt_sums.w_trig();
  corr.out_w_isr()             += wgt_sums.w_isr();
  corr.out_w_pu()              += wgt_sums.w_pu();
  corr.out_w_photon()          += wgt_sums.w_photon();
  corr.out_w_phshape()         += wgt_sums.w_phshape();
  corr.out_w_nnlo()            += wgt_sums.w_nnlo();

  VecAdd(wgt_sums.sys_el(),            corr.out_sys_el());
  VecAdd(wgt_sums.sys_mu(),            corr.out_sys_mu());
  VecAdd(wgt_sums.sys_lep(),           corr.out_sys_lep());
  VecAdd(wgt_sums.sys_photon(),        corr.out_sys_photon());
  VecAdd(wgt_sums.sys_photon_csev(),   corr.out_sys_photon_csev());
  VecAdd(wgt_sums.sys_fs_lep(),        corr.out_sys_fs_lep());
  VecAdd(wgt_sums.sys_bchig(),         corr.out_sys_bchig());
  VecAdd(wgt_sums.sys_udsghig(),       corr.out_sys_udsghig());
  VecAdd(wgt_sums.sys_fs_bchig(),      corr.out_sys_fs_bchig());
  VecAdd(wgt_sums.sys_fs_udsghig(),    corr.out_sys_fs_udsghig());
  VecAdd(wgt_sums.sys_isr(),           corr.out_sys_isr());
  VecAdd(wgt_sums.sys_trig(),          corr.out_sys_trig());
  VecAdd(wgt_sums.sys_pu(),            corr.out_sys_pu());
  VecAdd(wgt_sums.sys_murf(),          corr.out_sys_murf());
  VecAdd(wgt_sums.sys_ps(),            corr.out_sys_ps());
  // VecAdd(wgt_sums.w_pdf(),             corr.out_w_pdf());
  // VecAdd(wgt_sums.sys_pdf(),           corr.out_sys_pdf());
}

int GetHiggsinoMass(const string &path){
  string key = "_mChi-";
  // if (Contains(path, "T2tt")) key = "_mStop-"; 
  auto pos1 = path.rfind(key)+key.size();
  auto pos2 = path.find("_", pos1);
  string mass_string = path.substr(pos1, pos2-pos1);
  int unrounded_mass = stoi(mass_string);
  int rounded_mass = unrounded_mass;
  if (unrounded_mass != 127)
    rounded_mass = ((unrounded_mass+12)/25)*25;
  return rounded_mass;
}

// Tag examples: "_mChi-","-mGluino-" 
int GetSUSYMass(const string & path, const string & tag) {
  auto pos1 = path.rfind(tag)+tag.size();
  auto pos2 = path.find("_", pos1);
  string mass_string = path.substr(pos1, pos2-pos1);
  int unrounded_mass = stoi(mass_string);
  return unrounded_mass;
}

int GetGluinoMass(const string & path) {
  return GetSUSYMass(path, "-mGluino-");
}

void FixLumi(corrections_tree &corr, const string &corr_path, int year){
  double xsec(0.); const float lumi = 1000.;
  if (Contains(corr_path, "SMS-TChiHH_HToAll")){
    double exsec(0.);
    int mglu = GetHiggsinoMass(corr_path);
    xsec::higgsinoCrossSection(mglu, xsec, exsec);
    xsec = xsec / .5824/.5824; // Remove H to bb branch ratio
    exsec = exsec / .5824/.5824; // Remove H to bb branch ratio
  } else if (Contains(corr_path, "SMS-TChiHH")&&Contains(corr_path, "HToGG")) {
    double exsec(0.);
    int mglu = GetHiggsinoMass(corr_path);
    xsec::higgsinoCrossSection(mglu, xsec, exsec);
    xsec = xsec / .5824; // Remove H to bb branch ratio
    exsec = exsec / .5824; // Remove H to bb branch ratio
  } else if (Contains(corr_path, "SMS-TChi")){
    double exsec(0.);
    int mglu = GetHiggsinoMass(corr_path);
    xsec::higgsinoCrossSection(mglu, xsec, exsec);
  } else if (Contains(corr_path, "SMS-T5qqqqZH_HToBB")) {
    double exsec(0.);
    int mglu = GetGluinoMass(corr_path);
    xsec::gluinoCrossSection(mglu, xsec, exsec);
    xsec = xsec * .5824*.5824; // Add in H to bb branch ratio
    exsec = exsec * .5824*.5824; // Add in H to bb branch ratio
  } else if (Contains(corr_path, "SMS-T5qqqqZH-")) {
    double exsec(0.);
    int mglu = GetGluinoMass(corr_path);
    xsec::gluinoCrossSection(mglu, xsec, exsec);
  }else{
    xsec = xsec::crossSection(corr_path, year);  
  }

  corr.out_w_lumi() = xsec*lumi/corr.out_neff();
}

void FixISR(corrections_tree &corr, const string &corr_path, int year, bool is_zgamma){
  if (is_zgamma) {
    corr.out_weight() = 1.0;
    double nent_isr = corr.out_nent_isr();
    Normalize(corr.out_w_isr(), nent_isr);
    return;
  }

  double corr_w_isr(1.);
  vector<double> corr_sys_isr(2,1.);
  double tot_w_isr = corr.out_w_isr();
  if (Contains(corr_path,"TTJets_HT") || Contains(corr_path,"genMET-150")){ // Full sim
  // in this case take correction from inclusive since should not norm. to unity
  //values consistent within 0.001 between 2016 and 2017 amazingly...
    if (Contains(corr_path,"TTJets_DiLept")) {
      corr_w_isr = 1/0.997;
      corr_sys_isr[0] = 1/1.057;
      corr_sys_isr[1] = 1/0.938;
    } else {
      corr_w_isr = 1/1.017;
      corr_sys_isr[0] = 1/1.067;
      corr_sys_isr[1] = 1/0.967;        
    }
  }else{
    // out_nent is sum of entries, out_w_isr is sum of isr weights
    // corr_w_isr scales out_w_isr to nent.
    corr_w_isr = corr.out_w_isr() ? corr.out_nent()/corr.out_w_isr() : 1.;
    for(size_t i = 0; i<corr.out_sys_isr().size(); i++){
      corr_sys_isr[i] = corr.out_sys_isr()[i] ? corr.out_nent()/corr.out_sys_isr()[i] : 1.;
    }
  }
  
  corr.out_w_isr() = corr_w_isr;
  for(size_t i = 0; i<corr.out_sys_isr().size(); i++){
    corr.out_sys_isr()[i] = corr_sys_isr[i];
  }
  double nent = corr.out_nent();
  double nent_zlep = corr.out_nent_zlep();

  // Calculate correction to total weight whilst correcting zero lepton
  //----------------------------------------------------------------------
  double w_corr_l0 = 1.;
  if (corr.out_w_lep()) w_corr_l0 *= (nent-corr.out_w_lep())/nent_zlep;
  if (corr.out_w_fs_lep()) w_corr_l0 *= (nent-corr.out_w_fs_lep())/nent_zlep;
  if(nent_zlep==0) w_corr_l0 = 1.;
  // again normalize to total w_isr, not unity
  if(year==2016) 
    corr.out_weight() = (tot_w_isr*corr_w_isr)/(corr.out_tot_weight_l0()*w_corr_l0 + corr.out_tot_weight_l1());
  else
    corr.out_weight() = nent/(corr.out_tot_weight_l0()*w_corr_l0 + corr.out_tot_weight_l1());

  //cout<<"tot_w_isr: "<<tot_w_isr<<" corr_w_isr: "<<corr_w_isr<<" nent: "<<corr.out_nent()<<" out_tot_weight_l0: "<<corr.out_tot_weight_l0()<<" w_corr_l0: "<<w_corr_l0<<" out_tot_weight_l1: "<<corr.out_tot_weight_l1()<<" out_weight: "<<corr.out_weight()<<endl;
}

void Fix0L(corrections_tree &corr){
  double nent = corr.out_nent();
  double nent_zlep = corr.out_nent_zlep();

  // Lepton weights corrections to be applied only to 0-lep events
  //----------------------------------------------------------------
  corr.out_w_lep()           = corr.out_w_lep() ? (nent-corr.out_w_lep())/nent_zlep : 1.;
  corr.out_w_fs_lep()        = corr.out_w_fs_lep() ? (nent-corr.out_w_fs_lep())/nent_zlep : 1.;
  for(size_t i = 0; i<corr.out_sys_lep().size(); i++){
    corr.out_sys_lep()[i]    = corr.out_sys_lep()[i] ? (nent-corr.out_sys_lep()[i])/nent_zlep : 1.;
  }
  for(size_t i = 0; i<corr.out_sys_fs_lep().size(); i++){
    corr.out_sys_fs_lep()[i] = corr.out_sys_fs_lep()[i] ? (nent-corr.out_sys_fs_lep()[i])/nent_zlep : 1.;
  }
}

/*!\brief calculate additional correction needed for Zvtx issue in 2017, even in 
  UL, see https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018. 
  Normalization is performed here
 */
void FixZvtx(corrections_tree &corr, int year) {
  //in apply corrections, this will only be applied to events failing dielectron
  //triggers
  if (year != 2017) {
    corr.out_w_zvtx_pass() = 1.0;
    corr.out_w_zvtx_fail() = 1.0;
  }
  else {
    float neff_pass = corr.out_neff_pass_eltrigs();
    float neff_fail = corr.out_neff_el()-corr.out_neff_pass_eltrigs();
    corr.out_w_zvtx_pass() = 0.991;
    corr.out_w_zvtx_fail() = (corr.out_neff_el()-neff_pass*0.991)/neff_fail;
  }
}



void Normalize(corrections_tree &corr){
  double nent = corr.out_nent();

  // total weight fixed in FixISR
  // w_lep fixed in Fix0L
  
  Normalize(corr.out_w_el(), nent);
  Normalize(corr.out_w_mu(), nent);

  Normalize(corr.out_w_btag(), nent);
  Normalize(corr.out_w_btag_df(), nent);
  Normalize(corr.out_w_bhig(), nent);
  Normalize(corr.out_w_bhig_df(), nent);

  // w_isr done in FixISR()
  Normalize(corr.out_w_pu(), nent);

  Normalize(corr.out_w_photon(), nent);
  Normalize(corr.out_w_phshape(), nent);
  Normalize(corr.out_w_nnlo(), nent);

  Normalize(corr.out_w_trig(), nent);

  Normalize(corr.out_sys_el(), nent);
  Normalize(corr.out_sys_mu(), nent);
  Normalize(corr.out_sys_photon(), nent);
  Normalize(corr.out_sys_photon_csev(), nent);
  Normalize(corr.out_sys_bchig(), nent);
  Normalize(corr.out_sys_udsghig(), nent);
  Normalize(corr.out_sys_fs_bchig(), nent);
  Normalize(corr.out_sys_fs_udsghig(), nent);
  Normalize(corr.out_sys_trig(), nent);
  Normalize(corr.out_sys_pu(), nent);

  Normalize(corr.out_sys_murf(), nent);
  Normalize(corr.out_sys_ps(), nent);
  // Normalize(corr.out_w_pdf(), nent);
  // Normalize(corr.out_sys_pdf(), nent);
}


void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}

