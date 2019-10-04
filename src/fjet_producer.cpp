#include "fjet_producer.hpp"

#include "utilities.hpp"

using namespace std;

FatJetProducer::FatJetProducer(int year_){
    year = year_;
}

FatJetProducer::~FatJetProducer(){
}

void FatJetProducer::WriteFatJets(nano_tree &nano, pico_tree &pico){

  pico.out_nfjet() = 0; 
  for(int ijet(0); ijet<nano.nFatJet(); ++ijet){
    if (nano.FatJet_pt()[ijet] <= FatJetPtCut) continue;
    if (fabs(nano.FatJet_eta()[ijet]) > FatJetEtaCut) continue;

    pico.out_fjet_pt().push_back(nano.FatJet_pt()[ijet]);
    pico.out_fjet_eta().push_back(nano.FatJet_eta()[ijet]);
    pico.out_fjet_phi().push_back(nano.FatJet_phi()[ijet]);
    pico.out_fjet_m().push_back(nano.FatJet_mass()[ijet]);
    // Mass-decorrelated Deep Double B, H->bb vs QCD discriminator, endorsed by BTV
    pico.out_fjet_md_hbb_btv().push_back(nano.FatJet_btagDDBvL()[ijet]);
    // Mass-decorrelated DeepAk8, H->bb vs QCD discriminator, endorsed by JME
    pico.out_fjet_md_hbb_jme().push_back(nano.FatJet_deepTagMD_HbbvsQCD()[ijet]);

    //this will be filled in hig_producer or later
    pico.out_fjet_h1d().push_back(false);
    pico.out_fjet_h2d().push_back(false);

    pico.out_nfjet()++;
  }
  return;
}
