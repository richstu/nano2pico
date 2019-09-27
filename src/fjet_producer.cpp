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
      if (nano.FatJet_pt().at(ijet) > FatJetPtCut) {
        pico.out_fjet_pt().push_back(nano.FatJet_pt().at(ijet));
        pico.out_nfjet()++;
      }
    }

    return;
}
