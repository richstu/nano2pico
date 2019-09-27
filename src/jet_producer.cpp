#include "jet_producer.hpp"

#include "utilities.hpp"

using namespace std;

JetProducer::JetProducer(int year_){
    year = year_;
}

JetProducer::~JetProducer(){
}

void JetProducer::WriteJets(nano_tree &nano, pico_tree &pico){
    pico.out_njet() = 0;
    for(int ijet(0); ijet<nano.nJet(); ++ijet){
      if (nano.Jet_pt().at(ijet) > JetPtCut) {
        pico.out_jet_pt().push_back(nano.Jet_pt().at(ijet));
        pico.out_njet()++;
      }
    }

    return;
}
