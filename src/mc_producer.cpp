#include "mc_producer.hpp"

#include "utilities.hpp"

using namespace std;

GenParticleProducer::GenParticleProducer(int year_){
    year = year_;
}

GenParticleProducer::~GenParticleProducer(){
}

void GenParticleProducer::WriteGenParticles(nano_tree &nano, pico_tree &pico){

    for(int imc(0); imc<nano.nGenPart(); ++imc){
      if (nano.GenPart_status().at(imc) == 23) {
        pico.out_mc_id().push_back(nano.GenPart_pdgId().at(imc));
      }
    }

    return;
}
