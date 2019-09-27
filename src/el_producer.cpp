#include "el_producer.hpp"

#include "utilities.hpp"

using namespace std;

ElectronProducer::ElectronProducer(int year_){
    year = year_;
}

ElectronProducer::~ElectronProducer(){
}

void ElectronProducer::WriteElectrons(nano_tree &nano, pico_tree &pico){
    pico.out_nel() = 0;
    for(int iel(0); iel<nano.nElectron(); ++iel){
      if (nano.Electron_pt().at(iel) > SignalElectronPtCut) {
        pico.out_el_pt().push_back(nano.Electron_pt().at(iel));
        pico.out_nel()++;
      }
    }

    return;
}
