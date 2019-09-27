#include "tk_producer.hpp"

#include "utilities.hpp"

using namespace std;

IsoTrackProducer::IsoTrackProducer(int year_){
    year = year_;
}

IsoTrackProducer::~IsoTrackProducer(){
}

void IsoTrackProducer::WriteIsoTracks(nano_tree &nano, pico_tree &pico){
    pico.out_ntk() = 0;
    for(int itk(0); itk<nano.nIsoTrack(); ++itk){
        pico.out_tk_pt().push_back(nano.IsoTrack_pt()[itk]);
        pico.out_ntk()++;
    }

    return;
}
