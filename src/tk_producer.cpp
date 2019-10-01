#include "tk_producer.hpp"

#include "utilities.hpp"

using namespace std;

IsoTrackProducer::IsoTrackProducer(int year_){
    year = year_;
}

IsoTrackProducer::~IsoTrackProducer(){
}

vector<int> IsoTrackProducer::WriteIsoTracks(nano_tree &nano, pico_tree &pico){
  
  vector<int> sig_tk_nano_idx;
  pico.out_ntk() = 0;
  for (int itk(0); itk < nano.nIsoTrack(); itk++) {
    //if (nano.IsoTrack_charge()[itk]==0) continue; // charge not stored in NanoAOD, request?
    unsigned int id = abs(nano.IsoTrack_pdgId()[itk]);
    if (id!=11 && id!=13 && id!=211) continue;
    if (id==211) {
      if (nano.IsoTrack_pt()[itk] < HadronIsoTrackPtCut) continue;
      if (nano.IsoTrack_pfRelIso03_chg()[itk] < HadronIsoTrackRelIsoCut) continue;
    } else {
      if (nano.IsoTrack_pt()[itk] < LeptonIsoTrackPtCut) continue;
      if (nano.IsoTrack_pfRelIso03_chg()[itk] < LeptonIsoTrackRelIsoCut) continue;
    }

    if (fabs(nano.IsoTrack_eta()[itk]) > IsoTrackEtaCut) continue;
    if (fabs(nano.IsoTrack_dz()[itk]) > IsoTrackDzCut) continue;
    float mt_ = GetMT(nano.MET_pt(), nano.MET_phi(),  nano.IsoTrack_pt()[itk], nano.IsoTrack_phi()[itk]);
    if (IsoTrackMtCut>0 && mt_>IsoTrackMtCut) continue;

    pico.out_tk_pdgid().push_back(nano.IsoTrack_pdgId()[itk]);
    pico.out_tk_pt().push_back(nano.IsoTrack_pt()[itk]);
    pico.out_tk_eta().push_back(nano.IsoTrack_eta()[itk]);
    pico.out_tk_phi().push_back(nano.IsoTrack_phi()[itk]);
    pico.out_tk_d0().push_back(nano.IsoTrack_dxy()[itk]);
    pico.out_tk_miniso_chg().push_back(nano.IsoTrack_miniPFRelIso_chg()[itk]);
    pico.out_tk_reliso_chg().push_back(nano.IsoTrack_pfRelIso03_chg()[itk]);
    pico.out_tk_dz().push_back(nano.IsoTrack_dz()[itk]);
    pico.out_tk_mt().push_back(mt_);

    // filled by mc_producer
    pico.out_tk_tm().push_back(false);    

    sig_tk_nano_idx.push_back(itk);
    pico.out_ntk()++;
  }
  
  return sig_tk_nano_idx;
}
