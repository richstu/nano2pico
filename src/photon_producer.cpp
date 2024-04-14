#include "photon_producer.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "correction.hpp"
#include "utilities.hpp"

using namespace std;

PhotonProducer::PhotonProducer(int year_, bool isData_, bool preVFP, float nanoaod_version_){
  year = year_;
  isData = isData_;
  nanoaod_version = nanoaod_version_;
  if (year==2016 && preVFP) {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2016preVFP_UL/EGM_ScaleUnc.json");
    str_scale_syst_ = "2016preVFP";
  }
  else if (year==2016) {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2016postVFP_UL/EGM_ScaleUnc.json");
    str_scale_syst_ = "2016postVFP";
  }
  else if (year==2017) {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2017_UL/EGM_ScaleUnc.json");
    str_scale_syst_ = "2017";
  }
  else if (year==2018) {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2018_UL/EGM_ScaleUnc.json");
    str_scale_syst_ = "2018";
  }
  else {
    cs_scale_syst_ = correction::CorrectionSet::from_file(
        "data/zgamma/2018_UL/EGM_ScaleUnc.json");
    str_scale_syst_ = "2018";
    std::cout << "WARNING: No dedicated EGM scale/smearing JSONs, defaulting to 2018" << std::endl;
  }
  map_scale_syst_ = cs_scale_syst_->at("UL-EGM_ScaleUnc");
}

PhotonProducer::~PhotonProducer(){
}

vector<int> PhotonProducer::WritePhotons(nano_tree &nano, pico_tree &pico, vector<int> &jet_isphoton_nano_idx, vector<int> &sig_el_nano_idx, vector<int> &sig_mu_nano_idx, vector<int> &photon_el_pico_idx){
  pico.out_nphoton() = 0;
  // pico.out_nfsrphoton() = 0;
  vector<int> sig_photon_nano_idx;
  int nphotons(0), ndr(0), shift(0);

  vector<int> FsrPhoton_muonIdx;
  getFsrPhoton_muonIdx(nano, nanoaod_version, FsrPhoton_muonIdx);
  vector<int> Photon_jetIdx;
  getPhoton_jetIdx(nano, nanoaod_version, Photon_jetIdx);
  vector<int> Photon_cutBased;
  getPhoton_cutBased(nano, nanoaod_version, Photon_cutBased);
 
  for(int iph(0); iph<nano.nPhoton(); ++iph){
    float pt = nano.Photon_pt()[iph];
    float eta = nano.Photon_eta()[iph];
    float phi = nano.Photon_phi()[iph];
    float mva = nano.Photon_mvaID()[iph];
    bool eVeto = nano.Photon_electronVeto()[iph];

    if (pt <= PhotonPtCut) continue;
    if (!(nano.Photon_isScEtaEB()[iph] || nano.Photon_isScEtaEE()[iph])) continue;

    // Find min(dR) between photon and signal lepton
    double minLepDR(999.);
    double maxLepDR(0.);
    for(size_t iel(0); iel<sig_el_nano_idx.size(); iel++) {
      double tempDR = dR(eta, nano.Electron_eta()[sig_el_nano_idx.at(iel)],
                         phi, nano.Electron_phi()[sig_el_nano_idx.at(iel)]);
      if(tempDR < minLepDR) minLepDR = tempDR;
      if(tempDR > maxLepDR) maxLepDR = tempDR;
    }
    for(size_t imu(0); imu<sig_mu_nano_idx.size(); imu++) {
      double tempDR = dR(eta, nano.Muon_eta()[sig_mu_nano_idx.at(imu)],
                         phi, nano.Muon_phi()[sig_mu_nano_idx.at(imu)]);
      if(tempDR < minLepDR) minLepDR = tempDR;
      if(tempDR > maxLepDR) maxLepDR = tempDR;
    }

    bool isSignal = (nano.Photon_mvaID_WP80()[iph] &&
                    eVeto && minLepDR > 0.3 && 
                    pt > SignalPhotonPtCut &&
                    (photon_el_pico_idx[iph]==-1 || !(pico.out_el_sig()[photon_el_pico_idx[iph]])));

    // Photons passing the object selections are placed at the front
    if(isSignal) {
      shift = ndr;
      ndr++;
    }
    else
      shift = nphotons;

    float scale_syst_up = 1.0;
    float scale_syst_dn = 1.0;
    if (year <= 2018) {
      scale_syst_up = map_scale_syst_->evaluate({str_scale_syst_,"scaleup",eta,
          nano.Photon_seedGain()[iph]});
      scale_syst_dn = map_scale_syst_->evaluate({str_scale_syst_,"scaledown",eta,
          nano.Photon_seedGain()[iph]});
    }

    float origin_eta = 0.0;
    if (nano.Photon_isScEtaEB()[iph]) {
      float pv_tan_theta_over_2 = exp(-1.0*eta);
      float pv_tan_theta = 2.0*pv_tan_theta_over_2/(1.0-pv_tan_theta_over_2*pv_tan_theta_over_2);
      float photon_unit_x = cos(nano.Photon_phi()[iph]);
      float photon_unit_y = sin(nano.Photon_phi()[iph]);
      float pv_ecal_dr = 130.0 - (nano.PV_x()*photon_unit_x+nano.PV_y()*photon_unit_y);
      float pv_ecal_dz = pv_ecal_dr/pv_tan_theta;
      float origin_theta = atan(130.0/(nano.PV_z()+pv_ecal_dz));
      if (origin_theta < 0) origin_theta += M_PI;
      origin_eta = -1.0*log(tan(origin_theta/2.0));
    }
    else { //if (nano.Photon_isScEtaEE()[iph]) {
      float pv_tan_theta_over_2 = exp(-1.0*eta);
      float pv_tan_theta = 2.0*pv_tan_theta_over_2/(1.0-pv_tan_theta_over_2*pv_tan_theta_over_2);
      float photon_unit_x = cos(nano.Photon_phi()[iph]);
      float photon_unit_y = sin(nano.Photon_phi()[iph]);
      float pv_ecal_dz = 310.0-nano.PV_z(); //+ endcap
      if (eta < 0) pv_ecal_dz = 310.0+nano.PV_z(); //- endcap
      float pv_ecal_dr = pv_ecal_dz*pv_tan_theta;
      float origin_theta = atan(((photon_unit_x*nano.PV_x()+photon_unit_y*nano.PV_y())+pv_ecal_dr)/310.0);
      if (origin_theta < 0) origin_theta += M_PI;
      origin_eta = -1.0*log(tan(origin_theta/2.0));
    }
    
    switch(year) {
      case 2016:
      case 2017:
      case 2018:
        pico.out_photon_pt()    .insert(pico.out_photon_pt()    .begin()+shift, pt);
        pico.out_photon_eta()   .insert(pico.out_photon_eta()   .begin()+shift, eta);
        pico.out_photon_phi()   .insert(pico.out_photon_phi()   .begin()+shift, phi);
        pico.out_photon_ecorr() .insert(pico.out_photon_ecorr() .begin()+shift, nano.Photon_eCorr()[iph]);
        pico.out_photon_reliso().insert(pico.out_photon_reliso().begin()+shift, nano.Photon_pfRelIso03_all()[iph]);
        pico.out_photon_r9()    .insert(pico.out_photon_r9()    .begin()+shift, nano.Photon_r9()[iph]);
        pico.out_photon_sieie() .insert(pico.out_photon_sieie() .begin()+shift, nano.Photon_sieie()[iph]);
        pico.out_photon_energyErr().insert(pico.out_photon_energyErr() .begin()+shift, nano.Photon_energyErr()[iph]);
        pico.out_photon_hoe()   .insert(pico.out_photon_hoe()   .begin()+shift, nano.Photon_hoe()[iph]);
        pico.out_photon_elveto().insert(pico.out_photon_elveto().begin()+shift, eVeto);
        pico.out_photon_id()    .insert(pico.out_photon_id()    .begin()+shift, nano.Photon_mvaID_WP90()[iph]);
        pico.out_photon_id80()  .insert(pico.out_photon_id80()  .begin()+shift, nano.Photon_mvaID_WP80()[iph]);
        pico.out_photon_idmva() .insert(pico.out_photon_idmva() .begin()+shift, mva);
        pico.out_photon_idCutBased().insert(pico.out_photon_idCutBased().begin()+shift, Photon_cutBased[iph]);
        pico.out_photon_idCutBasedBitMap().insert(pico.out_photon_idCutBasedBitMap() .begin()+shift, nano.Photon_vidNestedWPBitmap()[iph]);
        pico.out_photon_origin_eta().insert(pico.out_photon_origin_eta().begin()+shift, origin_eta);
        pico.out_photon_isScEtaEB().insert(pico.out_photon_isScEtaEB().begin()+shift, nano.Photon_isScEtaEB()[iph]);
        pico.out_photon_isScEtaEE().insert(pico.out_photon_isScEtaEE().begin()+shift, nano.Photon_isScEtaEE()[iph]);
        pico.out_photon_sig()   .insert(pico.out_photon_sig()   .begin()+shift, isSignal);
        pico.out_photon_drmin() .insert(pico.out_photon_drmin() .begin()+shift, minLepDR);
        pico.out_photon_drmax() .insert(pico.out_photon_drmax() .begin()+shift, maxLepDR);
        pico.out_photon_elidx() .insert(pico.out_photon_elidx() .begin()+shift, photon_el_pico_idx[iph]);
        pico.out_photon_pixelseed().insert(pico.out_photon_pixelseed().begin()+shift,nano.Photon_pixelSeed()[iph]);
        pico.out_sys_photon_pt_resup().insert(pico.out_sys_photon_pt_resup().begin()+shift, pt*(1.0+nano.Photon_dEsigmaUp()[iph]));
        pico.out_sys_photon_pt_resdn().insert(pico.out_sys_photon_pt_resdn().begin()+shift, pt*(1.0+nano.Photon_dEsigmaDown()[iph]));
        pico.out_sys_photon_pt_scaleup().insert(pico.out_sys_photon_pt_scaleup().begin()+shift, pt*scale_syst_up);
        pico.out_sys_photon_pt_scaledn().insert(pico.out_sys_photon_pt_scaledn().begin()+shift, pt*scale_syst_dn);
        break;
      case 2022:
      case 2023:
        pico.out_photon_pt()      .insert(pico.out_photon_pt()      .begin()+shift, pt);
        pico.out_photon_eta()     .insert(pico.out_photon_eta()     .begin()+shift, eta);
        pico.out_photon_phi()     .insert(pico.out_photon_phi()     .begin()+shift, phi);
        pico.out_photon_reliso()  .insert(pico.out_photon_reliso()  .begin()+shift, nano.Photon_pfRelIso03_all_quadratic()[iph]);
        pico.out_photon_phiso()   .insert(pico.out_photon_phiso()   .begin()+shift, nano.Photon_pfPhoIso03()[iph]);
        pico.out_photon_chiso()   .insert(pico.out_photon_chiso()   .begin()+shift, nano.Photon_pfChargedIsoPFPV()[iph]);
        pico.out_photon_chiso_worst().insert(pico.out_photon_chiso_worst().begin()+shift,nano.Photon_pfChargedIsoWorstVtx()[iph]);
        pico.out_photon_r9()      .insert(pico.out_photon_r9()      .begin()+shift, nano.Photon_r9()[iph]);
        pico.out_photon_s4()      .insert(pico.out_photon_s4()      .begin()+shift, nano.Photon_s4()[iph]);
        pico.out_photon_sieie()   .insert(pico.out_photon_sieie()   .begin()+shift, nano.Photon_sieie()[iph]);
        pico.out_photon_sieip()   .insert(pico.out_photon_sieip()   .begin()+shift, nano.Photon_sieip()[iph]);
        pico.out_photon_etawidth().insert(pico.out_photon_etawidth().begin()+shift, nano.Photon_etaWidth()[iph]);
        pico.out_photon_phiwidth().insert(pico.out_photon_phiwidth().begin()+shift, nano.Photon_phiWidth()[iph]);
        pico.out_photon_energyErr().insert(pico.out_photon_energyErr() .begin()+shift, nano.Photon_energyErr()[iph]);
        pico.out_photon_hoe()     .insert(pico.out_photon_hoe()     .begin()+shift, nano.Photon_hoe()[iph]);
        pico.out_photon_energy_raw().insert(pico.out_photon_energy_raw().begin()+shift, nano.Photon_energyRaw()[iph]);
        pico.out_photon_esoversc().insert(pico.out_photon_esoversc().begin()+shift, nano.Photon_esEnergyOverRawE()[iph]);
        pico.out_photon_essigmarr().insert(pico.out_photon_essigmarr().begin()+shift,nano.Photon_esEffSigmaRR()[iph]);
        pico.out_photon_elveto()  .insert(pico.out_photon_elveto()  .begin()+shift, eVeto);
        pico.out_photon_id()      .insert(pico.out_photon_id()      .begin()+shift, nano.Photon_mvaID_WP90()[iph]);
        pico.out_photon_id80()    .insert(pico.out_photon_id80()    .begin()+shift, nano.Photon_mvaID_WP80()[iph]);
        pico.out_photon_idmva()   .insert(pico.out_photon_idmva()   .begin()+shift, mva);
        pico.out_photon_idCutBased().insert(pico.out_photon_idCutBased() .begin()+shift, Photon_cutBased[iph]);
        pico.out_photon_idCutBasedBitMap().insert(pico.out_photon_idCutBasedBitMap() .begin()+shift, nano.Photon_vidNestedWPBitmap()[iph]);
        pico.out_photon_origin_eta().insert(pico.out_photon_origin_eta().begin()+shift, origin_eta);
        pico.out_photon_isScEtaEB().insert(pico.out_photon_isScEtaEB().begin()+shift, nano.Photon_isScEtaEB()[iph]);
        pico.out_photon_isScEtaEE().insert(pico.out_photon_isScEtaEE().begin()+shift, nano.Photon_isScEtaEE()[iph]);
        pico.out_photon_sig()     .insert(pico.out_photon_sig()     .begin()+shift, isSignal);
        pico.out_photon_drmin()   .insert(pico.out_photon_drmin()   .begin()+shift, minLepDR);
        pico.out_photon_drmax()   .insert(pico.out_photon_drmax()   .begin()+shift, maxLepDR);
        pico.out_photon_elidx()   .insert(pico.out_photon_elidx()   .begin()+shift, photon_el_pico_idx[iph]);
        pico.out_photon_pixelseed().insert(pico.out_photon_pixelseed().begin()+shift,nano.Photon_pixelSeed()[iph]);
        break;
        //TODO: add correctionlib-based scale/smearing corrections for run 3 (not available yet?)
      default:
        std::cout<<"Need code for new year in getZGammaPhBr (in photon_producer.cpp)"<<endl;
        exit(1);
    }

    if (nanoaod_version > 9.49 && nanoaod_version < 9.51) { //custom NanoAOD production
      pico.out_photon_phiso()   .insert(pico.out_photon_phiso()   .begin()+shift, nano.Photon_pfPhoIso03()[iph]);
      pico.out_photon_chiso()   .insert(pico.out_photon_chiso()   .begin()+shift, nano.Photon_pfChargedIso()[iph]); //note different name from NanoAODv12
      pico.out_photon_chiso_worst().insert(pico.out_photon_chiso_worst().begin()+shift,nano.Photon_pfChargedIsoWorstVtx()[iph]);
      pico.out_photon_s4()      .insert(pico.out_photon_s4()      .begin()+shift, nano.Photon_s4()[iph]);
      pico.out_photon_sieip()   .insert(pico.out_photon_sieip()   .begin()+shift, nano.Photon_sieip()[iph]);
      pico.out_photon_etawidth().insert(pico.out_photon_etawidth().begin()+shift, nano.Photon_etaWidth()[iph]);
      pico.out_photon_phiwidth().insert(pico.out_photon_phiwidth().begin()+shift, nano.Photon_phiWidth()[iph]);
      pico.out_photon_energy_raw().insert(pico.out_photon_energy_raw().begin()+shift, nano.Photon_energyRaw()[iph]);
      pico.out_photon_esoversc().insert(pico.out_photon_esoversc().begin()+shift, nano.Photon_esEnergyOverRawE()[iph]);
      pico.out_photon_essigmarr().insert(pico.out_photon_essigmarr().begin()+shift,nano.Photon_esEffSigmaRR()[iph]);
    }

    nphotons++;
    if (!isData)
      pico.out_photon_pflavor().insert(pico.out_photon_pflavor().begin()+shift, nano.Photon_genPartFlav()[iph]);

    // All photons with pt > 15 GeV are considered for creating the ZGamma candidates
    if(isSignal) {
      pico.out_nphoton()++;
      sig_photon_nano_idx.push_back(iph);
      // save indices of matching jets
      if (Photon_jetIdx[iph]>=0)
        jet_isphoton_nano_idx.push_back(Photon_jetIdx[iph]);
      else
        for (int ijet(0); ijet<nano.nJet(); ijet++)
          if (dR(eta, nano.Jet_eta()[ijet], phi, nano.Jet_phi()[ijet])<0.4f)
            jet_isphoton_nano_idx.push_back(ijet);
    }
  }

  //Saves indices for FSR photons
  pico.out_nfsrphoton() = 0;
  for(int iph(0); iph < nano.nFsrPhoton(); ++iph){
    //These set of selections require the photon be close to the corresponding muon and relatively soft
    if (nano.FsrPhoton_pt()[iph] <= FsrPhotonPtCut) continue;
    if (fabs(nano.FsrPhoton_eta()[iph]) > FsrPhotonEtaCut) continue;
    if (nano.FsrPhoton_relIso03()[iph] > FsrPhotonIsoCut) continue;
    if (nano.FsrPhoton_dROverEt2()[iph] > FsrPhotondRCut) continue;

    //Check for separation between fsrphoton and first photon candidate
    if ( dR(pico.out_photon_eta()[0],nano.FsrPhoton_eta()[iph],pico.out_photon_phi()[0],nano.FsrPhoton_phi()[iph]) < FsrSeparationReq) continue;
    
    //Add the values to the pico trees
    pico.out_fsrphoton_pt().push_back(nano.FsrPhoton_pt()[iph]);
    pico.out_fsrphoton_eta().push_back(nano.FsrPhoton_eta()[iph]);
    pico.out_fsrphoton_phi().push_back(nano.FsrPhoton_phi()[iph]);
    pico.out_fsrphoton_reliso().push_back(nano.FsrPhoton_relIso03()[iph]);
    pico.out_fsrphoton_muonidx().push_back(FsrPhoton_muonIdx[iph]);
    pico.out_fsrphoton_droveret2().push_back(nano.FsrPhoton_dROverEt2()[iph]);
    pico.out_nfsrphoton()++;
  } 

  return sig_photon_nano_idx;
}

bool PhotonProducer::idPhoton(int bitmap, int level){
  // decision for each cut represented by 1 bit
  //0 - MinPtCut
  //1 - PhoSCEtaMultiRangeCut
  //2 - PhoSingleTowerHadOverEmCut
  //3 - PhoFull5x5SigmaIEtaIEtaCut
  //4 - PhoGenericRhoPtScaledCut
  //5 - PhoGenericRhoPtScaledCut
  //6 - PhoGenericRhoPtScaledCut
  bool pass = true;
  for (int i(0); i<7; i++){
    if (((bitmap >> i*2) & 0x3) < level) pass = false;
  }
  return pass;
}

