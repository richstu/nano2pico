#include "zgamma_producer.hpp"
#include "utilities.hpp"
#include "KinZfitter.hpp"
#include "TLorentzVector.h"

using namespace std;

ZGammaVarProducer::ZGammaVarProducer(int year_){
    year = year_;
    kinZfitter = new KinZfitter();
}

ZGammaVarProducer::~ZGammaVarProducer(){
    delete kinZfitter;

}

void ZGammaVarProducer::WriteZGammaVars(nano_tree &nano, pico_tree &pico, vector<int> sig_jet_nano_idx){
  if(sig_jet_nano_idx.size() > 1) {
    TLorentzVector j1, j2, dijet;
    j1.SetPtEtaPhiM(nano.Jet_pt()[sig_jet_nano_idx[0]], nano.Jet_eta()[sig_jet_nano_idx[0]], nano.Jet_phi()[sig_jet_nano_idx[0]],nano.Jet_mass()[sig_jet_nano_idx[0]]);
    j2.SetPtEtaPhiM(nano.Jet_pt()[sig_jet_nano_idx[1]], nano.Jet_eta()[sig_jet_nano_idx[1]], nano.Jet_phi()[sig_jet_nano_idx[1]],nano.Jet_mass()[sig_jet_nano_idx[1]]);
    dijet = j1 + j2;
    pico.out_dijet_pt()   = dijet.Pt();
    pico.out_dijet_eta()  = dijet.Eta();
    pico.out_dijet_phi()  = dijet.Phi();
    pico.out_dijet_m()    = dijet.M();
    pico.out_dijet_dr()   = j1.DeltaR(j2);
    pico.out_dijet_dphi() = j1.DeltaPhi(j2);
    pico.out_dijet_deta() = fabs(j1.Eta() - j2.Eta());
  }


//Generate pT cut summary branches here so they are useful even for unskimmed files
  pico.out_trig_el_pt() = false;
  pico.out_trig_mu_pt() = false;
  if (year==2016) {
    if(pico.out_el_pt().size() > 1){
      if((pico.out_trig_double_el() && pico.out_el_pt().at(0)>25.f && pico.out_el_pt().at(1)>15.f) || (pico.out_trig_single_el() && pico.out_el_pt().at(0)>30.f)){
        pico.out_trig_el_pt() = true;
      }
    } else if(pico.out_el_pt().size() > 0){
      if(pico.out_trig_single_el() && pico.out_el_pt().at(0)>30.f){
        pico.out_trig_el_pt() = true;
      }
    }

    if(pico.out_mu_pt().size() > 1){
      if((pico.out_trig_double_mu() && pico.out_mu_pt().at(0)>20.f && pico.out_mu_pt().at(1)>10.f) || (pico.out_trig_single_mu() && pico.out_mu_pt().at(0)>25.f)){
        pico.out_trig_mu_pt() = true;
      }
    }else if(pico.out_mu_pt().size() > 0){
      if(pico.out_trig_single_mu() && pico.out_mu_pt().at(0)>25.f){
        pico.out_trig_mu_pt() = true;
      }
    }
  }

  if (year==2017) {
    if(pico.out_el_pt().size() > 1){
      if((pico.out_trig_double_el() && pico.out_el_pt().at(0)>25.f && pico.out_el_pt().at(1)>15.f) || (pico.out_trig_single_el() && pico.out_el_pt().at(0)>35.f)){
        pico.out_trig_el_pt() = true;
      }
    } else if(pico.out_el_pt().size() > 0){
      if(pico.out_trig_single_el() && pico.out_el_pt().at(0)>35.f){
        pico.out_trig_el_pt() = true;
      }
    }

    if(pico.out_mu_pt().size() > 1){
      if((pico.out_trig_double_mu() && pico.out_mu_pt().at(0)>20.f && pico.out_mu_pt().at(1)>10.f) || (pico.out_trig_single_mu() && pico.out_mu_pt().at(0)>28.f)){
        pico.out_trig_mu_pt() = true;
      }
    }else if(pico.out_mu_pt().size() > 0){
      if(pico.out_trig_single_mu() && pico.out_mu_pt().at(0)>28.f){
        pico.out_trig_mu_pt() = true;
      }
    }
  }

  if (year==2018 || year==2022 || year==2023) {
    if(pico.out_el_pt().size() > 1){
      if((pico.out_trig_double_el() && pico.out_el_pt().at(0)>25.f && pico.out_el_pt().at(1)>15.f) || (pico.out_trig_single_el() && pico.out_el_pt().at(0)>35.f)){
        pico.out_trig_el_pt() = true;
      }
    } else if(pico.out_el_pt().size() > 0){
      if(pico.out_trig_single_el() && pico.out_el_pt().at(0)>35.f){
        pico.out_trig_el_pt() = true;
      }
    }

    if(pico.out_mu_pt().size() > 1){
      if((pico.out_trig_double_mu() && pico.out_mu_pt().at(0)>20.f && pico.out_mu_pt().at(1)>10.f) || (pico.out_trig_single_mu() && pico.out_mu_pt().at(0)>25.f)){
        pico.out_trig_mu_pt() = true;
      }
    }else if(pico.out_mu_pt().size() > 0){
      if(pico.out_trig_single_mu() && pico.out_mu_pt().at(0)>25.f){
        pico.out_trig_mu_pt() = true;
      }
    }
  }


  pico.out_nllphoton() = 0;
  int baseBit = 0b000000000000;//If 0, no selections applied. Placed here to get 0 for cases where no dilep or photon
  int categoryBit = 0b00000000;


  if(pico.out_ll_pt().size()!=0 && pico.out_nphoton()==0){//Quick hacky fix to allow bitmaps for cases without a photon (unskimmed picos). Please consolidate bitmap tools to a separate file for future productions
    if(pico.out_nel()>=2 && pico.out_ll_lepid().at(0)==11){
      baseBit += 0b110000000000;
      if(pico.out_trig_single_el() || pico.out_trig_double_el()){baseBit +=0b000100000000;}
      if(pico.out_trig_el_pt()){baseBit += 0b000010000000;}
    }
    if(pico.out_nmu()>=2 && pico.out_ll_lepid().at(0)==13){
      baseBit += 0b101000000000;
      if(pico.out_trig_single_mu() || pico.out_trig_double_mu()){baseBit +=0b000100000000;}
      if(pico.out_trig_mu_pt()){baseBit +=0b000010000000;}
    }
    pico.out_zg_cutBitMap() = baseBit;
  }

  if (pico.out_ll_pt().size() == 0 || pico.out_nphoton() == 0){
 
    pico.out_zg_cutBitMap() = baseBit;
    pico.out_zg_categorizationBitMap() = categoryBit;
    return;
  }



  for(size_t ill(0); ill < pico.out_ll_pt().size(); ill++){
    for(size_t igamma(0); igamma < pico.out_photon_pt().size(); igamma++) {
      TLorentzVector dilep, dilep_refit, photon, llg;
      dilep.SetPtEtaPhiM(pico.out_ll_pt()[ill], pico.out_ll_eta()[ill],
                         pico.out_ll_phi()[ill],pico.out_ll_m()[ill]);

      photon.SetPtEtaPhiM(pico.out_photon_pt()[igamma], pico.out_photon_eta()[igamma], pico.out_photon_phi()[igamma], 0);
      llg = dilep + photon;
      pico.out_nllphoton()++;
      pico.out_llphoton_pt().push_back(llg.Pt());
      pico.out_llphoton_eta().push_back(llg.Eta());
      pico.out_llphoton_phi().push_back(llg.Phi());
      pico.out_llphoton_m().push_back(llg.M());
      pico.out_llphoton_dr().push_back(dilep.DeltaR(photon));
      pico.out_llphoton_dphi().push_back(dilep.DeltaPhi(photon));
      pico.out_llphoton_deta().push_back(fabs(dilep.Eta() - photon.Eta()));
      pico.out_llphoton_iph().push_back(igamma);
      pico.out_llphoton_ill().push_back(ill);


      if(sig_jet_nano_idx.size() > 1) {
        TLorentzVector j1, j2, dijet;
        j1.SetPtEtaPhiM(nano.Jet_pt()[sig_jet_nano_idx[0]], nano.Jet_eta()[sig_jet_nano_idx[0]], nano.Jet_phi()[sig_jet_nano_idx[0]], nano.Jet_mass()[sig_jet_nano_idx[0]]);
        j2.SetPtEtaPhiM(nano.Jet_pt()[sig_jet_nano_idx[1]], nano.Jet_eta()[sig_jet_nano_idx[1]], nano.Jet_phi()[sig_jet_nano_idx[1]], nano.Jet_mass()[sig_jet_nano_idx[1]]);
        dijet = j1 + j2;
        pico.out_llphoton_dijet_dphi().push_back(llg.DeltaPhi(dijet));
        pico.out_llphoton_dijet_balance().push_back((dilep+photon+j1+j2).Pt()/(dilep.Pt()+photon.Pt()+j1.Pt()+j2.Pt()));
        pico.out_llphoton_dijet_dr().push_back(llg.DeltaR(dijet));
        pico.out_photon_jet_mindr().push_back(min(photon.DeltaR(j1), photon.DeltaR(j2)));
        pico.out_photon_jet_maxdr().push_back(max(photon.DeltaR(j1), photon.DeltaR(j2)));
        pico.out_photon_jet1_dr().push_back(photon.DeltaR(j1));
        pico.out_photon_jet2_dr().push_back(photon.DeltaR(j2));
        pico.out_photon_zeppenfeld().push_back(abs(photon.Eta() - (j1.Eta() + j2.Eta())/2));

      }

      TVector3 g_pT = photon.Vect();
      TVector3 h_pT = llg.Vect();
      TVector3 z_pT = dilep.Vect();
      g_pT.SetZ(0); h_pT.SetZ(0); z_pT.SetZ(0);
      pico.out_llphoton_pTt().push_back( h_pT.Cross((z_pT-g_pT).Unit()).Mag() );
      pico.out_llphoton_pTt_an_hig019014().push_back( h_pT.Unit().Cross(z_pT-g_pT).Mag() );

      TLorentzVector lminus, lplus;
      TLorentzVector lep1, lep2;
      TLorentzVector l1err, l2err, pherr;
      double ptl1err, ptl2err, ptpherr;
      double dml1, dml2, dmph;
      ptpherr = pico.out_photon_energyErr()[igamma] * photon.Pt() / photon.P();
      pherr.SetPtEtaPhiM(pico.out_photon_pt()[igamma] + ptpherr,
                         pico.out_photon_eta()[igamma], pico.out_photon_phi()[igamma], 0);
      if(pico.out_ll_lepid()[ill] == 11) {
        int iel1 = pico.out_ll_i1()[ill];
        int iel2 = pico.out_ll_i2()[ill];
        if(pico.out_el_charge()[iel1] < 0) {
          lminus.SetPtEtaPhiM(pico.out_el_pt()[iel1], pico.out_el_eta()[iel1],
                              pico.out_el_phi()[iel1], 0.000511);
          lplus.SetPtEtaPhiM(pico.out_el_pt()[iel2], pico.out_el_eta()[iel2],
                             pico.out_el_phi()[iel2], 0.000511);
        }
        else {
          lplus.SetPtEtaPhiM(pico.out_el_pt()[iel1], pico.out_el_eta()[iel1],
                             pico.out_el_phi()[iel1], 0.000511);
          lminus.SetPtEtaPhiM(pico.out_el_pt()[iel2], pico.out_el_eta()[iel2],
                              pico.out_el_phi()[iel2], 0.000511);
        }
        lep1.SetPtEtaPhiM(pico.out_el_pt()[iel1], pico.out_el_eta()[iel1], pico.out_el_phi()[iel1], 0.000511);
        lep2.SetPtEtaPhiM(pico.out_el_pt()[iel2], pico.out_el_eta()[iel2], pico.out_el_phi()[iel2], 0.000511);
        ptl1err = pico.out_el_energyErr()[iel1] * lep1.Pt() / lep1.P();
        ptl2err = pico.out_el_energyErr()[iel2] * lep2.Pt() / lep2.P();
        l1err.SetPtEtaPhiM(pico.out_el_pt()[iel1] + ptl1err,
                           pico.out_el_eta()[iel1], pico.out_el_phi()[iel1], 0.000511);
        l2err.SetPtEtaPhiM(pico.out_el_pt()[iel2] + ptl2err,
                           pico.out_el_eta()[iel2], pico.out_el_phi()[iel2], 0.000511);
        dml1 = (l1err + lep2 + photon).M() - (lep1 + lep2 + photon).M();
        dml2 = (lep1 + l2err + photon).M() - (lep1 + lep2 + photon).M();
        dmph = (lep1 + lep2 + pherr).M() - (lep1 + lep2 + photon).M();
      }
      else {
        int imu1 = pico.out_ll_i1()[ill];
        int imu2 = pico.out_ll_i2()[ill];
        if(pico.out_mu_charge()[imu1] < 0) {
          lminus.SetPtEtaPhiM(pico.out_mu_pt()[imu1], pico.out_mu_eta()[imu1],
                              pico.out_mu_phi()[imu1], 0.10566);
          lplus.SetPtEtaPhiM(pico.out_mu_pt()[imu2], pico.out_mu_eta()[imu2],
                             pico.out_mu_phi()[imu2], 0.10566);
        }
        else {
          lplus.SetPtEtaPhiM(pico.out_mu_pt()[imu1], pico.out_mu_eta()[imu1],
                             pico.out_mu_phi()[imu1], 0.10566);
          lminus.SetPtEtaPhiM(pico.out_mu_pt()[imu2], pico.out_mu_eta()[imu2],
                              pico.out_mu_phi()[imu2], 0.10566);
        }
        lep1.SetPtEtaPhiM(pico.out_mu_pt()[imu1], pico.out_mu_eta()[imu1], pico.out_mu_phi()[imu1], 0.10566);
        lep2.SetPtEtaPhiM(pico.out_mu_pt()[imu2], pico.out_mu_eta()[imu2], pico.out_mu_phi()[imu2], 0.10566);
        ptl1err = pico.out_mu_ptErr()[imu1];
        ptl2err = pico.out_mu_ptErr()[imu2];
        l1err.SetPtEtaPhiM(pico.out_mu_pt()[imu1] + ptl1err,
                           pico.out_mu_eta()[imu1], pico.out_mu_phi()[imu1], 0.10566);
        l2err.SetPtEtaPhiM(pico.out_mu_pt()[imu2] + ptl2err,
                           pico.out_mu_eta()[imu2], pico.out_mu_phi()[imu2], 0.10566);
        dml1 = (l1err + lep2 + photon).M() - (lep1 + lep2 + photon).M();
        dml2 = (lep1 + l2err + photon).M() - (lep1 + lep2 + photon).M();
        dmph = (lep1 + lep2 + pherr).M() - (lep1 + lep2 + photon).M();
      }
      pico.out_llphoton_l1_masserr().push_back(dml1);
      pico.out_llphoton_l2_masserr().push_back(dml2);
      pico.out_llphoton_ph_masserr().push_back(dmph);

      // Variables used for defining kinematic angles presented in https://arxiv.org/pdf/1108.2274.pdf
      double M = llg.M(), mll = dilep.M();
      double lZ = sqrt(pow(llg.Dot(dilep)/M,2)-pow(mll,2));
      TVector3 hBoost = llg.BoostVector();

      // Cosine of angle between lepton 1 and parent Z
      double costheta = llg.Dot(lminus-lplus)/(M*lZ);
      pico.out_llphoton_costheta().push_back(costheta);

      // 4-momenta of q1/q2 (quarks from gluon-gluon fusion)
      TLorentzVector q, qBar;
      TVector3 hTransverseBoost = llg.BoostVector();
      hTransverseBoost.SetZ(0);
      TLorentzVector hH = llg;
      hH.Boost(-1*hTransverseBoost);
      double hPz = hH.Pz(), hE = hH.E();
      q.SetPxPyPzE(0, 0, (hPz+hE)/2, (hE+hPz)/2);
      qBar.SetPxPyPzE(0, 0, (hPz-hE)/2, (hE-hPz)/2);
      q.Boost(hTransverseBoost);
      qBar.Boost(hTransverseBoost);

      // Cosine of angle between incoming quarks and outgoing Zs in Higgs frame
      double cosTheta = (qBar-q).Dot(dilep)/(M*lZ);
      double sinTheta = sqrt(1 - pow(cosTheta, 2));
      pico.out_llphoton_cosTheta().push_back(cosTheta);

      // Angle phi
      dilep.Boost(-1*hBoost);
      TVector3 zBoost = dilep.BoostVector();
      q.Boost(-1*hBoost);
      lminus.Boost(-1*hBoost);
      lplus.Boost(-1*hBoost);
      TVector3 l1 = lminus.Vect(), l2 = lplus.Vect(), Z = dilep.Vect();
      TVector3 qvec = q.Vect();
      double cospsi = -1*l1.Cross(l2).Dot(qvec.Cross(Z))/l1.Cross(l2).Mag()/qvec.Cross(Z).Mag();
      double sinpsi = -1*l1.Cross(l2).Dot(qvec)/l1.Cross(l2).Mag()/qvec.Mag()/sinTheta;
      if (cospsi > 1) cospsi = 1;
      else if (cospsi < -1) cospsi = -1;
      double psi(0);
      if(sinpsi < 0) psi = -1*acos(cospsi);
      else           psi = acos(cospsi);
      pico.out_llphoton_psi().push_back(psi);
      // // Cosine of angle between incoming quarks and outgoing Zs in Higgs frame (Alternate formulation)
      // double cosThetaT = cos(dilep.Angle(q.Vect()));

      // // Angle phi (Alternate formulation)
      // TVector3 zAxis = dilep.Vect().Unit();
      // TVector3 yAxis = qvec.Cross(zAxis.Unit()).Unit();
      // TVector3 xAxis = (yAxis.Unit().Cross(zAxis.Unit())).Unit();
      // TRotation rotation;
      // rotation = rotation.RotateAxes(xAxis,yAxis,zAxis).Inverse();
      // lminus.Transform(rotation);
      // double psiT = lminus.Phi();

      // // Cosine of angle between lepton 1 and parent Z (Alternate formulation)
      // lminus.Boost(-1*zBoost);
      // lplus.Boost(-1*zBoost);
      // double costhetaT = cos(dilep.Angle(lplus.Vect()));

      // cout << "Phi: " << psi << " " << psiT << endl;
      // cout << "costheta: " << costheta << " " << costhetaT << endl;
      // cout << "cosTheta: " << cosTheta << " " << cosThetaT << endl;

      //pico.out_llphoton_costhj().push_back(cosThetaJeff(lminus,lplus,photon));
    }
  }

  //Code Executing the kinematic refit
  if(pico.out_nllphoton() > 0){

    //Definitions of variables for the refit
    TLorentzVector ll_refit(0,0,0,0);
    std::vector<TLorentzVector> refit_leptons{ll_refit,ll_refit};
    std::map<unsigned int, TLorentzVector> leptons_map;
    std::map<unsigned int, double> leptons_pterr_map;
    std::map<unsigned int, TLorentzVector> fsrphotons_map;
    TLorentzVector fsrphoton1, fsrphoton2,l1,l2,photon;
    int status, covmatstatus = -5;
    float minnll = -5;
    int idx_l1 = pico.out_ll_i1()[0];
    int idx_l2 = pico.out_ll_i2()[0];

    if(pico.out_ll_lepid()[0] == 13){
       int idx_fsr1, idx_fsr2;

       l1.SetPtEtaPhiM(pico.out_mu_pt()[idx_l1], pico.out_mu_eta()[idx_l1],
                       pico.out_mu_phi()[idx_l1], 0.10566);
       l2.SetPtEtaPhiM(pico.out_mu_pt()[idx_l2], pico.out_mu_eta()[idx_l2],
                       pico.out_mu_phi()[idx_l2], 0.10566);

      leptons_map[0] = l1;
      leptons_map[1] = l2;
      leptons_pterr_map[0] = pico.out_mu_ptErr()[idx_l1];
      leptons_pterr_map[1] = pico.out_mu_ptErr()[idx_l2];
     
      //Loops through FSRphotons to find which leptons are associated with to be used for constrained fit
      idx_fsr1 = -1; idx_fsr2 = -1;
      for(size_t idx_fsr = 0; idx_fsr < static_cast<unsigned int>(pico.out_nfsrphoton()); idx_fsr++){
        if(pico.out_fsrphoton_muonidx()[idx_fsr]==idx_l1){
          if(idx_fsr1 != -1 && (pico.out_fsrphoton_droveret2()[idx_fsr1] <  pico.out_fsrphoton_droveret2()[idx_fsr])){ continue;}
            idx_fsr1 = idx_fsr; continue;
        }
        if(pico.out_fsrphoton_muonidx()[idx_fsr]==idx_l2){
          if(idx_fsr2 != -1 && (pico.out_fsrphoton_droveret2()[idx_fsr2] <  pico.out_fsrphoton_droveret2()[idx_fsr])){ continue;}
            idx_fsr2 = idx_fsr; continue;
          }
        }

      //Adds the FSR photons to a map
      fsrphotons_map.clear();
      int fsr_cnt = 0;
      if(idx_fsr1!=-1){
        fsrphoton1.SetPtEtaPhiM( pico.out_fsrphoton_pt()[idx_fsr1],pico.out_fsrphoton_eta()[idx_fsr1],pico.out_fsrphoton_phi()[idx_fsr1],0 );
        fsrphotons_map[fsr_cnt] = fsrphoton1;
        fsr_cnt++;
      }
      if(idx_fsr2!=-1){
        fsrphoton2.SetPtEtaPhiM( pico.out_fsrphoton_pt()[idx_fsr2],pico.out_fsrphoton_eta()[idx_fsr2],pico.out_fsrphoton_phi()[idx_fsr2],0 ); 
        fsrphotons_map[fsr_cnt] = fsrphoton2;     
      } 
    } else {
      l1.SetPtEtaPhiM(pico.out_el_pt()[idx_l1], pico.out_el_eta()[idx_l1],
                      pico.out_el_phi()[idx_l1], 0.000511);
      l2.SetPtEtaPhiM(pico.out_el_pt()[idx_l2], pico.out_el_eta()[idx_l2],
                      pico.out_el_phi()[idx_l2], 0.000511);

      leptons_map[0] = l1;
      leptons_map[1] = l2; 
      leptons_pterr_map[0] = pico.out_el_energyErr()[idx_l1]*l1.Pt()/l1.P();
      leptons_pterr_map[1] = pico.out_el_energyErr()[idx_l2]*l2.Pt()/l2.P();
    }

      kinZfitter->Setup(leptons_map, fsrphotons_map, leptons_pterr_map);
      kinZfitter->KinRefitZ1();
      refit_leptons = kinZfitter->GetRefitP4s();
      ll_refit     = refit_leptons[0] + refit_leptons[1];
      status       = kinZfitter -> GetStatus();
      covmatstatus = kinZfitter -> GetCovMatStatus();
      minnll       = kinZfitter -> GetMinNll();

      //debug statement
      //cnt_refit++; std::cout << cnt_refit << std::endl;

      photon.SetPtEtaPhiM(pico.out_photon_pt()[0], pico.out_photon_eta()[0], pico.out_photon_phi()[0], 0);
      TLorentzVector llg_refit = ll_refit + photon;

      //Assigns the refit lepton/dilepton quantities
      pico.out_ll_refit_pt()    = ll_refit.Pt();
      pico.out_ll_refit_eta()   = ll_refit.Eta();
      pico.out_ll_refit_phi()   = ll_refit.Phi();
      pico.out_ll_refit_m()     = ll_refit.M();
      pico.out_ll_refit_l1_pt() = refit_leptons[0].Pt();
      pico.out_ll_refit_l2_pt() = refit_leptons[1].Pt();

      pico.out_ll_refit_status()        = status;
      pico.out_ll_refit_covmat_status() = covmatstatus;
      pico.out_ll_refit_minnll()        = minnll;

      pico.out_llphoton_refit_pt()   = llg_refit.Pt();
      pico.out_llphoton_refit_eta()  = llg_refit.Eta();
      pico.out_llphoton_refit_phi()  = llg_refit.Phi();
      pico.out_llphoton_refit_m()    = llg_refit.M();
      pico.out_llphoton_refit_dr()   = ll_refit.DeltaR(photon);
      pico.out_llphoton_refit_dphi() = ll_refit.DeltaPhi(photon);
      pico.out_llphoton_refit_deta() = fabs(ll_refit.Eta() - photon.Eta());

      if(sig_jet_nano_idx.size() > 1) {
        TLorentzVector j1, j2, dijet;
        j1.SetPtEtaPhiM(nano.Jet_pt()[sig_jet_nano_idx[0]], nano.Jet_eta()[sig_jet_nano_idx[0]], nano.Jet_phi()[sig_jet_nano_idx[0]], nano.Jet_mass()[sig_jet_nano_idx[0]]);
        j2.SetPtEtaPhiM(nano.Jet_pt()[sig_jet_nano_idx[1]], nano.Jet_eta()[sig_jet_nano_idx[1]], nano.Jet_phi()[sig_jet_nano_idx[1]], nano.Jet_mass()[sig_jet_nano_idx[1]]);
        dijet = j1 + j2;
        pico.out_llphoton_refit_dijet_dphi()    = llg_refit.DeltaPhi(dijet);
        pico.out_llphoton_refit_dijet_balance() = (ll_refit+photon+j1+j2).Pt()/(ll_refit.Pt()+photon.Pt()+j1.Pt()+j2.Pt());
        pico.out_llphoton_refit_dijet_dr()      = llg_refit.DeltaR(dijet);

      }


      //Refit angles and various variables using the kinematic refit quantities
      TVector3 g_pT = photon.Vect();
      TVector3 h_pT = llg_refit.Vect();
      TVector3 z_pT = ll_refit.Vect();
      g_pT.SetZ(0); h_pT.SetZ(0); z_pT.SetZ(0);
      pico.out_llphoton_refit_pTt() = h_pT.Cross((z_pT-g_pT).Unit()).Mag();
      pico.out_llphoton_refit_pTt_an_hig019014() = h_pT.Unit().Cross(z_pT-g_pT).Mag();

      TLorentzVector lminus, lplus;
      TLorentzVector lep1, lep2;
      TLorentzVector l1err, l2err, pherr;
      double ptl1err, ptl2err, ptpherr;
      double dml1, dml2, dmph;
      ptpherr = pico.out_photon_energyErr()[0] * photon.Pt() / photon.P();
      pherr.SetPtEtaPhiM(pico.out_photon_pt()[0] + ptpherr,
                         pico.out_photon_eta()[0], pico.out_photon_phi()[0], 0);
      lep1 = refit_leptons[0];
      lep2 = refit_leptons[1];
      if(pico.out_ll_lepid()[0] == 11) {
        int iel1 = pico.out_ll_i1()[0];
        int iel2 = pico.out_ll_i2()[0];
        if(pico.out_el_charge()[iel1] < 0) {
          lminus = lep1;
          lplus  = lep2;
        }
        else {
          lminus = lep2;
          lplus  = lep1;
        }
        ptl1err = pico.out_el_energyErr()[iel1] * lep1.Pt() / lep1.P();
        ptl2err = pico.out_el_energyErr()[iel2] * lep2.Pt() / lep2.P();
        l1err.SetPtEtaPhiM(lep1.Pt() + ptl1err, lep1.Eta(), lep1.Phi(), 0.000511);
        l2err.SetPtEtaPhiM(lep2.Pt() + ptl2err, lep2.Eta(), lep2.Phi(), 0.000511);
        dml1 = (l1err + lep2 + photon).M() - (lep1 + lep2 + photon).M();
        dml2 = (lep1 + l2err + photon).M() - (lep1 + lep2 + photon).M();
        dmph = (lep1 + lep2 + pherr).M() - (lep1 + lep2 + photon).M();
      }
      else {
        int imu1 = pico.out_ll_i1()[0];
        int imu2 = pico.out_ll_i2()[0];
        if(pico.out_mu_charge()[imu1] < 0) {
          lminus = lep1;
          lplus  = lep2;
        }
        else {
          lminus = lep2;
          lplus  = lep1;
        }
        ptl1err = pico.out_mu_ptErr()[imu1];
        ptl2err = pico.out_mu_ptErr()[imu2];
        l1err.SetPtEtaPhiM(lep1.Pt() + ptl1err, lep1.Eta(), lep1.Phi(), 0.1056);
        l2err.SetPtEtaPhiM(lep2.Pt() + ptl2err, lep2.Eta(), lep2.Phi(), 0.1056);
        dml1 = (l1err + lep2 + photon).M() - (lep1 + lep2 + photon).M();
        dml2 = (lep1 + l2err + photon).M() - (lep1 + lep2 + photon).M();
        dmph = (lep1 + lep2 + pherr).M() - (lep1 + lep2 + photon).M();
      }
      pico.out_llphoton_refit_l1_masserr() = dml1;
      pico.out_llphoton_refit_l2_masserr() = dml2;
      pico.out_llphoton_refit_ph_masserr() = dmph;

      // Variables used for defining kinematic angles presented in https://arxiv.org/pdf/1108.2274.pdf
      double M = llg_refit.M(), mll_refit = ll_refit.M();
      double lZ = sqrt(pow(llg_refit.Dot(ll_refit)/M,2)-pow(mll_refit,2));
      TVector3 hBoost = llg_refit.BoostVector();

      // Cosine of angle between lepton 1 and parent Z
      double costheta = llg_refit.Dot(lminus-lplus)/(M*lZ);
      pico.out_llphoton_refit_costheta() = costheta;

      // 4-momenta of q1/q2 (quarks from gluon-gluon fusion)
      TLorentzVector q, qBar;
      TVector3 hTransverseBoost = llg_refit.BoostVector();
      hTransverseBoost.SetZ(0);
      TLorentzVector hH = llg_refit;
      hH.Boost(-1*hTransverseBoost);
      double hPz = hH.Pz(), hE = hH.E();
      q.SetPxPyPzE(0, 0, (hPz+hE)/2, (hE+hPz)/2);
      qBar.SetPxPyPzE(0, 0, (hPz-hE)/2, (hE-hPz)/2);
      q.Boost(hTransverseBoost);
      qBar.Boost(hTransverseBoost);

      // Cosine of angle between incoming quarks and outgoing Zs in Higgs frame
      double cosTheta = (qBar-q).Dot(ll_refit)/(M*lZ);
      double sinTheta = sqrt(1 - pow(cosTheta, 2));
      pico.out_llphoton_refit_cosTheta() = cosTheta;

      // Angle phi
      ll_refit.Boost(-1*hBoost);
      TVector3 zBoost = ll_refit.BoostVector();
      q.Boost(-1*hBoost);
      lminus.Boost(-1*hBoost);
      lplus.Boost(-1*hBoost);
      TVector3 l1_vec = lminus.Vect(), l2_vec = lplus.Vect(), Z_vec = ll_refit.Vect();
      TVector3 qvec = q.Vect();
      double cospsi = -1*l1_vec.Cross(l2_vec).Dot(qvec.Cross(Z_vec))/l1_vec.Cross(l2_vec).Mag()/qvec.Cross(Z_vec).Mag();
      double sinpsi = -1*l1_vec.Cross(l2_vec).Dot(qvec)/l1_vec.Cross(l2_vec).Mag()/qvec.Mag()/sinTheta;
      if (cospsi > 1) cospsi = 1;
      else if (cospsi < -1) cospsi = -1;
      double psi(0);
      if(sinpsi < 0) psi = -1*acos(cospsi);
      else           psi = acos(cospsi);
      pico.out_llphoton_refit_psi() = psi;


    }

  //================Bitmap for zgamma cut flow================//


  //Read from left. Bit 1 n_ll>=2, Bit 2 n_ee>=2, Bit 3 n_mumu>=2, Bit 4 trigs single lep, Bit 5 trigs dilep, Bit 6 leading lepton pT, Bit 7 subleading lepton pT
  // Bit 8 nphoton>=1, Bit 9 photon_id80, Bit 10 m_ll cut, Bit 11 m_lly cut, Bit 12 15/110, Bit 13 is 1 if outside m_lly signal region
  //New:
  //Read from left. Bit 1 n_ll>=2, Bit 2 n_ee>=2, Bit 3 n_mumu>=2, Bit 4 passes triggers Bit 5 passes lep pT cut, Bit 6 nphoton>=1, Bit 7 m_ll cut, Bit 8 15/110 ratio cut, bit 9 m_lly+mll cut, Bit 10 m_lly cut, Bit 11 pass filters. Bit 12 signal blinding window
  if(pico.out_nel()>=2 && pico.out_ll_lepid().at(0)==11){
    baseBit += 0b110000000000; 
    if(pico.out_trig_single_el() || pico.out_trig_double_el()){baseBit +=0b000100000000;}
    if(pico.out_trig_el_pt()){baseBit += 0b000010000000;}
  }
  if(pico.out_nmu()>=2 && pico.out_ll_lepid().at(0)==13){
    baseBit += 0b101000000000;
    if(pico.out_trig_single_mu() || pico.out_trig_double_mu()){baseBit +=0b000100000000;}
    if(pico.out_trig_mu_pt()){baseBit +=0b000010000000;}
  }
  if(pico.out_nphoton()>=1){baseBit+= 0b000001000000;}
  if(pico.out_ll_m().at(pico.out_llphoton_ill().at(0))>=80.f && pico.out_ll_m().at(pico.out_llphoton_ill().at(0))<=100.f){baseBit+= 0b000000100000;}
  if(pico.out_photon_pt().at(pico.out_llphoton_iph().at(0))/pico.out_llphoton_m().at(0) >=15.0f/110.f){baseBit+= 0b000000010000;}
  if(pico.out_ll_m().at(pico.out_llphoton_ill().at(0))+pico.out_llphoton_m().at(0) > 185.f){baseBit+= 0b000000001000;}
  if(pico.out_llphoton_m().at(0)>=100.f && pico.out_llphoton_m().at(0)<=180.f){baseBit+= 0b000000000100;}
  if(pico.out_pass()){baseBit+=0b000000000010;}
  if(pico.out_llphoton_m().at(0)<=120.f || pico.out_llphoton_m().at(0)>=130.f){baseBit+= 0b000000000001;}
  pico.out_zg_cutBitMap() = baseBit;
  //End Bitmap for zgamma cut flow
            
  //Category bitmap - checks if the event matches one of the categories and their baselines
  //Bit 0: ggF, bit 1: VBF, bit 2: ttH leptonic, bit 3: VH 3l, bit 4: ttH hadronic, bit 5: ZH ptmiss, bit 6: untagged, bit 7: category specific baseline selection
  //ggF
  if(pico.out_nlep()==2 && pico.out_njet()<=1 && pico.out_met()<90.f){ categoryBit += 0b10000000; } 
  if(pico.out_nlep()==2 && pico.out_njet()>=2 && pico.out_nbdfm()==0){ categoryBit += 0b01000000; } //VBF
  
  //ttH leptonic
  if((pico.out_nlep()==3 && pico.out_njet()>=3 && pico.out_nbdfm()>=1) || (pico.out_nlep()>=4 && pico.out_nbdfm() >= 1)){
    categoryBit+=0b00100000; 

    float mll = pico.out_ll_m().at(pico.out_llphoton_ill().at(0));
    bool pass_miniso = check_miniso(pico,0.1);
    //Category selections
    if(pass_miniso && mll > 85.f && mll < 95.f){categoryBit+=0b00000001;}
  }

  //VH 3l
  if(pico.out_nlep()>=3 && pico.out_nbdfm()==0){
    categoryBit += 0b00010000;

    //Category selections
    float ptom_llgamma = pico.out_llphoton_pt().at(0)/pico.out_llphoton_m().at(0);
    float mll = pico.out_ll_m().at(pico.out_llphoton_ill().at(0));
    bool pass_miniso = check_miniso(pico,0.15);
    if(pass_miniso && pico.out_met() > 30.0f && ptom_llgamma > 0.3f && mll > 85.f && mll < 95.f){categoryBit+=0b00000001;}
  }

  //ttH hadronic
  if(pico.out_nlep()==2 && pico.out_njet()>=5 && pico.out_nbdfm()>=1){ 
    categoryBit += 0b00001000;

    //Category selections
    float mll = pico.out_ll_m().at(pico.out_llphoton_ill().at(0));
    if(mll > 85.f && mll < 95.f){categoryBit+=0b00000001;}
  }

  //ZH ptmiss
  if(pico.out_nlep()==2 && pico.out_njet()<=1 && pico.out_met()>90.f){ 
    categoryBit += 0b00000100;

    //Category selections
    float mll = pico.out_ll_m().at(pico.out_llphoton_ill().at(0));
    float ptom_llgamma = pico.out_llphoton_pt().at(0)/pico.out_llphoton_m().at(0);
    if(mll > 85.f && mll < 95.f && ptom_llgamma > 0.4f){categoryBit+=0b00000001;}
  }

  //"Untagged" aka events that pass baseline but arent put into a category
  if(categoryBit==0b00000000){categoryBit = 0b00000010;};

  //set categorization BitMap in picos
  pico.out_zg_categorizationBitMap() = categoryBit;
  //End of category bitmap

  return;
}

