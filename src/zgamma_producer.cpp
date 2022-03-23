#include "zgamma_producer.hpp"

#include "utilities.hpp"
#include "TLorentzVector.h"

using namespace std;

ZGammaVarProducer::ZGammaVarProducer(int year_){
    year = year_;
}

ZGammaVarProducer::~ZGammaVarProducer(){
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
  pico.out_nllphoton() = 0;
  if (pico.out_ll_pt().size() == 0 || pico.out_nphoton() == 0) return;
  for(size_t ill(0); ill < pico.out_ll_pt().size(); ill++)
    for(size_t igamma(0); igamma < pico.out_photon_pt().size(); igamma++) {
      if(pico.out_photon_sig()[igamma]) {
        TLorentzVector dilep, photon, llg;
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
	  pico.out_photon_jet_mindr().push_back(min(photon.DeltaR(j1), photon.DeltaR(j2)));
	  pico.out_photon_zeppenfeld().push_back(abs(photon.Eta() - (j1.Eta() + j2.Eta())/2));
	  TVector3 g = photon.Vect();
	  TVector3 h = llg.Vect();
	  TVector3 z = dilep.Vect();
	  g.SetZ(0); h.SetZ(0); z.SetZ(0);
	  pico.out_llphoton_pTt2().push_back(h.Cross(z-g).Mag()/h.Mag());
	}
        TLorentzVector lminus, lplus;
        TLorentzVector lep1, lep2;
        TLorentzVector l1err, l2err, pherr;
        double ptl1err, ptl2err, ptpherr;
        double dml1, dml2, dmph;
        ptpherr = pico.out_photon_pterr()[igamma] * photon.Pt() / photon.P();
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

        pico.out_llphoton_costhj().push_back(cosThetaJeff(lminus,lplus,photon));
      }
    }
  return;
}

