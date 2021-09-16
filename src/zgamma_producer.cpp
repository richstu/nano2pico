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
    j1.SetPtEtaPhiM(nano.Jet_pt()[sig_jet_nano_idx[0]], nano.Jet_eta()[sig_jet_nano_idx[0]], 
                    nano.Jet_phi()[sig_jet_nano_idx[0]],nano.Jet_mass()[sig_jet_nano_idx[0]]);
    j2.SetPtEtaPhiM(nano.Jet_pt()[sig_jet_nano_idx[1]], nano.Jet_eta()[sig_jet_nano_idx[1]], 
                    nano.Jet_phi()[sig_jet_nano_idx[1]],nano.Jet_mass()[sig_jet_nano_idx[1]]);
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
        photon.SetPtEtaPhiM(pico.out_photon_pt()[igamma], pico.out_photon_eta()[igamma],
                            pico.out_photon_phi()[igamma],0);
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
        TLorentzVector lminus, lplus;
        if(pico.out_ll_lepid()[ill] == 11) {
          int iel1 = pico.out_ll_i1()[ill];
          int iel2 = pico.out_ll_i2()[ill];
          if(pico.out_el_charge()[iel1] < 0) {
            lminus.SetPtEtaPhiM(pico.out_el_pt()[iel1], pico.out_el_eta()[iel1],
                                pico.out_el_phi()[iel1],0.0005);
            lplus.SetPtEtaPhiM(pico.out_el_pt()[iel2], pico.out_el_eta()[iel2],
                               pico.out_el_phi()[iel2],0.0005);
          }
          else {
            lplus.SetPtEtaPhiM(pico.out_el_pt()[iel1], pico.out_el_eta()[iel1],
                               pico.out_el_phi()[iel1],0.0005);
            lminus.SetPtEtaPhiM(pico.out_el_pt()[iel2], pico.out_el_eta()[iel2],
                                pico.out_el_phi()[iel2],0.0005);
          }
        }
        else {
          int imu1 = pico.out_ll_i1()[ill];
          int imu2 = pico.out_ll_i2()[ill];
          if(pico.out_mu_charge()[imu1] < 0) {
            lminus.SetPtEtaPhiM(pico.out_mu_pt()[imu1], pico.out_mu_eta()[imu1],
                                pico.out_mu_phi()[imu1],0.105);
            lplus.SetPtEtaPhiM(pico.out_mu_pt()[imu2], pico.out_mu_eta()[imu2],
                               pico.out_mu_phi()[imu2],0.105);
          }
          else {
            lplus.SetPtEtaPhiM(pico.out_mu_pt()[imu1], pico.out_mu_eta()[imu1],
                               pico.out_mu_phi()[imu1],0.105);
            lminus.SetPtEtaPhiM(pico.out_mu_pt()[imu2], pico.out_mu_eta()[imu2],
                                pico.out_mu_phi()[imu2],0.105);
          }
        }

        // Variables used for defining kinematic angles presented in https://arxiv.org/pdf/1108.2274.pdf
        TVector3 l1 = lminus.Vect(), l2 = lplus.Vect(), Z = dilep.Vect();
        double M = llg.M(), mll = dilep.M();
        double lZ = sqrt(pow(llg.Dot(dilep)/M,2)-pow(mll,2));
        TLorentzVector q1, q2;
        TVector3 hTransverseBoost = llg.BoostVector();
        hTransverseBoost.SetZ(0);
        TLorentzVector hZ = llg;
        hZ.Boost(-1*hTransverseBoost);
        double hPz = hZ.Pz(), hE = hZ.E();

        // 4-momenta of q1/q2 (quarks from gluon-gluon fusion)
        //  Defined in Equation 4/5
        q1.SetPxPyPzE(0,0,(hPz+hE)/2,(hE+hPz)/2);
        q2.SetPxPyPzE(0,0,(hPz-hE)/2,(hE-hPz)/2);
        q1.Boost(hTransverseBoost);
        q2.Boost(hTransverseBoost);
        TVector3 q1vec = q1.Vect();

        // Cosine of angle between lepton 1 and parent Z in Higgs frame 
        //  Defined in Equation 13
        double costheta = llg.Dot(lminus-lplus)/(M*lZ);
        pico.out_llphoton_costheta().push_back(costheta);

        // Cosine of angle between incoming quarks and outgoing Zs in higgs frame 
        //  Defined in Equation 8
        double cosTheta = dilep.Dot(q1-q2)/(M*lZ);
        pico.out_llphoton_cosTheta().push_back(cosTheta);

        // Angle of the Z decay plane from the z-axis (defined in Equation 1) in the higgs frame
        //  Defined in Equation 21+22 (it's called phi in the paper, but psi is
        //  used here to distinguish it from the detector angle phi)
        double cospsi = -1*l1.Cross(l2).Dot(q1vec.Cross(Z))/l1.Cross(l2).Mag()/q1vec.Cross(Z).Mag();
        double sinpsi = -1*l1.Cross(l2).Dot(q1vec)/l1.Cross(l2).Mag()/q1vec.Mag();
        double psi(0);
        if(cospsi > 1) cospsi = 1;
        if(cospsi < -1) cospsi = -1;
        if(sinpsi < 0) psi = -1*acos(cospsi);
        else           psi = acos(cospsi);
        pico.out_llphoton_psi().push_back(psi);

        pico.out_llphoton_costhj().push_back(cosThetaJeff(lminus,lplus,photon));
      }
    }
  return;
}

