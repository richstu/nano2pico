/*************************************************************************
 *  Authors:   Tongguang Cheng
 *  Adapted for H-->Zgamma by: Prasanna Siddireddy and Anders Barzdukas 
 *  Contact at abarzdukas@ucsb.edu
 *************************************************************************/
#ifndef KinZfitter_cpp
#define KinZfitter_cpp

/// KinFitter header
#include "KinZfitter.hpp"

KinZfitter::KinZfitter() {

  //Default values drawn from HZg_Crystal_ball_and_3Gaussian_fit.txt
  PDFName_ = "./txt/constrained_fit_input/HZg_Crystal_ball_and_3Gaussian_fit.txt";
  meanCB_      = 90.8919;
  sigmaCB_     = 4.007;
  alphaCB_     = 1.1981;
  nCB_         = 3.25604;
  meanGauss1_  = 96.4278;
  sigmaGauss1_ = 6.17509;
  f1_          = 0.86449;
  meanGauss2_  = 91.1649;
  sigmaGauss2_ = 0.856305;
  f2_          = 0.514371;
  meanGauss3_  = 91.1513;
  sigmaGauss3_ = 1.77148;
  f3_          = 0.648225;
  threegauss_  = true;

  //This flag will bypass the section that reads a text file for the input values
  debug_ = false;
  if(debug_) std::cout << "KinZfitter. The debug flag is ON with "<<PDFName_<< std::endl;
}

//This constructor takes a .txt file with the probability distribution function
KinZfitter::KinZfitter(TString pdf_txtfile){
  PDFName_ = pdf_txtfile;
  debug_ = false;

  std::ifstream input(PDFName_);
  std::string line;
  if(debug_) cout<<"PDFName_ in "<<PDFName_<<endl;

  if(PDFName_.Contains("3G")){
    threegauss_ = true;
    while (!input.eof() && std::getline(input,line))
      {
        std::istringstream iss(line);
        string p; double val;
        if(iss >> p >> val) {
          if(p=="meanCB")      { meanCB_ = val;}
          if(p=="sigmaCB")     { sigmaCB_ = val;}
          if(p=="alphaCB")     { alphaCB_ = val;}
          if(p=="nCB")         { nCB_ = val;}
          if(p=="meanGauss1")  { meanGauss1_ = val;}
          if(p=="sigmaGauss1") { sigmaGauss1_ = val;}
          if(p=="f1")          { f1_ = val;}
          if(p=="meanGauss2")  { meanGauss2_ = val; }
          if(p=="sigmaGauss2") { sigmaGauss2_ = val;}
          if(p=="f2")          { f2_ = val;}
          if(p=="meanGauss3")  { meanGauss3_ = val;}
          if(p=="sigmaGauss3") { sigmaGauss3_ = val;}
          if(p=="f3")          { f3_ = val;}
        }
      }
  }
  else {
    threegauss_ = false;
    while (!input.eof() && std::getline(input,line))
      {
      std::istringstream iss(line);
      string p; double val;
      if(iss >> p >> val) {
        if(p=="bwMean")  { BWmean_ = val; }
        if(p=="bwGamma" ){ BWgamma_ = val;  }
        if(p=="Gsigma" ) { sigmaValG_ = val; }
      }
    }
  }
  input.close(); 

  if(debug_) std::cout << "KinZfitter. The debug flag is ON with "<<PDFName_<< std::endl;
}


double KinZfitter::masserrorFullCov(std::vector<TLorentzVector> p4s, TMatrixDSym covMatrix){

  int ndim = 3*p4s.size();
  if(debug_) cout<<""<<endl;
  TMatrixD jacobian(1,ndim);
  double e = 0; double mass = 0;
  double px = 0; double py = 0; double pz = 0;
  for (unsigned int ip = 0; ip < p4s.size(); ip++)
    {
      e = e + p4s[ip].E();
      px = px + p4s[ip].Px();
      py = py + p4s[ip].Py();
      pz = pz + p4s[ip].Pz();
    }
  mass = TMath::Sqrt(e*e-px*px-py*py-pz*pz);
  for (unsigned int i = 0, o = 0; i < p4s.size(); i++, o += 3)
    {
      double pxi = p4s[i].Px();
      double pyi = p4s[i].Py();
      double pzi = p4s[i].Pz();
      double ei = p4s[i].E();
      jacobian(0, o+0) = (e*(pxi/ei) - px)/mass;
      jacobian(0, o+1) = (e*(pyi/ei) - py)/mass;
      jacobian(0, o+2) = (e*(pzi/ei) - pz)/mass;
    }
  TMatrixDSym massCov = covMatrix.Similarity(jacobian);
  double dm2 = massCov(0,0);
  return (dm2 > 0 ? std::sqrt(dm2) : 0.0);

}


double KinZfitter::masserror( std::vector<TLorentzVector> Lep, std::vector<double> pterr){

  TLorentzVector compositeParticle ;
  for(unsigned int i=0; i<Lep.size(); i++)
    {
      compositeParticle+=Lep[i];
    }
  double mass  =  compositeParticle.M();
  double masserr = 0;
  for(unsigned int i=0; i<Lep.size(); i++)
    {
      TLorentzVector variedLep; // = Lep[i];
      variedLep.SetPtEtaPhiM(Lep[i].Pt()+ pterr[i], Lep[i].Eta(), Lep[i].Phi(), Lep[i].M());
      TLorentzVector compositeParticleVariation ;
      for(unsigned int j=0; j<Lep.size(); j++)
        {
          if(i!=j)compositeParticleVariation+=Lep[j];
          else compositeParticleVariation+=variedLep;
        }
     masserr += (compositeParticleVariation.M()-mass)*(compositeParticleVariation.M()-mass);
    }
  return sqrt(masserr);

}


double KinZfitter::pterr(TLorentzVector ph){

  double C, S, N;
  if (abs(ph.Eta()) < 1.48) {
    C = 0.35 / 100;
    S = 5.51 / 100;
    N = 98. / 1000.;
  } else {
    C = 0;
    S = 12.8 / 100;
    N = 440. / 1000.;
  }
  double pherr = sqrt(C * C * ph.Energy() * ph.Energy() + S * S * ph.Energy() + N * N);
  return pherr;

}


void KinZfitter::Setup(std::map<unsigned int, TLorentzVector> selectedLeptons, std::map<unsigned int, TLorentzVector> selectedFsrPhotons, std::map<unsigned int, double> errorLeptons) {

  // reset everything for each event
  p4sZ1_.clear();
  p4sZ1ph_.clear();
  p4sZ1REFIT_.clear();
  p4sZ1phREFIT_.clear();

  pTerrsZ1_.clear();
  pTerrsZ1ph_.clear();
  pTerrsZ1REFIT_.clear();
  pTerrsZ1phREFIT_.clear();

  gErrorIgnoreLevel = kWarning;
  RooMsgService::instance().setStreamStatus(1,false);
  initZs(selectedLeptons, selectedFsrPhotons, errorLeptons);
  if(debug_){ cout << "Setup complete" << endl;} 
}

///----------------------------------------------------------------------------------------------
void KinZfitter::initZs(std::map<unsigned int, TLorentzVector> selectedLeptons, std::map<unsigned int, TLorentzVector> selectedFsrPhotons, std::map<unsigned int, double> errorLeptons) {

  if(debug_) cout<<"init leptons"<<endl;

  for(unsigned int il = 0; il < selectedLeptons.size(); il++)
    {
      double pTerr = 0;

      TLorentzVector p4 = selectedLeptons[il];
      if(debug_) cout << "lep_error before corrections = " << errorLeptons[il] << endl; 

      pTerr = errorLeptons[il];
      if(debug_) cout<<" pt err is "<<pTerr<<endl;

      pTerrsZ1_.push_back(pTerr);
      p4sZ1_.push_back(p4);

    }

  if(debug_) cout<<"init fsr photons"<<endl;

    TLorentzVector p4;
    for(unsigned int ifsr = 0; ifsr < selectedFsrPhotons.size(); ifsr++)
      {
        p4 = selectedFsrPhotons[ifsr];
        if(selectedFsrPhotons[ifsr].Pt()==0) continue;

        if(debug_) cout<<"ifsr "<<ifsr<<endl;

        double pTerr = 0;

        pTerr = pterr(p4);

        if(debug_) cout<<" pt err is "<<pTerr<<endl;
        if(debug_) cout<<"for fsr Z1 photon"<<endl;

        pTerrsZ1ph_.push_back(pTerr);
        p4sZ1ph_.push_back(p4);
    }
//  }
  if(debug_) cout<<"p4sZ1ph_ "<<p4sZ1ph_.size()<<endl;

}

void KinZfitter::SetZ1Result(double l1, double l2, double lph1, double lph2) {

  if(debug_) cout<<"start set Z1 result"<<endl;

  // pT scale after refitting w.r.t. reco pT
  lZ1_l1_ = l1; lZ1_l2_ = l2;
  if(debug_) cout<<"l1 "<<l1<<" l2 "<<l2<<endl;
  lZ1_ph1_ = lph1; lZ1_ph2_ = lph2;

  TLorentzVector Z1_1 = p4sZ1_[0]; TLorentzVector Z1_2 = p4sZ1_[1];

  TLorentzVector Z1_1_True(0,0,0,0);
  Z1_1_True.SetPtEtaPhiM(lZ1_l1_*Z1_1.Pt(),Z1_1.Eta(),Z1_1.Phi(),Z1_1.M());
  TLorentzVector Z1_2_True(0,0,0,0);
  Z1_2_True.SetPtEtaPhiM(lZ1_l2_*Z1_2.Pt(),Z1_2.Eta(),Z1_2.Phi(),Z1_2.M());

  p4sZ1REFIT_.push_back(Z1_1_True); p4sZ1REFIT_.push_back(Z1_2_True);

  TLorentzVector Z1ph;
  TLorentzVector Z1phTrue(0,0,0,0);
  for(unsigned int ifsr1 = 0; ifsr1<p4sZ1ph_.size(); ifsr1++) {
    Z1ph = p4sZ1ph_[ifsr1];

    double l = 1.0;
    if(ifsr1==0){ l = lZ1_ph1_;}
    if(ifsr1==1){ l = lZ1_ph2_;}

    Z1phTrue.SetPtEtaPhiM(l*Z1ph.Pt(),Z1ph.Eta(),Z1ph.Phi(),Z1ph.M());

    p4sZ1phREFIT_.push_back(Z1phTrue);
  }

  if(debug_) cout<<"end set Z1 result"<<endl;

}


double KinZfitter::GetRefitMZ1()
{

  vector<TLorentzVector> p4s = GetRefitP4s();

  TLorentzVector pZ1(0,0,0,0);

  pZ1 = p4s[0] + p4s[1];

  return pZ1.M();

}


double KinZfitter::GetMZ1Err()
{

  vector<TLorentzVector> p4s;
  vector<double> pTErrs;

  p4s.push_back(p4sZ1_[0]);
  p4s.push_back(p4sZ1_[1]);

  pTErrs.push_back(pTerrsZ1_[0]);
  pTErrs.push_back(pTerrsZ1_[1]);

    for(unsigned int ifsr1 = 0; ifsr1<p4sZ1ph_.size(); ifsr1++) {
      p4s.push_back(p4sZ1ph_[ifsr1]);
      pTErrs.push_back(pTerrsZ1ph_[ifsr1]);
    }
  return masserror(p4s,pTErrs);

}


vector<TLorentzVector> KinZfitter::GetRefitP4s()
{

  TLorentzVector Z1_1 = p4sZ1REFIT_[0];
  TLorentzVector Z1_2 = p4sZ1REFIT_[1];

  /// fsr photons
  TLorentzVector Z1ph;
  for(unsigned int ifsr1 = 0; ifsr1<p4sZ1phREFIT_.size(); ifsr1++) {
    Z1ph = p4sZ1phREFIT_[ifsr1];
    if(ifsr1==0) Z1_1 = Z1_1 + Z1ph;
    if(ifsr1==1) Z1_2 = Z1_2 + Z1ph;
  }

  vector<TLorentzVector> p4s;
  p4s.push_back(Z1_1); p4s.push_back(Z1_2);

  return p4s;

}

vector<TLorentzVector> KinZfitter::GetP4s()
{

  TLorentzVector Z1_1 = p4sZ1_[0];
  TLorentzVector Z1_2 = p4sZ1_[1];

  // fsr photons
  TLorentzVector Z1ph; 
  for(unsigned int ifsr1 = 0; ifsr1<p4sZ1ph_.size(); ifsr1++) {
    Z1ph = p4sZ1ph_[ifsr1];
    if(ifsr1==0) Z1_1 = Z1_1 + Z1ph;
    if(ifsr1==1) Z1_2 = Z1_2 + Z1ph;
  }

  vector<TLorentzVector> p4s;
  p4s.push_back(Z1_1);
  p4s.push_back(Z1_2);

  return p4s;

}

void KinZfitter::KinRefitZ1()
{

  double l1, l2, lph1, lph2;
  l1 = 1.0; l2 = 1.0; lph1 = 1.0; lph2 = 1.0;

  PerZ1Likelihood(l1, l2, lph1, lph2);
  if(debug_) cout<<"l1 "<<l1<<"; l2 "<<l2<<" lph1 "<<lph1<<" lph2 "<<lph2<<endl;
  SetZ1Result(l1, l2, lph1, lph2);
  if(debug_) cout<<"Z1 refit done"<<endl;

}

int KinZfitter::GetStatus()
{
  return status_;
}

int KinZfitter::GetCovMatStatus()
{
  return covmat_status_;
}

float KinZfitter::GetMinNll()
{
  return minnll_;
}

int KinZfitter::PerZ1Likelihood(double & l1, double & l2, double & lph1, double & lph2)
{

  l1= 1.0; l2 = 1.0;
  lph1 = 1.0; lph2 = 1.0;

  //Declaring start time to help time which part of the refit is the longest
  //clock_t time_start; time_start = static_cast<float>(clock())/CLOCKS_PER_SEC;
  //This code is used to time the refit
  //cout << static_cast<float>(clock())/CLOCKS_PER_SEC - time_start<< endl;


  if(debug_) cout<<"start Z1 refit"<<endl;

  TLorentzVector Z1_1 = p4sZ1_[0];
  TLorentzVector Z1_2 = p4sZ1_[1];

  double RECOpT1 = Z1_1.Pt();
  double RECOpT2 = Z1_2.Pt();

  double pTerrZ1_1 = pTerrsZ1_[0];
  double pTerrZ1_2 = pTerrsZ1_[1];

  if(debug_)cout<<"pT1 "<<RECOpT1<<" pTerrZ1_1 "<<pTerrZ1_1<<endl;
  if(debug_)cout<<"pT2 "<<RECOpT2<<" pTerrZ1_2 "<<pTerrZ1_2<<endl;

  //////////////

  TLorentzVector Z1_ph1, Z1_ph2;
  double pTerrZ1_ph1, pTerrZ1_ph2;
  double RECOpTph1, RECOpTph2;

  TLorentzVector nullFourVector(0, 0, 0, 0);
  Z1_ph1=nullFourVector; Z1_ph2=nullFourVector;
  RECOpTph1 = 0; RECOpTph2 = 0;
  pTerrZ1_ph1 = 0; pTerrZ1_ph2 = 0;


  if(p4sZ1ph_.size()>=1){ // && (idsZ1_[0]==13) ) {
    Z1_ph1 = p4sZ1ph_[0]; pTerrZ1_ph1 = pTerrsZ1ph_[0];
    RECOpTph1 = Z1_ph1.Pt();
    if(debug_) cout<<"put in Z1 fsr photon 1 pT "<<RECOpTph1<<" pT err "<<pTerrZ1_ph1<<endl;
  }
  if(p4sZ1ph_.size()==2){// && (idsZ1_[0]==13) ) {
    if(debug_) cout<<"put in Z1 fsr photon 2"<<endl;
    Z1_ph2 = p4sZ1ph_[1]; pTerrZ1_ph2 = pTerrsZ1ph_[1];
    RECOpTph2 = Z1_ph2.Pt();
  }
  double RECOpTph1min = max(0.0, RECOpTph1-3*pTerrZ1_ph1);
  double RECOpTph2min = max(0.0, RECOpTph2-3*pTerrZ1_ph2);
  double RECOpTph1max = RECOpTph1 < 2 ? RECOpTph1min : RECOpTph1+3*pTerrZ1_ph1;
  double RECOpTph2max = RECOpTph2 < 2 ? RECOpTph2min : RECOpTph2+3*pTerrZ1_ph2;

  RooRealVar* pTph1RECO = new RooRealVar("pTph1RECO", "pTph1RECO", RECOpTph1, 2, 1200);
  RooRealVar* pTph2RECO = new RooRealVar("pTph2RECO", "pTph2RECO", RECOpTph2, 2, 1200);

  RooRealVar* pTph1 = new RooRealVar("pTph1", "pTph1FIT", RECOpTph1, RECOpTph1min, RECOpTph1max);
  RooRealVar* pTph2 = new RooRealVar("pTph2", "pTph2FIT", RECOpTph2, RECOpTph2min, RECOpTph2max);

  RooRealVar* pT1RECO = new RooRealVar("pT1RECO", "pT1RECO", RECOpT1, 5, 1200);
  RooRealVar* pT2RECO = new RooRealVar("pT2RECO", "pT2RECO", RECOpT2, 5, 1200);

  double RECOpT1min = max(5.0, RECOpT1-3*pTerrZ1_1);
  double RECOpT2min = max(5.0, RECOpT2-3*pTerrZ1_2);

  // observables pT1,2,ph1,ph2
  RooRealVar* pT1 = new RooRealVar("pT1", "pT1FIT", RECOpT1, RECOpT1min, RECOpT1+3*pTerrZ1_1 );
  RooRealVar* pT2 = new RooRealVar("pT2", "pT2FIT", RECOpT2, RECOpT2min, RECOpT2+3*pTerrZ1_2 );

  RooRealVar* m1 = new RooRealVar("m1", "m1", Z1_1.M());
  RooRealVar* m2 = new RooRealVar("m2", "m2", Z1_2.M());

  if(debug_) cout<<"m1 "<<m1->getVal()<<" m2 "<<m2->getVal()<<endl;

  double Vtheta1, Vphi1, Vtheta2, Vphi2;
  Vtheta1 = (Z1_1).Theta(); Vtheta2 = (Z1_2).Theta();
  Vphi1 = (Z1_1).Phi(); Vphi2 = (Z1_2).Phi();

  RooRealVar* theta1 = new RooRealVar("theta1", "theta1", Vtheta1);
  RooRealVar* phi1   = new RooRealVar("phi1", "phi1", Vphi1);
  RooRealVar* theta2 = new RooRealVar("theta2", "theta2", Vtheta2);
  RooRealVar* phi2   = new RooRealVar("phi2", "phi2", Vphi2);

  // dot product to calculate (p1+p2+ph1+ph2).M()
  RooFormulaVar E1("E1", "TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1)))+@2*@2)",
       RooArgList(*pT1, *theta1, *m1));
  RooFormulaVar E2("E2", "TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1)))+@2*@2)",
       RooArgList(*pT2, *theta2, *m2));
  if(debug_) cout<<"E1 "<<E1.getVal()<<"; E2 "<<E2.getVal()<<endl;

  /////

  //This code is used to time the refit
  //cout << "After variable def: " << static_cast<float>(clock())/CLOCKS_PER_SEC - time_start<< endl;

  double Vthetaph1, Vphiph1, Vthetaph2, Vphiph2;
  Vthetaph1 = (Z1_ph1).Theta(); Vthetaph2 = (Z1_ph2).Theta();
  Vphiph1 = (Z1_ph1).Phi(); Vphiph2 = (Z1_ph2).Phi();

  RooRealVar* thetaph1 = new RooRealVar("thetaph1", "thetaph1", Vthetaph1);
  RooRealVar* phiph1   = new RooRealVar("phiph1", "phiph1", Vphiph1);
  RooRealVar* thetaph2 = new RooRealVar("thetaph2", "thetaph2", Vthetaph2);
  RooRealVar* phiph2   = new RooRealVar("phiph2", "phi2", Vphiph2);

  RooFormulaVar Eph1("Eph1", "TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1))))",
         RooArgList(*pTph1, *thetaph1));
  RooFormulaVar Eph2("Eph2", "TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1))))",
         RooArgList(*pTph2, *thetaph2));

  //// dot products of 4-vectors

  // 3-vector DOT
  RooFormulaVar* p1v3D2 = new RooFormulaVar("p1v3D2",
              "@0*@1*( ((TMath::Cos(@2))*(TMath::Cos(@3)))/((TMath::Sin(@2))*(TMath::Sin(@3)))+(TMath::Cos(@4-@5)))",
              RooArgList(*pT1, *pT2, *theta1, *theta2, *phi1, *phi2));
  if(debug_) cout<<"p1 DOT p2 is "<<p1v3D2->getVal()<<endl;
  // 4-vector DOT metric 1 -1 -1 -1
  RooFormulaVar p1D2("p1D2", "@0*@1-@2", RooArgList(E1, E2, *p1v3D2));

  //lep DOT fsrPhoton1

  // 3-vector DOT
  RooFormulaVar* p1v3Dph1 = new RooFormulaVar("p1v3Dph1",
                "@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
                RooArgList(*pT1, *pTph1, *theta1, *thetaph1, *phi1, *phiph1));

  // 4-vector DOT metric 1 -1 -1 -1
  RooFormulaVar p1Dph1("p1Dph1", "@0*@1-@2", RooArgList(E1, Eph1, *p1v3Dph1));

  // 3-vector DOT
  RooFormulaVar* p2v3Dph1 = new RooFormulaVar("p2v3Dph1",
                "@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
                RooArgList(*pT2, *pTph1, *theta2, *thetaph1, *phi2, *phiph1));
  // 4-vector DOT metric 1 -1 -1 -1
  RooFormulaVar p2Dph1("p2Dph1", "@0*@1-@2", RooArgList(E2, Eph1, *p2v3Dph1));

  // lep DOT fsrPhoton2

  // 3-vector DOT
  RooFormulaVar* p1v3Dph2 = new RooFormulaVar("p1v3Dph2",
                "@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
                RooArgList(*pT1, *pTph2, *theta1, *thetaph2, *phi1, *phiph2));

  // 4-vector DOT metric 1 -1 -1 -1
  RooFormulaVar p1Dph2("p1Dph2", "@0*@1-@2", RooArgList(E1, Eph2, *p1v3Dph2));

  // 3-vector DOT
  RooFormulaVar* p2v3Dph2 = new RooFormulaVar("p2v3Dph2",
                "@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
                RooArgList(*pT2, *pTph2, *theta2, *thetaph2, *phi2, *phiph2));
  // 4-vector DOT metric 1 -1 -1 -1
  RooFormulaVar p2Dph2("p2Dph2", "@0*@1-@2", RooArgList(E2, Eph2, *p2v3Dph2));

  // fsrPhoton1 DOT fsrPhoton2

  // 3-vector DOT
  RooFormulaVar* ph1v3Dph2 = new RooFormulaVar("ph1v3Dph2",
                 "@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
                 RooArgList(*pTph1, *pTph2, *thetaph1, *thetaph2, *phiph1, *phiph2));
  // 4-vector DOT metric 1 -1 -1 -1
  RooFormulaVar ph1Dph2("ph1Dph2", "@0*@1-@2", RooArgList(Eph1, Eph2, *ph1v3Dph2));

  //This code is used to time the refit
  //cout << "After dot products: " << static_cast<float>(clock())/CLOCKS_PER_SEC - time_start<< endl;


  // mZ1
  RooFormulaVar* mZ1;
  if(p4sZ1ph_.size()==0){
    mZ1 = new RooFormulaVar("mZ1", "TMath::Sqrt(2*@0+@1*@1+@2*@2)", RooArgList(p1D2, *m1, *m2));
  } else if(p4sZ1ph_.size()==1) {
    mZ1 = new RooFormulaVar("mZ1", "TMath::Sqrt(2*@0+2*@1+2*@2+@3*@3+@4*@4)",
          RooArgList(p1D2, p1Dph1, p2Dph1, *m1, *m2));
  } else {
    mZ1 = new RooFormulaVar("mZ1", "TMath::Sqrt(2*@0+2*@1+2*@2+2*@3+2*@4+2*@5+@6*@6+@7*@7)",
          RooArgList(p1D2, p1Dph1, p2Dph1, p1Dph2, p2Dph2, ph1Dph2, *m1, *m2));
  }

  //If mll is outside the bounds of the fit return mll without a fit
  if(mZ1 -> getVal() < 60 || mZ1 -> getVal() > 120){return mZ1 -> getVal();}


  if(debug_) cout<<"mZ1 is "<<mZ1->getVal()<<endl;

  // pTerrs, 1, 2, ph1, ph2
  RooRealVar sigmaZ1_1("sigmaZ1_1", "sigmaZ1_1", pTerrZ1_1);
  RooRealVar sigmaZ1_2("sigmaZ1_2", "sigmaZ1_2", pTerrZ1_2);

  RooRealVar sigmaZ1_ph1("sigmaZ1_ph1", "sigmaZ1_ph1", pTerrZ1_ph1);
  RooRealVar sigmaZ1_ph2("sigmaZ1_ph2", "sigmaZ1_ph2", pTerrZ1_ph2);

  // resolution for decay products
  RooGaussian gauss1("gauss1", "gaussian PDF", *pT1RECO, *pT1, sigmaZ1_1);
  RooGaussian gauss2("gauss2", "gaussian PDF", *pT2RECO, *pT2, sigmaZ1_2);

  RooGaussian gaussph1("gaussph1", "gaussian PDF", *pTph1RECO, *pTph1, sigmaZ1_ph1);
  RooGaussian gaussph2("gaussph2", "gaussian PDF", *pTph2RECO, *pTph2, sigmaZ1_ph2);

  RooProdPdf *model;

  //3 Gaussian Declarations
  RooCBShape* singleCB; RooGaussian* gaussShape1; RooAddPdf* CBplusGauss; RooGaussian* gaussShape2; 
  RooAddPdf* CBplusGaussplusGauss; RooGaussian* gaussShape3; RooAddPdf* CBplusGaussplusGaussplusGauss;

  //Voigtian Declarations
  RooVoigtian* VoiToFit; RooAddPdf* VoiPDF;

  //Variable declarations here to prevent scope issues from if statements below
  RooRealVar meanCB("meanCB", "", meanCB_);
  RooRealVar sigmaCB("sigmaCB", "", sigmaCB_);
  RooRealVar alphaCB("alphaCB", "", alphaCB_);
  RooRealVar nCB("nCB", "", nCB_);
  RooRealVar meanGauss1("meanGauss1", "", meanGauss1_);
  RooRealVar sigmaGauss1("sigmaGauss1", "", sigmaGauss1_);
  RooRealVar f1("f1", "", f1_);
  RooRealVar meanGauss2("meanGauss2", "", meanGauss2_);
  RooRealVar sigmaGauss2("sigmaGauss2", "", sigmaGauss2_);
  RooRealVar f2("f2", "", f2_);
  RooRealVar meanGauss3("meanGauss3", "", meanGauss3_);
  RooRealVar sigmaGauss3("sigmaGauss3", "", sigmaGauss3_);
  RooRealVar f3("f3", "", f3_);
  RooRealVar BWMean("BWMean", "", BWmean_);
  RooRealVar BWGamma("BWGamma", "", BWgamma_);
  RooRealVar meanGauss("meanGauss", "", sigmaValG_);

  //This code is used to time the refit
  //cout << "Before model creation: " << static_cast<float>(clock())/CLOCKS_PER_SEC - time_start<< endl;


  if(threegauss_==true){

  singleCB = new RooCBShape("singleCB", "", *mZ1, meanCB, sigmaCB, alphaCB, nCB);
  gaussShape1 = new RooGaussian("gaussShape1", "", *mZ1, meanGauss1, sigmaGauss1);
  CBplusGauss = new RooAddPdf("CBplusGauss", "", *singleCB, *gaussShape1, f1);
  gaussShape2 = new RooGaussian("gaussShape2", "", *mZ1, meanGauss2, sigmaGauss2);
  CBplusGaussplusGauss = new RooAddPdf("CBplusGaussplusGauss", "", *CBplusGauss, *gaussShape2, f2);
  gaussShape3 = new RooGaussian("gaussShape3", "", *mZ1, meanGauss3, sigmaGauss3);
  CBplusGaussplusGaussplusGauss = new RooAddPdf("CBplusGaussplusGaussplusGauss", "", *CBplusGaussplusGauss, *gaussShape3, f3);

  if(p4sZ1ph_.size()==0){
    model = new RooProdPdf("model", "model", RooArgList(gauss1, gauss2, *CBplusGaussplusGaussplusGauss) );
  } else if(p4sZ1ph_.size()==1) {
    model = new RooProdPdf("model", "model", RooArgList(gauss1, gauss2, gaussph1, *CBplusGaussplusGaussplusGauss) );
  } else {
    model = new RooProdPdf("model", "model", RooArgList(gauss1, gauss2, gaussph1, gaussph2, *CBplusGaussplusGaussplusGauss) );
  }
  VoiToFit = nullptr; VoiPDF = nullptr;


  } else {

  VoiToFit = new RooVoigtian("VoiToFit","Voigtian Fit to m_{ll}", *mZ1, BWMean, BWGamma, meanGauss);
  VoiPDF = new RooAddPdf("VoiShape","",*VoiToFit);
  
  if(p4sZ1ph_.size()==0){
    model = new RooProdPdf("model", "model", RooArgList(gauss1, gauss2, *VoiPDF) );
  } else if(p4sZ1ph_.size()==1) {
    model = new RooProdPdf("model", "model", RooArgList(gauss1, gauss2, gaussph1, *VoiPDF) );
  } else {
    model = new RooProdPdf("model", "model", RooArgList(gauss1, gauss2, gaussph1, gaussph2, *VoiPDF) );
  }


  singleCB = nullptr; gaussShape1 = nullptr; CBplusGauss = nullptr; gaussShape2 = nullptr; 
  CBplusGaussplusGauss = nullptr; gaussShape3 = nullptr; CBplusGaussplusGaussplusGauss = nullptr;

  }


  //This code is used to time the refit
  //cout << "After Model creation: " << static_cast<float>(clock())/CLOCKS_PER_SEC - time_start<< endl;


  // observable set
  RooArgSet *rastmp;
  if(p4sZ1ph_.size()==0){
    rastmp = new RooArgSet(*pT1RECO, *pT2RECO);
  } else if(p4sZ1ph_.size()==1) {
    rastmp = new RooArgSet(*pT1RECO, *pT2RECO, *pTph1RECO);
  } else {
    rastmp = new RooArgSet(*pT1RECO, *pT2RECO, *pTph1RECO, *pTph2RECO);
  }

  RooDataSet* pTs = new RooDataSet("pTs", "pTs", *rastmp);
  pTs->add(*rastmp);

  //RooAbsReal* nll;
  //nll = model->createNLL(*pTs);
  //RooMinuit(*nll).migrad();

  //This code is used to time the refit
  //cout << "Before fit: " << static_cast<float>(clock())/CLOCKS_PER_SEC - time_start<< endl;


  RooFitResult* r = model->fitTo(*pTs, RooFit::Save(), RooFit::PrintLevel(-1));//,RooFit::Timer(true));
  const TMatrixDSym& covMatrix = r->covarianceMatrix();
  status_ = r->status();
  covmat_status_ = r->covQual();
  minnll_ = r -> minNll();

  //This code is used to time the refit
  //cout << "After fit: " << static_cast<float>(clock())/CLOCKS_PER_SEC - time_start<< endl;


  const RooArgList& finalPars = r->floatParsFinal();
  for (int i=0 ; i<finalPars.getSize(); i++) {
//    TString name = dynamic_cast<TString>( dynamic_cast<RooRealVar*>(finalPars.at(i)->GetName()) );
    TString name = static_cast<TString>( finalPars.at(i)->GetName() );

    if(debug_) cout<<"name list of RooRealVar for covariance matrix "<<name<<endl;

  }

  int size = covMatrix.GetNcols();
  //TMatrixDSym covMatrixTest_(size);
  covMatrixZ1_.ResizeTo(size, size);
  covMatrixZ1_ = covMatrix;

  if(debug_) cout<<"save the covariance matrix"<<endl;

  l1 = pT1->getVal()/RECOpT1;
  l2 = pT2->getVal()/RECOpT2;
  double pTerrZ1REFIT1 = pT1->getError();
  double pTerrZ1REFIT2 = pT2->getError();

  pTerrsZ1REFIT_.push_back(pTerrZ1REFIT1);
  pTerrsZ1REFIT_.push_back(pTerrZ1REFIT2);

  double pTerrZ1phREFIT1;
  double pTerrZ1phREFIT2;
  if(p4sZ1ph_.size()>=1) {

    if(debug_) cout<<"set refit result for Z1 fsr photon 1"<<endl;

    lph1 = pTph1->getVal()/RECOpTph1;
    pTerrZ1phREFIT1 = pTph1->getError();
    if(debug_) cout<<"scale "<<lph1<<" pterr "<<pTerrZ1phREFIT1<<endl;

    pTerrsZ1phREFIT_.push_back(pTerrZ1phREFIT1);

  }
  if(p4sZ1ph_.size()==2){

    lph2 = pTph2->getVal()/RECOpTph2;
    pTerrZ1phREFIT2 = pTph2->getError();
    pTerrsZ1phREFIT_.push_back(pTerrZ1phREFIT2);

  }

  //This code is used to time the refit
  //cout << "Before delete statements: " << static_cast<float>(clock())/CLOCKS_PER_SEC - time_start<< endl;


 // delete pT1; delete pT2;delete pTph1; delete pTph2;
 //delete pT1RECO; delete pT2RECO; delete pTph1RECO; delete pTph2RECO;  delete m1; delete m2; delete theta1; 
 // delete p1v3D2; delete p1v3Dph1; delete p2v3Dph2; delete p2v3Dph2; delete ph1v3Dph2; delete mZ1;
  //delete model; delete rastmp; delete pTs; delete r;
   
  //delete nll;
  delete pTs;
  delete rastmp;
  delete r;
  delete singleCB; delete gaussShape1; delete CBplusGauss; delete gaussShape2; delete CBplusGaussplusGauss; delete gaussShape3; delete CBplusGaussplusGaussplusGauss;
  delete VoiToFit; delete VoiPDF;
  delete model;
  delete mZ1;
  delete ph1v3Dph2; delete p2v3Dph2; delete p1v3Dph2; delete p2v3Dph1; delete p1v3Dph1; delete p1v3D2;
  delete thetaph1; delete phiph1; delete thetaph2; delete phiph2; delete pTph1; delete pTph2;

  delete theta1; delete phi1; delete theta2; delete phi2;
  delete m1; delete m2; delete pTph1RECO; delete pTph2RECO;
  delete pT1; delete pT2;  delete pT1RECO; delete pT2RECO;

  //This code is used to time the refit
  //cout << "After delete statements: " << static_cast<float>(clock())/CLOCKS_PER_SEC - time_start<< endl;


  if(debug_) cout<<"end Z1 refit"<<endl;

  return 0;

}

#endif
