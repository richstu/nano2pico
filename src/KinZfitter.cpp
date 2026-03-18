/*************************************************************************
 *  Authors:   Tongguang Cheng
 *  Adapted for H-->Zgamma by: Prasanna Siddireddy and Anders Barzdukas 
 *  Contact at abarzdukas@ucsb.edu
 *************************************************************************/
#ifndef KinZfitter_cpp
#define KinZfitter_cpp

/// KinFitter header
#include "KinZfitter.hpp"
#include "RooHelpers.h"
#include <chrono>

void KinZfitter::set_consts(double ptl1, double phil1, double etal1, double sigmal1, double ml1,
                            double ptl2, double phil2, double etal2, double sigmal2, double ml2,
                            unsigned int nfsrph,
                            double ptg3, double phig3, double etag3, double sigmag3,
                            double ptg4, double phig4, double etag4, double sigmag4){
  pTl1_      = ptl1;
  pTl2_      = ptl2;
  phil1_     = phil1;
  phil2_     = phil2;
  etal1_     = etal1;
  etal2_     = etal2;
  if(lepid_ == 11){//Left in place in case we want to investigate lepton split refit in future.
    sigmal1_   = sigmal1;
    sigmal2_   = sigmal2;
  }else {
    sigmal1_   = sigmal1;
    sigmal2_   = sigmal2;
  }
  ml1_       = ml1;
  ml2_       = ml2;
  if(nfsrph>0){
    pTg3_    = ptg3;
    phig3_   = phig3;
    etag3_   = etag3;
    sigmag3_ = sigmag3;
  }
  if(nfsrph>1){
    pTg4_    = ptg4;
    phig4_   = phig4;
    etag4_   = etag4;
    sigmag4_ = sigmag4;
  }
}

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


void KinZfitter::setEs(double pT1, double pT2, unsigned int nfsrph, double pT3, double pT4){
  En1_ = sqrt(pow(pT1,2)*(1 + pow(sinh(etal1_),2)) + pow(ml1_,2));
  En2_ = sqrt(pow(pT2,2)*(1 + pow(sinh(etal2_),2)) + pow(ml2_,2));
  En3_ = 0;
  En4_ = 0;
  if(nfsrph>0) En3_ = sqrt(pow(pT3,2)*(1 + pow(sinh(etag3_),2)));
  if(nfsrph>1) En4_ = sqrt(pow(pT4,2)*(1 + pow(sinh(etag4_),2)));
}

void KinZfitter::setmZ(double pT1, double pT2, unsigned int nfsrph, double pT3, double pT4){
  mll_ = 0;
  if(nfsrph == 0){
    mll_ = sqrt(pow((En1_ + En2_),2) - pow((pT1*cos(phil1_) + pT2*cos(phil2_)),2) 
          - pow((pT1*sin(phil1_) + pT2*sin(phil2_)),2) - pow((pT1*sinh(etal1_) + pT2*sinh(etal2_)),2));
  } else if (nfsrph == 1){
    mll_ = sqrt(pow((En1_ + En2_ + En3_),2) - pow((pT1*cos(phil1_) + pT2*cos(phil2_) + pT3*cos(phig3_)),2)
                                            - pow((pT1*sin(phil1_) + pT2*sin(phil2_) + pT3*sin(phig3_)),2)
                                            - pow((pT1*sinh(etal1_) + pT2*sinh(etal2_) + pT3*sinh(etag3_)),2));
  } else if (nfsrph ==2){
    mll_ = sqrt(pow((En1_ + En2_ + En3_ + En4_),2) - pow((pT1*cos(phil1_) + pT2*cos(phil2_) + pT3*cos(phig3_) + pT4*cos(phig4_)),2)
                                                   - pow((pT1*sin(phil1_) + pT2*sin(phil2_) + pT3*sin(phig3_) + pT4*sin(phig4_)),2)
                                                   - pow((pT1*sinh(etal1_) + pT2*sinh(etal2_) + pT3*sinh(etag3_) + pT4*sinh(etag4_)),2));
  }
}

double KinZfitter::gaussian(double x, double mu, double sigma){
  return exp(-0.5*pow(x - mu,2)/pow(sigma,2))/(sigma*sqrt(2*PI));
}

void KinZfitter::evaluateShape(double pT1, double pT2, unsigned int nfsrph, double pT3, double pT4){
  if(nfsrph == 0){
    setEs(pT1, pT2, 0);
    setmZ(pT1, pT2, 0);
  } else if (nfsrph == 1){
    setEs(pT1, pT2, 1, pT3);
    setmZ(pT1, pT2, 1, pT3);
  } else if (nfsrph == 2){
    setEs(pT1, pT2, 2, pT3, pT4);
    setmZ(pT1, pT2, 2, pT3, pT4);
  }
  double mll = mll_;
  double gausst1, gausst2, gausst3;
  gausst1 = gaussian(mll, meanGauss1_, sigmaGauss1_);
  gausst2 = gaussian(mll, meanGauss2_, sigmaGauss2_);
  gausst3 = gaussian(mll, meanGauss3_, sigmaGauss3_);
  double CBA = pow(nCB_/alphaCB_,nCB_)*exp(-pow(alphaCB_,2)/2);
  double CBB = nCB_/alphaCB_ - alphaCB_; 
  double CBC = (nCB_/alphaCB_)*(1/(nCB_-1))*exp(-pow(alphaCB_,2)/2);
  double CBD = sqrt(PI/2)*(1+erf(alphaCB_/sqrt(2)));

  double CBt = 0;
  if(mll - meanCB_ > -1*alphaCB_ * sigmaCB_){
    CBt = 1*exp(-0.5*pow(mll - meanCB_,2)/pow(sigmaCB_,2));
  } else if(mll - meanCB_ <= -1*alphaCB_ * sigmaCB_){
    CBt = 1*CBA*pow((CBB - (mll-meanCB_)/sigmaCB_),-1*nCB_);
  }
  CBt = CBt/(sigmaCB_*(CBC+CBD));
  shapeEval_ = (((f1_*CBt + (1-f1_)*gausst1)*f2_ + (1-f2_)*gausst2)*f3_ + (1-f3_)*gausst3);
}

double KinZfitter::NLL_0(const double *pTs){
  double pTl1r = pTs[0];
  double pTl2r = pTs[1];
  double gauss1, gauss2, full, NLL;
  setEs(pTl1r, pTl2r, 0);
  setmZ(pTl1r, pTl2r, 0);

  gauss1 = gaussian(pTl1r, pTl1_, sigmal1_);
  gauss2 = gaussian(pTl2r, pTl2_, sigmal2_);
  
  evaluateShape(pTl1r, pTl2r, 0); 

  full = gauss1*gauss2*shapeEval_;//Normalized
  NLL = -log(full);
  return NLL;
}

double KinZfitter::NLL_1(const double *pTs){
  double pTl1r = pTs[0];
  double pTl2r = pTs[1];
  double pTg3r = pTs[2];

  double gauss1, gauss2, gauss3, full, NLL;
  setEs(pTl1r, pTl2r, 1, pTg3r);
  setmZ(pTl1r, pTl2r, 1, pTg3r);

  gauss1 = gaussian(pTl1r, pTl1_, sigmal1_);
  gauss2 = gaussian(pTl2r, pTl2_, sigmal2_);
  gauss3 = gaussian(pTg3r, pTg3_, sigmag3_);

  evaluateShape(pTl1r, pTl2r, 1, pTg3r);

  full = gauss1*gauss2*gauss3*shapeEval_;//Normalized
  NLL = -log(full);
  return NLL;
}

double KinZfitter::NLL_2(const double *pTs){
  double pTl1r = pTs[0];
  double pTl2r = pTs[1];
  double pTg3r = pTs[2];
  double pTg4r = pTs[3];

  double gauss1, gauss2, gauss3, gauss4, full, NLL;
  setEs(pTl1r, pTl2r, 2, pTg3r, pTg4r);
  setmZ(pTl1r, pTl2r, 2, pTg3r, pTg4r);

  gauss1 = gaussian(pTl1r, pTl1_, sigmal1_);
  gauss2 = gaussian(pTl2r, pTl2_, sigmal2_);
  gauss3 = gaussian(pTg3r, pTg3_, sigmag3_);
  gauss4 = gaussian(pTg4r, pTg4_, sigmag4_);

  evaluateShape(pTl1r, pTl2r, 2, pTg3r, pTg4r);
  
  full = gauss1*gauss2*gauss3*gauss4*shapeEval_;//Normalized
  NLL = -log(full);
  return NLL;
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


void KinZfitter::Setup(std::map<unsigned int, TLorentzVector> selectedLeptons, std::map<unsigned int, TLorentzVector> selectedFsrPhotons, std::map<unsigned int, double> errorLeptons, int lepid) {

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
  lepid_ = lepid;
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
  RooHelpers::LocalChangeMsgLevel changeMsgLvl(RooFit::ERROR);
  //Temporarily suppressing warnings to test n2p overall
  //KinZFitter needs some reworking to avoid sigma going negative
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
  double RECOpT1min = max(5.0, RECOpT1-3*pTerrZ1_1);
  double RECOpT2min = max(5.0, RECOpT2-3*pTerrZ1_2);
  double RECOpT1max = RECOpT1+3*pTerrZ1_1;//why??
  double RECOpT2max = RECOpT2+3*pTerrZ1_2;

  //auto startTime = std::chrono::steady_clock::now();
  //Set global fit constants
  if(p4sZ1ph_.size()==0){
    set_consts(Z1_1.Pt(), Z1_1.Phi(), Z1_1.Eta(), pTerrZ1_1, Z1_1.M(),
               Z1_2.Pt(), Z1_2.Phi(), Z1_2.Eta(), pTerrZ1_2, Z1_2.M(),
               p4sZ1ph_.size());
    setEs(Z1_1.Pt(), Z1_2.Pt(), 0);
    setmZ(Z1_1.Pt(), Z1_2.Pt(), 0);
  } else if(p4sZ1ph_.size()==1){
    set_consts(Z1_1.Pt(), Z1_1.Phi(), Z1_1.Eta(), pTerrZ1_1, Z1_1.M(),
               Z1_2.Pt(), Z1_2.Phi(), Z1_2.Eta(), pTerrZ1_2, Z1_2.M(),
               p4sZ1ph_.size(),
               Z1_ph1.Pt(), Z1_ph1.Phi(), Z1_ph1.Eta(), pTerrZ1_ph1);
    setEs(Z1_1.Pt(), Z1_2.Pt(), 1, Z1_ph1.Pt());
    setmZ(Z1_1.Pt(), Z1_2.Pt(), 1, Z1_ph1.Pt());
  } else if(p4sZ1ph_.size()==2){
    set_consts(Z1_1.Pt(), Z1_1.Phi(), Z1_1.Eta(), pTerrZ1_1, Z1_1.M(),
               Z1_2.Pt(), Z1_2.Phi(), Z1_2.Eta(), pTerrZ1_2, Z1_2.M(),
               p4sZ1ph_.size(),
               Z1_ph1.Pt(), Z1_ph1.Phi(), Z1_ph1.Eta(), pTerrZ1_ph1,
               Z1_ph2.Pt(), Z1_ph2.Phi(), Z1_ph2.Eta(), pTerrZ1_ph2);
    setEs(Z1_1.Pt(), Z1_2.Pt(), 1, Z1_ph1.Pt(), Z1_ph2.Pt());
    setmZ(Z1_1.Pt(), Z1_2.Pt(), 1, Z1_ph1.Pt(), Z1_ph2.Pt());
  }

  //If mll is outside the bounds of the fit return mll without a fit
  if(mll_ < 60 || mll_ > 120) return mll_;

  int status_my = -1;
  int covstatus_my = -1;
  double minnll_my = 1;
  const double *xs;
  const double *errs;
  xs = 0; //dodging -Werror=maybe-uninitialized, these will not be used unless redefined later
  errs = 0;
  ROOT::Math::Minimizer* minimum =
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

  //Do the fits
  if(p4sZ1ph_.size()==0){
    minimum->SetMaxFunctionCalls(1000); // for Minuit/Minuit2
    minimum->SetTolerance(1);
    minimum->SetPrintLevel(-1);
 
    // create gradfunctor
    //ROOT::Math::GradFunctor f(this,&KinZfitter::NLL_0,&KinZfitter::gradNLL_0,2);
    ROOT::Math::Functor f(this, &KinZfitter::NLL_0,2);
    double step[2] = {0.5,0.5};
    double variable[2] = {Z1_1.Pt(),Z1_2.Pt()};
    minimum->SetFunction(f);
 
    minimum->SetVariable(0,"pt1",variable[0], step[0]);
    minimum->SetVariable(1,"pt2",variable[1], step[1]);
    minimum->SetVariableLimits(0,RECOpT1min,RECOpT1max);
    minimum->SetVariableLimits(1,RECOpT2min,RECOpT2max);
    minimum->Minimize();
    minimum->Hesse();
    xs = minimum->X();
    errs = minimum->Errors();
    status_my = minimum->Status();
    covstatus_my = minimum->CovMatrixStatus();
    minnll_my = minimum->MinValue();
  } else if(p4sZ1ph_.size()==1){
    minimum->SetMaxFunctionCalls(1000); // for Minuit/Minuit2
    minimum->SetTolerance(1);
    minimum->SetPrintLevel(-1);
    // create gradfunctor
    //ROOT::Math::GradFunctor f(this,&KinZfitter::NLL_1,&KinZfitter::gradNLL_1,3);
    ROOT::Math::Functor f(this, &KinZfitter::NLL_1,3);
    double step[3] = {0.5,0.5,0.5};
    double variable[3] = {Z1_1.Pt(),Z1_2.Pt(),Z1_ph1.Pt()};
    minimum->SetFunction(f);

    minimum->SetVariable(0,"pt1",variable[0], step[0]);
    minimum->SetVariable(1,"pt2",variable[1], step[1]);
    minimum->SetVariable(2,"pt3",variable[2], step[2]);
    minimum->SetVariableLimits(0,RECOpT1min,RECOpT1max);
    minimum->SetVariableLimits(1,RECOpT2min,RECOpT2max);
    minimum->SetVariableLimits(2,RECOpTph1min,RECOpTph1max);
    
    minimum->Minimize();
    minimum->Hesse();
    xs = minimum->X();
    errs = minimum->Errors();
    status_my = minimum->Status();
    covstatus_my = minimum->CovMatrixStatus();
    minnll_my = minimum->MinValue();
  } else if(p4sZ1ph_.size()==2){
    minimum->SetMaxFunctionCalls(1000); // for Minuit/Minuit2
    minimum->SetTolerance(1);
    minimum->SetPrintLevel(-1);
  
    // create gradfunctor
    //ROOT::Math::GradFunctor f(this,&KinZfitter::NLL_2,&KinZfitter::gradNLL_2,4);
    ROOT::Math::Functor f(this, &KinZfitter::NLL_2,4);
    double step[4] = {0.5,0.5,0.5,0.5};
    double variable[4] = {Z1_1.Pt(),Z1_2.Pt(),Z1_ph1.Pt(),Z1_ph2.Pt()};
    minimum->SetFunction(f);
  
    minimum->SetVariable(0,"pt1",variable[0], step[0]);
    minimum->SetVariable(1,"pt2",variable[1], step[1]);
    minimum->SetVariable(2,"pt3",variable[2], step[2]);
    minimum->SetVariable(3,"pt4",variable[3], step[3]);
    minimum->SetVariableLimits(0,RECOpT1min,RECOpT1max);
    minimum->SetVariableLimits(1,RECOpT2min,RECOpT2max);
    minimum->SetVariableLimits(2,RECOpTph1min,RECOpTph1max);
    minimum->SetVariableLimits(3,RECOpTph2min,RECOpTph2max);

    minimum->Minimize();
    minimum->Hesse();
    xs = minimum->X();
    errs = minimum->Errors();
    status_my = minimum->Status();
    covstatus_my = minimum->CovMatrixStatus();
    minnll_my = minimum->MinValue();
  }

  status_ = status_my;
  covmat_status_ = covstatus_my;
  minnll_ = minnll_my;
  if(debug_) cout<<"save the covariance matrix"<<endl;
  double pTerrZ1REFIT1;
  double pTerrZ1REFIT2;
  l1 = xs[0]/RECOpT1;
  l2 = xs[1]/RECOpT2;
  pTerrZ1REFIT1 = errs[0];
  pTerrZ1REFIT2 = errs[1];

  //cout<<"l1,l2: "<<xs[0]<<", "<<xs[1]<<endl;
  pTerrsZ1REFIT_.push_back(pTerrZ1REFIT1);
  pTerrsZ1REFIT_.push_back(pTerrZ1REFIT2);

  double pTerrZ1phREFIT1;
  double pTerrZ1phREFIT2;
  if(p4sZ1ph_.size()>=1) {

    if(debug_) cout<<"set refit result for Z1 fsr photon 1"<<endl;
    lph1 = xs[2]/RECOpTph1;
    pTerrZ1phREFIT1 = errs[2];
    if(debug_) cout<<"scale "<<lph1<<" pterr "<<pTerrZ1phREFIT1<<endl;

    pTerrsZ1phREFIT_.push_back(pTerrZ1phREFIT1);

  }
  if(p4sZ1ph_.size()==2){
    lph2 = xs[3]/RECOpTph2;
    pTerrZ1phREFIT2 = errs[3];
    pTerrsZ1phREFIT_.push_back(pTerrZ1phREFIT2);

  }

  //This code is used to time the refit
  //cout << "After delete statements: " << static_cast<float>(clock())/CLOCKS_PER_SEC - time_start<< endl;


  if(debug_) cout<<"end Z1 refit"<<endl;

  return 0;

}

#endif
