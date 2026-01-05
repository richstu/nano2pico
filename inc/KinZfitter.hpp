/*************************************************************************
*  Authors:   Tongguang CHeng(IHEP, Beijing) Hualin Mei(UF)
*************************************************************************/
#ifndef KinZfitter_h
#define KinZfitter_h

// C++ includes
#include <iostream>
#include <complex>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>

// ROOT includes
#include "TString.h"
#include "TLorentzVector.h"

// ROOFIT
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooCBShape.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooVoigtian.h"
#include "RooMsgService.h"


// fit result covariance matrix
#include <TMatrixDSym.h>

// CMSSW related pT error calculator
//#include "HelperFunction/interface/HelperFunction.h"
//#include "DataFormats/Candidate/interface/Candidate.h"


using namespace std;

class KinZfitter {
 public:

  //Constructors
  KinZfitter();
  KinZfitter(TString pdf_filename);

  /// Kinematic fit of lepton momenta
  void Setup(std::map<unsigned int, TLorentzVector> selectedLeptons, std::map<unsigned int, TLorentzVector> selectedFsrPhotons, std::map<unsigned int, double> errorLeptons);

  ///
  void KinRefitZ1();

  int  PerZ1Likelihood(double & l1, double & l2, double & lph1, double & lph2);
  void SetZ1Result(double l1, double l2, double lph1, double lph2);
  double lep_pterr(TLorentzVector lepton, double leptonError);
  double pterr(TLorentzVector fsrPhoton);
  double masserror(std::vector<TLorentzVector> p4s, std::vector<double> pTErrs);
  double masserrorFullCov(std::vector<TLorentzVector> p4s, TMatrixDSym covMatrix);

  // result wrappers
  double GetRefitM4l();
  double GetM4l();
  double GetRefitMZ1();

  double GetMZ1Err();
  double GetRefitM4lErr();
  double GetM4lErr();
  double GetRefitM4lErrFullCov();

  float GetMinNll();
  int   GetStatus();
  int   GetCovMatStatus();

  std::vector<TLorentzVector> GetRefitP4s();
  std::vector<TLorentzVector> GetP4s();

  ////////////////

  // Need to check in deep whether this function is doable
  // RooFormulaVar p1DOTp2(RooRealVar pT1, RooRealVar theta1, RooRealVar phi1, RooRealVar m1, TString index1, RooRealVar pT2, RooRealVar theta2, RooRealVar phi2, RooRealVar m2, TString index2);

 private:

  /// True mZ/mZ1 shape
  TString PDFName_;

  /// debug flag
  bool debug_;

  /// whether use correction for pT error
  bool isCorrPTerr_;
  /// whether use data or mc correction
  bool isData_;

  void initZs(std::map<unsigned int, TLorentzVector> selectedLeptons, std::map<unsigned int, TLorentzVector> selectedFsrPhoton, std::map<unsigned int, double> errorLeptons);

  /// lepton ids that fsr photon associated to
  std::vector<int> idsFsrZ1_;
  /// (Four) TLorentzVectors that form the Higgs Candidate
  std::vector<TLorentzVector> p4sZ1_, p4sZ1ph_;
  std::vector<TLorentzVector> p4sZ1REFIT_, p4sZ1phREFIT_;

  /// pTerr vector
  std::vector<double> pTerrsZ1_, pTerrsZ1ph_;
  std::vector<double> pTerrsZ1REFIT_, pTerrsZ1phREFIT_;

  // covariance matrix
  // what directly coming from Refit
  TMatrixDSym covMatrixZ1_;

  //status' returned for fits
  int status_;
  int covmat_status_;
  double minnll_;

  // refit energy scale with respect to reco pT
  double lZ1_l1_, lZ1_l2_;
  double lZ1_ph1_, lZ1_ph2_;

  // True mZ1 shape parameters
  //double sgVal_, aVal_, nVal_, fVal_, meanVal_, sigmaVal_, f1Val_;
  //double sgValCB_, aValCB_, nValCB_, meanValG_, sigmaValG_, f1Val_, f2Val_, meanBW_, gammaBW_;
  double sigmaValG_, BWmean_, BWgamma_;

  bool threegauss_;
  //double aVal_, meanVal_, nVal_, sgVal_, bwsigVal_;
  double meanCB_, sigmaCB_, alphaCB_, nCB_, meanGauss1_, sigmaGauss1_, f1_, meanGauss2_, sigmaGauss2_, f2_, meanGauss3_, sigmaGauss3_, f3_;

};

#endif
