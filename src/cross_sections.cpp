// Lists the MC cross sections

#include <iostream>
#include "cross_sections.hpp"

using namespace std;

namespace xsec{

  float crossSection(const TString &file, int year){
    float xsec(-999999.), Htobb(0.5824), HToZG(0.001533), ZToLL(0.100974);
    float HToGG(0.00227), HToMM(0.000218), HToTT(0.06256);
    float HToZZ(0.02641), ZToQQ(0.69911), ZToNuNu(0.2);
    float WToLNu(0.1086), HToWW(0.2152);

    if (year == 2016) {
        // Cross-section taken from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
        // Alternative option: https://twiki.cern.ch/twiki/bin/view/Sandbox/FullNNLOcrossSections#Top_cross_section_for_13_TeV
        if(file.Contains("TT_TuneCUETP8M2T4")  xsec = 815.96;
        if(file.Contains("TTJets_HT")){//LO cross sections with k-factor of 1.625 already applied
          if(file.Contains("2500toInf")) xsec = 0.0023234211;
          if(file.Contains("1200to2500")) xsec = 0.194972521;
          if(file.Contains("800to1200")) xsec = 1.07722318;
          if(file.Contains("600to800")) xsec = 2.61537118;
        }
	if(file.Contains("TTJets_DiLept")) xsec = 831.8*0.105; //source: PDG XS*BR http://pdg.lbl.gov/2019/reviews/rpp2018-rev-top-quark.pdf ; 85.66 in Humboldtv3+earlier
        if(file.Contains("TTJets_SingleLept")) xsec = 831.8*0.219; //source: PDG XS*BR http://pdg.lbl.gov/2019/reviews/rpp2018-rev-top-quark.pdf ; 178.7 in Humboldtv3+earlier
        
        if(file.Contains("TTJets_DiLept_genMET-150")) xsec = 3.638*1.627; //XSDB*XSDB/GenXSecAnalyzer K-Factor; 0.0676543*85.66 in Humboldtv3+earlier
        if(file.Contains("TTJets_SingleLept") && file.Contains("genMET-150")) xsec = 5.952*1.627; //GenXSecAnalyzer*XSB K-Factor; 0.0568246*178.7 in Humboldtv3+earlier
        

        // cross sections from mcm
        if(file.Contains("TTG")) xsec = 3.697;                
        if(file.Contains("TTTT_Tune")) xsec = 0.009103; //twiki table
        
        // From https://cms-pdmv.cern.ch/mcm
        // k-factors from https://mangano.web.cern.ch/mangano/public/MECCA/samples_50ns_miniaod.txt
        // k-factors are ratio of https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
        // NLO/NNLO cross-sections to that of an inclusive sample in mcm at lower order (LO/NLO)

        if(file.Contains("WJetsToLNu_TuneCUETP8M1")) xsec=61526.7; //NNLO from Lesya's summary table 

        //cross-section per slice changed due to change in genHT definition
        if(file.Contains("WJetsToLNu_HT-70To100_TuneCUETP8M1"))    xsec = 1353.*1.224; //XSDB*XSDB K-factor; 1372.*1.21 in Humboldtv3+earlier
        if(file.Contains("WJetsToLNu_HT-100To200_TuneCUETP8M1"))   xsec = 1346.*1.224; //XSDB*XSDB K-factor; 1347.*1.21 in Humboldtv3+earlier
        if(file.Contains("WJetsToLNu_HT-200To400_TuneCUETP8M1"))   xsec = 360.1*1.224; //XSDB*XSDB K-factor; 360.*1.21 in Humboldtv3+earlier
        if(file.Contains("WJetsToLNu_HT-400To600_TuneCUETP8M1"))   xsec = 48.8*1.224; //XSDB*XSDB K-factor; 48.98*1.21 in Humboldtv3+earlier
        if(file.Contains("WJetsToLNu_HT-600ToInf_TuneCUETP8M1"))   xsec = 18.77*1.224; //XSDB*XSDB K-factor; 18.77*1.2 in Humboldtv3+earlier
        if(file.Contains("WJetsToLNu_HT-600To800_TuneCUETP8M1"))   xsec = 12.07*1.224; //XSDB*XSDB K-factor; 12.05*1.2 in Humboldtv3+earlier
        if(file.Contains("WJetsToLNu_HT-800To1200_TuneCUETP8M1"))  xsec = 5.497*1.224; //XSDB*XSDB K-factor; 5.501*1.2 in Humboldtv3+earlier
        if(file.Contains("WJetsToLNu_HT-1200To2500_TuneCUETP8M1")) xsec = 1.329*1.224; //XSDB*XSDB K-factor; 1.329*1.2 in Humboldtv3+earlier
        if(file.Contains("WJetsToLNu_HT-2500ToInf_TuneCUETP8M1"))  xsec = 0.03209*1.224; //XSDB*XSDB K-factor; 0.03216*1.2 in Humboldtv3+earlier

        //updated 02-2019 with XSDB
        if(file.Contains("QCD_HT100to200_TuneCUETP8M1"))   xsec = 28060000;
        if(file.Contains("QCD_HT200to300_TuneCUETP8M1"))   xsec = 1710000;
        if(file.Contains("QCD_HT300to500_TuneCUETP8M1"))   xsec = 347500;
        if(file.Contains("QCD_HT500to700_TuneCUETP8M1"))   xsec = 32060;
        if(file.Contains("QCD_HT700to1000_TuneCUETP8M1"))  xsec = 6829;
        if(file.Contains("QCD_HT1000to1500_TuneCUETP8M1")) xsec = 1207;
        if(file.Contains("QCD_HT1500to2000_TuneCUETP8M1")) xsec = 120.0;
        if(file.Contains("QCD_HT2000toInf_TuneCUETP8M1"))  xsec = 25.25;

        // Cross sections from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec
        // multiplied by BF(W->mu,e,tau) = 0.324
        if (file.Contains("ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1"))     xsec = 3.365; //from XSDB, 3.4 in Humboldtv3+earlier
        if (file.Contains("ST_t-channel_antitop_4f_InclusiveDecays")) xsec = 67.91;
        if (file.Contains("ST_t-channel_top_4f_InclusiveDecays")) xsec = 113.3;
        if (file.Contains("ST_tW_antitop_5f_NoFullyHadronicDecays"))     xsec = 38.06*0.530775; //from XSDB*PDG BR; 35.85*0.543 in Humboldtv3+earlier
        if (file.Contains("ST_tW_top_5f_NoFullyHadronicDecays"))     xsec = 38.09*0.530775; //from XSB*PDG BR; 35.85*0.543 in Humboldtv3+earlier

        if(file.Contains("DYJetsToLL_M-10to50_TuneCUETP8M1")) xsec = 18610*1.23;
	if(file.Contains("DYJetsToLL_M-50_TuneCUETP8M1"))     xsec = 4895*1.23;
	if(file.Contains("DYJetsToLL_M-50_TuneCP5"))          xsec = 6077.22;

        if(file.Contains("DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1"))    xsec = 175.3*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1"))   xsec = 139.4*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1"))   xsec = 42.75*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1"))   xsec = 5.497*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1"))   xsec = 1.363*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1"))  xsec = 0.6759*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1")) xsec = 0.116*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1"))  xsec = 0.002592*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1"))   xsec = 2.21*1.23;

        if(file.Contains("ZJetsToNuNu_HT-100To200"))       xsec = 280.35*1.2245; //twiki table * XSDB K-factor; 280.35*1.27 in Humboldtv3+earlier
        if(file.Contains("ZJetsToNuNu_HT-200To400"))       xsec = 77.67*1.2245; // twiki table * XSDB K-factor; 77.67*1.27 in Humboldtv3+earlier
        if(file.Contains("ZJetsToNuNu_HT-400To600"))       xsec = 10.73*1.2245; // twiki table * XSDB K-factor; 10.73*1.27 in Humboldtv3+earlier
        if(file.Contains("ZJetsToNuNu_HT-600To800"))       xsec = 2.559*1.2245; // twiki table * XSDB K-factor; 2.536*1.27 in Humboldtv3+earlier
        if(file.Contains("ZJetsToNuNu_HT-800To1200"))      xsec = 1.796*1.2245; // twiki table * XSDB K-factor; 1.161*1.27 in Humboldtv3+earlier
        if(file.Contains("ZJetsToNuNu_HT-1200To2500"))     xsec = 0.28833*1.2245; // twiki table * XSDB K-factor; 0.2824*1.27 in Humboldtv3+earlier
        if(file.Contains("ZJetsToNuNu_HT-2500ToInf"))      xsec = 0.006945*1.2245; // twiki table * XSDB K-factor; 0.006459*1.27 in Humboldtv3+earlier
        if(file.Contains("ZJetsToNuNu_HT-600ToInf"))       xsec = 3.986*1.2245; // twiki table * XSDB K-factor; 3.986*1.27 in Humboldtv3+earlier

        if(file.Contains("TTZToQQ_TuneCUETP8M1"))          xsec = 0.5297; //XSDB
        if(file.Contains("TTZToLLNuNu_M-10_TuneCUETP8M1")) xsec = 0.2529; //XSDB
        if(file.Contains("TTWJetsToQQ_TuneCUETP8M1"))      xsec = 0.4062; //TWiki table (NLO)
        if(file.Contains("TTWJetsToLNu_TuneCUETP8M1"))     xsec = 0.2043; //TWiki table (NLO)
       
        //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
        if(file.Contains("WWTo2L2Nu_13TeV-powheg"))        xsec = 12.178; //NNLO
        if(file.Contains("WWToLNuQQ_13TeV-powheg"))        xsec = 49.997; //NNLO; RA2b tuples use WWTo1L1Nu2Q

        if(file.Contains("WZTo1L3Nu"))             xsec = 3.054; // from XSDB(LO); 3.05 in Humboldtv3+earlier
        if(file.Contains("WZTo1L1Nu2Q"))           xsec = 10.73; // from XSDB(LO); 10.96 in Humbdoltv3+earlier
        if(file.Contains("WZTo2L2Q"))              xsec = 5.595;
        if(file.Contains("WZTo3LNu_TuneCUETP8M1")) xsec = 4.42965;
        if(file.Contains("ZZ_TuneCUETP8M1"))       xsec = 16.523;
	if(file.Contains("ZZ_TuneCP5"))            xsec = 16.91;

        // Calculated at 13 TeV in
        // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV
        // Higgs branching ratios from
        // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR
        if(file.Contains("ZH_HToBB_ZToNuNu_M125"))  xsec = 0.883*Htobb*0.2;
        if(file.Contains("WH_HToBB_WToLNu_M125"))   xsec = 1.373*Htobb*(0.1071+0.1063+0.1138);

    } else {

        //  Cross-section taken from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
        // Alternative option: https://twiki.cern.ch/twiki/bin/view/Sandbox/FullNNLOcrossSections#Top_cross_section_for_13_TeV
        if(file.Contains("TTJets_DiLept")) xsec = 831.8*0.105; //source: PDG XS*BR http://pdg.lbl.gov/2019/reviews/rpp2018-rev-top-quark.pdf ; 85.66 in Humboldtv3+earlier
        if(file.Contains("TTJets_SingleLept")) xsec = 831.8*0.219; //source: PDG XS*BR http://pdg.lbl.gov/2019/reviews/rpp2018-rev-top-quark.pdf ; 178.7 in Humboldtv3+earlier
        
        if(file.Contains("TTJets_DiLept_genMET-150")) xsec = 3.655*1.679; //XSDB*XSDB/GenXSecAnalyzer K-Factor; 0.0676543*85.66 in Humboldtv3+earlier
        if(file.Contains("TTJets_SingleLept") && file.Contains("genMET-150")) xsec = 6.179*1.679; //GenXSecAnalyzer*XSB K-Factor; 0.0568246*178.7 in Humboldtv3+earlier

        if(file.Contains("TTJets_DiLept_genMET-80")) xsec = 22.46*1.677; // GenXSecAnalyzer*XSDB K-Factor; 0.412882*85.66 in Humboldtv3+earlier
        if(file.Contains("TTJets_SingleLept") && file.Contains("genMET-80")) xsec = 32.23*1.677; // GenXSecAnalyze*XSDB K-factor; 0.293137*178.7 in Humboldtv3+earlier

        if(file.Contains("TTJets_HT")){//LO cross sections with k-factor of 1.625 already applied
          if(file.Contains("2500toInf")) xsec = 0.0015; // 0.0023234211;TTJets_HT-1200to2500
          if(file.Contains("1200to2500")) xsec = 0.1326; // 0.194972521;
          if(file.Contains("800to1200")) xsec = 0.7378; // 1.07722318;
          if(file.Contains("600to800")) xsec = 1.7966; // 2.61537118;
        }

        // from cross XSDB
        if(file.Contains("TTG")) xsec = 3.697; //from twiki summary table; 4.078 in Humboldtv3+earlier               
        if(file.Contains("TTTT_Tune")) xsec = 0.008213; //XSDB/GenXSecAnalyzer
        
        // From https://cms-pdmv.cern.ch/mcm
        // k-factors from https://mangano.web.cern.ch/mangano/public/MECCA/samples_50ns_miniaod.txt
        // k-factors are ratio of https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
        // NLO/NNLO cross-sections to that of an inclusive sample in mcm at lower order (LO/NLO)

        if(file.Contains("WJetsToLNu_TuneCP5")) xsec=20508.9*3; //NNLO from Lesya's summary table 

        //cross-section per slice based on inclusive sample, roughly 10% higher than 2016, less in extreme tail
        //for corrections, see https://twiki.cern.ch/twiki/bin/viewauth/CMS/MCKnownIssues#WJetsToLNu_HT_and_DYJets_HT_LO_M
        if(file.Contains("WJetsToLNu_HT-70To100_TuneCP5")) {
          if (year==2017)  xsec = 1292.0*1.162; //XSDB * XSDB K-Factor 0.0243795*20508.9*3 in Humboldtv3+earlier
          else  xsec = 1292.0*1.164; //GenXSecAnalyzer * XSDB K-Factor 0.0243795*20508.9*3 in Humboldtv3+earlier
        }
        if(file.Contains("WJetsToLNu_HT-100To200_TuneCP5")) {
          if (year==2017)  xsec = 1395.0*1.162*0.993; //XSDB * XSDB K-Factor * Correction 0.0262096*20508.9*3 in Humboldtv3+earlier
          else  xsec = 1393.0*1.164*0.993; //GenXSecAnalyzer * XSDB K-Factor * Correction 0.0262096*20508.9*3 in Humboldtv3+earlier
        }
        if(file.Contains("WJetsToLNu_HT-200To400_TuneCP5")) {
          if (year==2017)  xsec = 407.9*1.162*1.002; //XSDB * XSDB K-Factor * Correction 0.00772818*20508.9*3 in Humboldtv3+earlier
          else  xsec = 409.9*1.164*1.002; //GenXSecAnalyzer * XSDB K-Factor * Correction 0.00772818*20508.9*3 in Humboldtv3+earlier
        }
        if(file.Contains("WJetsToLNu_HT-400To600_TuneCP5")) {
          if (year==2017)  xsec = 57.48*1.162*1.009; //XSDB * XSDB K-Factor * Correction 0.00109366*20508.9*3 in Humboldtv3+earlier
          else  xsec = 57.80*1.164*1.009; //GenXSecAnalyzer * XSDB K-Factor * Correction 0.00109366*20508.9*3 in Humboldtv3+earlier
        }
        if(file.Contains("WJetsToLNu_HT-600To800_TuneCP5")) {
          if (year==2017)  xsec = 12.87*1.162*1.120; //XSDB * XSDB K-Factor * Correction 0.000272388*20508.9*3 in Humboldtv3+earlier
          else  xsec = 12.94*1.164*1.120; //GenXSecAnalyzer * XSDB K-Factor * Correction 0.000272388*20508.9*3 in Humboldtv3+earlier
        }
        if(file.Contains("WJetsToLNu_HT-800To1200_TuneCP5")) {
          if (year==2017)  xsec = 5.366*1.162*1.202; //XSDB * XSDB K-Factor * Correction 0.000122233*20508.9*3 in Humboldtv3+earlier
          else  xsec = 5.451*1.164*1.202; //GenXSecAnalyzer * XSDB K-Factor * Correction 0.000122233*20508.9*3 in Humboldtv3+earlier
        }
        if(file.Contains("WJetsToLNu_HT-1200To2500_TuneCP5")) {
          if (year==2017)  xsec = 1.074*1.162*1.332; //XSDB * XSDB K-Factor * Correction 2.71060e-5*20508.9*3 in Humboldtv3+earlier
          else  xsec = 1.085*1.164*1.332; //GenXSecAnalyzer * XSDB K-Factor * Correction 2.71060e-5*20508.9*3 in Humboldtv3+earlier
        }
        if(file.Contains("WJetsToLNu_HT-2500ToInf_TuneCP5")) {
          if (year==2017)  xsec = 0.008001*1.162*4.200; //XSDB * XSDB K-Factor * Correction 3.94174e-07*20508.9*3 in Humboldtv3+earlier
          else  xsec = 0.008060*1.164*4.200; //GenXSecAnalyzer * XSDB K-Factor * Correction 3.94174e-07*20508.9*3 in Humboldtv3+earlier
        }

        if(file.Contains("QCD_HT100to200_TuneCP5")) xsec = 23700000;
        if(file.Contains("QCD_HT200to300_TuneCP5")) {
          if (year==2017) xsec = 1546000; //XSDB (LO)
          else xsec=1557000; //XSDB (LO), 1546000 in Humboldtv3+earlier
        }        
        if(file.Contains("QCD_HT300to500_TuneCP5")) {
          if (year==2017) xsec = 322600; //XSDB (LO)
          else xsec=323400; //XSDB (LO), 323400 in Humboldtv3+earlier
        }        
        if(file.Contains("QCD_HT500to700_TuneCP5")) {
          if (year==2017) xsec = 29980; //XSDB (LO)
          else xsec=30140; //XSDB (LO), 30140 in Humboldtv3+earlier
        }        
        if(file.Contains("QCD_HT700to1000_TuneCP5")) {
          if (year==2017) xsec = 6334; //XSDB (LO)
          else xsec=6310; //XSDB (LO), 6310 in Humboldtv3+earlier
        }        
        if(file.Contains("QCD_HT1000to1500_TuneCP5")) {
          if (year==2017) xsec = 1088; //XSDB (LO)
          else xsec=1094; //XSDB (LO), 1094 in Humboldtv3+earlier
        }        
        if(file.Contains("QCD_HT1500to2000_TuneCP5")) {
          if (year==2017) xsec = 99.11; //XSDB (LO)
          else xsec=99.38; //XSDB (LO), 99.38 in Humboldtv3+earlier
        }        
        if(file.Contains("QCD_HT2000toInf_TuneCP5")) {
          if (year==2017) xsec = 20.23; //XSDB (LO)
          else xsec=20.20; //XSDB (LO), 20.23 in Humboldtv3+earlier
        }        

        // Cross sections from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec
        // multiplied by BF(W->mu,e,tau) = 0.324
        if (file.Contains("ST_s-channel_4f_leptonDecays_TuneCP5"))     xsec = 3.74; //GenXSecAnalyzer, 3.34 in Humboldtv3+earlier
        if (file.Contains("ST_t-channel_antitop_4f_inclusiveDecays") ||
            file.Contains("ST_t-channel_antitop_4f_InclusiveDecays")) {
          if (year == 2017) xsec = 67.91; //XSDB (NLO) ; 80.95 in Humboldtv3+earlier
          else xsec = 69.09; //GenXSecAnalyzer ; 80.95 in Humboldtv3+earlier
        }
        if (file.Contains("ST_t-channel_top_4f_inclusiveDecays") || 
            file.Contains("ST_t-channel_top_4f_InclusiveDecays")) {
          if (year == 2017) xsec = 113.3; //XSDB (NLO) ; 136.02 in Humboldtv3+earlier
          else xsec = 115.3; //GenXSecAnalyzer ; 136.02 in Humboldtv3+earlier
        }
        if (file.Contains("ST_tW_antitop_5f_NoFullyHadronicDecays")) xsec = 34.97*0.530775; //XSDB/GenXSecAnalyzer * PDG BR; 35.85*0.543 in Humboldtv3+earlier
        if (file.Contains("ST_tW_top_5f_NoFullyHadronicDecays"))     xsec = 34.91*0.530775; //XSDB/GenXSecAnalyzer * PDG BR; 35.85*0.543 in Humboldtv3+earlier

        if(file.Contains("DYJetsToLL_M-50_TuneCP5"))               xsec = 6077.22;

        if(file.Contains("DYJetsToLL_M-50_HT-70to100_TuneCP5"))    xsec = 0.0275032*6077.22;
        if(file.Contains("DYJetsToLL_M-50_HT-100to200_TuneCP5"))   xsec = 0.0302083*6077.22;
        if(file.Contains("DYJetsToLL_M-50_HT-200to400_TuneCP5"))   xsec = 0.00907651*6077.22;
        if(file.Contains("DYJetsToLL_M-50_HT-400to600_TuneCP5"))   xsec = 0.00129238*6077.22;
        if(file.Contains("DYJetsToLL_M-50_HT-600to800_TuneCP5"))   xsec = 0.000316039*6077.22;
        if(file.Contains("DYJetsToLL_M-50_HT-800to1200_TuneCP5"))  xsec = 0.000137432*6077.22;
        if(file.Contains("DYJetsToLL_M-50_HT-1200to2500_TuneCP5")) xsec = 3.09368e-05*6077.22;
        if(file.Contains("DYJetsToLL_M-50_HT-2500toInf_TuneCP5"))  xsec = 4.39860e-07*6077.22;

        // k-factor from DYJets 1.165
        if(file.Contains("ZJetsToNuNu_HT-100To200")) {
          if (year==2017)  xsec = 302.8*1.1374*0.994; //XSDB * XSDB K-Factor * Correction 302.8*1.165 in Humboldtv3+earlier
          else  xsec = 304.0*1.1421*0.994; //GenXSecAnalyzer * XSDB K-Factor * Correction 302.8*1.165 in Humboldtv3+earlier
        }
        if(file.Contains("ZJetsToNuNu_HT-200To400")) {
          if (year==2017)  xsec = 92.59*1.1374*0.981; //XSDB * XSDB K-Factor * Correction 92.59*1.165 in Humboldtv3+earlier
          else  xsec = 91.68*1.1421*0.981; //GenXSecAnalyzer * XSDB K-Factor * Correction 92.59*1.165 in Humboldtv3+earlier
        }
        if(file.Contains("ZJetsToNuNu_HT-400To600")) {
          if (year==2017)  xsec = 13.18*1.1374*0.977; //XSDB * XSDB K-Factor * Correction 13.18*1.165 in Humboldtv3+earlier
          else  xsec = 13.11*1.1421*0.977; //GenXSecAnalyzer * XSDB K-Factor * Correction 13.18*1.165 in Humboldtv3+earlier
        }
        if(file.Contains("ZJetsToNuNu_HT-600To800")) {
          if (year==2017)  xsec = 3.257*1.1374*0.975; //XSDB * XSDB K-Factor * Correction 3.257*1.165 in Humboldtv3+earlier
          else  xsec = 3.245*1.1421*0.975; //GenXSecAnalyzer * XSDB K-Factor * Correction 3.257*1.165 in Humboldtv3+earlier
        }
        if(file.Contains("ZJetsToNuNu_HT-800To1200")) {
          if (year==2017)  xsec = 1.49*1.1374*0.916; //XSDB * XSDB K-Factor * Correction 1.49*1.165 in Humboldtv3+earlier
          else  xsec = 1.497*1.1421*0.916; //GenXSecAnalyzer * XSDB K-Factor * Correction 1.49*1.165 in Humboldtv3+earlier
        }
        if(file.Contains("ZJetsToNuNu_HT-1200To2500")) {
          if (year==2017)  xsec = 0.3419*1.1374*0.880; //XSDB * XSDB K-Factor * Correction 0.3419*1.165 in Humboldtv3+earlier
          else  xsec = 0.3425*1.1421*0.880; //GenXSecAnalyzer * XSDB K-Factor * Correction 0.3419*1.165 in Humboldtv3+earlier
        }
        if(file.Contains("ZJetsToNuNu_HT-2500ToInf")) {
          if (year==2017)  xsec = 0.005146*1.1374*1.276; //XSDB * XSDB K-Factor * Correction 0.3419*1.165 in Humboldtv3+earlier
          else  xsec = 0.005263*1.1421*1.276; //GenXSecAnalyzer * XSDB K-Factor * Correction 0.005146*1.165 in Humboldtv3+earlier
        }

        if(file.Contains("TTZToQQ_TuneCP5"))           xsec = 0.5104; //XSDB/GenXSecAnalyzer
        if(file.Contains("TTZToLLNuNu_M-10_TuneCP5"))  xsec = 0.2432; //XSDB/GenXSecAnalyzer
        if(file.Contains("TTWJetsToQQ_TuneCP5"))       xsec = 0.4062; //twiki table; 0.4316 in Humboldtv3+earlier
        if(file.Contains("TTWJetsToLNu_TuneCP5"))      xsec = 0.2043; //twiki table; 0.2149 in Humboldtv3+earlier
       
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
        if(file.Contains("WWTo2L2Nu_NNPDF31_TuneCP5")) xsec = 12.178; //NNLO
        if(file.Contains("WWToLNuQQ_NNPDF31_TuneCP5")) xsec = 49.997; //NNLO RA2b uses WWTo1L1Nu2Q

        if(file.Contains("WZTo1L3Nu")) {
          if (year==2017)  xsec = 3.294; //XSDB(LO) ; 3.05 in Humboldtv3+earlier
          else  xsec = 3.322; //GenXSecAnalyzer ; 3.05 in Humboldtv3+earlier
        }
        if(file.Contains("WZTo1L1Nu2Q")) {
          if (year==2017)  xsec = 11.66; //XSDB(LO) ; 10.73 in Humboldtv3+earlier
          else  xsec = 11.76; //GenXSecAnalyzer ; 10.73 in Humboldtv3+earlier
        }
        if(file.Contains("WZTo2L2Q"))          xsec = 5.606;
        if(file.Contains("WZTo3LNu_TuneCP5"))  xsec = 4.42965;
        if(file.Contains("ZZ_TuneCP5"))        xsec = 16.91;

        // Calculated at 13 TeV in
        // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV
        // Higgs branching ratios from
        // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR
        if(file.Contains("ZH_HToBB_ZToNuNu_M125"))    xsec = 0.883*Htobb*0.2;
        if(file.Contains("WHiggs0PMToBB"))     xsec = 1.373*Htobb*(0.1071+0.1063+0.1138); //using same as 2016

    }
    if(file.Contains("ttHTobb_M125")) xsec = 0.2934;

    // Zgamma cross sections at 13 TeV
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
    if(file.Contains("EWKZ2Jets"))            xsec = 6.215; // from XSDB
    if(file.Contains("ZGTo2LG"))              xsec = 117.864; // 123.8 in XSDB
    if(file.Contains("ZGToLLG"))              xsec = 55.48; // Other sample name
    if(file.Contains("TTTo2L2Nu"))            xsec = 88.29;
    if(file.Contains("TGJets") && !file.Contains("TTGJets")) xsec = 2.967;

    if(file.Contains("LLAJJ"))                xsec = 0.1084; // from XSDB
    if(file.Contains("WZ_Tune"))              xsec = 51.11;
    if(file.Contains("WW_Tune"))              xsec = 118.7;
    if(file.Contains("WWG_Tune"))             xsec = 0.2147;
    if(file.Contains("WZG_Tune"))             xsec = 0.04345; // from XSDB
    
    if(file.Contains("GluGluHToGG"))          xsec = HToGG * 48.58 ;
    if(file.Contains("GluGluHToTauTau"))      xsec = HToTT * 48.58 ;
    if(file.Contains("GluGluHToWWTo2L2Nu"))   xsec = HToWW * WToLNu * WToLNu * 48.58 ;

    if(file.Contains("GluGluHToMuMu"))        xsec = HToMM * 48.58 ;
    if(file.Contains("VBFHToMuMu"))           xsec = HToMM * 3.782 ;
    if(file.Contains("WminusH_HToMuMu"))      xsec = HToMM * 0.527 ;
    if(file.Contains("WplusH_HToMuMu"))       xsec = HToMM * 0.831 ;
    if(file.Contains("ZH_HToMuMu"))           xsec = HToMM * 0.8839 ;
    if(file.Contains("ttHToMuMu"))            xsec = HToMM * 0.5071 ;
    
    if(file.Contains("GluGluHToZZTo2L2Nu"))   xsec = HToZZ * ZToLL * ZToNuNu * 48.58 ;
    if(file.Contains("GluGluHToZZTo2L2Q"))    xsec = HToZZ * ZToLL * ZToQQ * 48.58 ;
    if(file.Contains("GluGluHToZZTo4L"))      xsec = HToZZ * ZToLL * ZToLL * 48.58 ;
    
    // Zgamma signal
    if(file.Contains("GluGluHToZG"))          xsec = HToZG * ZToLL * 48.58 ;
    if(file.Contains("VBFHToZG"))             xsec = HToZG * ZToLL * 3.782 ;
    if(file.Contains("WplusH_HToZG"))         xsec = HToZG * ZToLL * 0.831 ;
    if(file.Contains("WminusH_HToZG"))        xsec = HToZG * ZToLL * 0.527 ;
    if(file.Contains("ZH_HToZG"))             xsec = HToZG * ZToLL * 0.8839;
    if(file.Contains("ZH_ZToAll_HToZG"))      xsec = HToZG * ZToLL * 0.8839;
    if(file.Contains("ttHToZG"))              xsec = HToZG * ZToLL * 0.5071;
    
    if(xsec<=0) std::cout<<"ERROR:: Cross section not found for "<<file<<std::endl;

    return xsec;
  }

  float fractionNegWeights(const TString &file){
    float fneg(0.);

    // ttH, ttZ, ttW, ttGamma
    if(file.Contains("ttHJetTobb_M125_13TeV_amcatnloFXFX"))     fneg = 0.3515;
    if(file.Contains("TTZToQQ"))                                fneg = 0.2657;
    if(file.Contains("TTZToLLNuNu_M-10"))                       fneg = 0.2676;
    if(file.Contains("TTWJetsToQQ"))                            fneg = 0.2412;
    if(file.Contains("TTWJetsToLNu"))                           fneg = 0.2433;
    if(file.Contains("TTG"))                                    fneg = 0.342; // from MCM

    if(file.Contains("TTTT_TuneCUETP8M1_13TeV-amcatnlo"))       fneg = 0.41; // from MCM

    // Single top
    if (file.Contains("ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8")) fneg = 0.1884;

    return fneg;
  }

  // xsec is pb, xsec_unc is relative uncertainty
  void gluinoCrossSection(int glu_mass, double &xsec, double &xsec_unc) {
    if (glu_mass == 500) { xsec = 33.8; xsec_unc = 0.07869822485207102; return;}
    else if (glu_mass == 505) { xsec = 31.9; xsec_unc = 0.07899686520376176; return;}
    else if (glu_mass == 510) { xsec = 30.1; xsec_unc = 0.07906976744186046; return;}
    else if (glu_mass == 515) { xsec = 28.4; xsec_unc = 0.07957746478873239; return;}
    else if (glu_mass == 520) { xsec = 26.8; xsec_unc = 0.07985074626865672; return;}
    else if (glu_mass == 525) { xsec = 25.3; xsec_unc = 0.08023715415019762; return;}
    else if (glu_mass == 530) { xsec = 24.0; xsec_unc = 0.08083333333333333; return;}
    else if (glu_mass == 535) { xsec = 22.7; xsec_unc = 0.08105726872246696; return;}
    else if (glu_mass == 540) { xsec = 21.4; xsec_unc = 0.08130841121495327; return;}
    else if (glu_mass == 545) { xsec = 20.3; xsec_unc = 0.08177339901477831; return;}
    else if (glu_mass == 550) { xsec = 19.2; xsec_unc = 0.08229166666666668; return;}
    else if (glu_mass == 555) { xsec = 18.2; xsec_unc = 0.08241758241758242; return;}
    else if (glu_mass == 560) { xsec = 17.2; xsec_unc = 0.08313953488372093; return;}
    else if (glu_mass == 565) { xsec = 16.3; xsec_unc = 0.0834355828220859; return;}
    else if (glu_mass == 570) { xsec = 15.5; xsec_unc = 0.08387096774193549; return;}
    else if (glu_mass == 575) { xsec = 14.7; xsec_unc = 0.08435374149659865; return;}
    else if (glu_mass == 580) { xsec = 13.9; xsec_unc = 0.0841726618705036; return;}
    else if (glu_mass == 585) { xsec = 13.2; xsec_unc = 0.08484848484848487; return;}
    else if (glu_mass == 590) { xsec = 12.6; xsec_unc = 0.08492063492063492; return;}
    else if (glu_mass == 595) { xsec = 11.9; xsec_unc = 0.08571428571428572; return;}
    else if (glu_mass == 600) { xsec = 11.3; xsec_unc = 0.08619469026548672; return;}
    else if (glu_mass == 605) { xsec = 10.8; xsec_unc = 0.08657407407407407; return;}
    else if (glu_mass == 610) { xsec = 10.2; xsec_unc = 0.0869607843137255; return;}
    else if (glu_mass == 615) { xsec = 9.74; xsec_unc = 0.08737166324435318; return;}
    else if (glu_mass == 620) { xsec = 9.26; xsec_unc = 0.08779697624190064; return;}
    else if (glu_mass == 625) { xsec = 8.81; xsec_unc = 0.08819523269012486; return;}
    else if (glu_mass == 630) { xsec = 8.39; xsec_unc = 0.08855780691299164; return;}
    else if (glu_mass == 635) { xsec = 7.99; xsec_unc = 0.08898623279098873; return;}
    else if (glu_mass == 640) { xsec = 7.61; xsec_unc = 0.08935611038107753; return;}
    else if (glu_mass == 645) { xsec = 7.25; xsec_unc = 0.08979310344827586; return;}
    else if (glu_mass == 650) { xsec = 6.9; xsec_unc = 0.09014492753623188; return;}
    else if (glu_mass == 655) { xsec = 6.58; xsec_unc = 0.09057750759878419; return;}
    else if (glu_mass == 660) { xsec = 6.27; xsec_unc = 0.09106858054226476; return;}
    else if (glu_mass == 665) { xsec = 5.98; xsec_unc = 0.09147157190635452; return;}
    else if (glu_mass == 670) { xsec = 5.71; xsec_unc = 0.09194395796847636; return;}
    else if (glu_mass == 675) { xsec = 5.44; xsec_unc = 0.09227941176470587; return;}
    else if (glu_mass == 680) { xsec = 5.2; xsec_unc = 0.09269230769230768; return;}
    else if (glu_mass == 685) { xsec = 4.96; xsec_unc = 0.09314516129032259; return;}
    else if (glu_mass == 690) { xsec = 4.74; xsec_unc = 0.09345991561181434; return;}
    else if (glu_mass == 695) { xsec = 4.52; xsec_unc = 0.09380530973451329; return;}
    else if (glu_mass == 700) { xsec = 4.32; xsec_unc = 0.09421296296296296; return;}
    else if (glu_mass == 705) { xsec = 4.13; xsec_unc = 0.09467312348668282; return;}
    else if (glu_mass == 710) { xsec = 3.95; xsec_unc = 0.09493670886075949; return;}
    else if (glu_mass == 715) { xsec = 3.77; xsec_unc = 0.09549071618037135; return;}
    else if (glu_mass == 720) { xsec = 3.61; xsec_unc = 0.09584487534626038; return;}
    else if (glu_mass == 725) { xsec = 3.45; xsec_unc = 0.09623188405797102; return;}
    else if (glu_mass == 730) { xsec = 3.3; xsec_unc = 0.09636363636363637; return;}
    else if (glu_mass == 735) { xsec = 3.16; xsec_unc = 0.09683544303797467; return;}
    else if (glu_mass == 740) { xsec = 3.02; xsec_unc = 0.09735099337748344; return;}
    else if (glu_mass == 745) { xsec = 2.89; xsec_unc = 0.09757785467128026; return;}
    else if (glu_mass == 750) { xsec = 2.77; xsec_unc = 0.09783393501805054; return;}
    else if (glu_mass == 755) { xsec = 2.65; xsec_unc = 0.09811320754716982; return;}
    else if (glu_mass == 760) { xsec = 2.54; xsec_unc = 0.09881889763779528; return;}
    else if (glu_mass == 765) { xsec = 2.43; xsec_unc = 0.09917695473251027; return;}
    else if (glu_mass == 770) { xsec = 2.33; xsec_unc = 0.09957081545064378; return;}
    else if (glu_mass == 775) { xsec = 2.23; xsec_unc = 0.1; return;}
    else if (glu_mass == 780) { xsec = 2.14; xsec_unc = 0.09999999999999999; return;}
    else if (glu_mass == 785) { xsec = 2.05; xsec_unc = 0.10048780487804879; return;}
    else if (glu_mass == 790) { xsec = 1.97; xsec_unc = 0.10101522842639594; return;}
    else if (glu_mass == 795) { xsec = 1.88; xsec_unc = 0.10106382978723405; return;}
    else if (glu_mass == 800) { xsec = 1.81; xsec_unc = 0.1016574585635359; return;}
    else if (glu_mass == 805) { xsec = 1.73; xsec_unc = 0.10173410404624277; return;}
    else if (glu_mass == 810) { xsec = 1.66; xsec_unc = 0.10240963855421688; return;}
    else if (glu_mass == 815) { xsec = 1.6; xsec_unc = 0.1025; return;}
    else if (glu_mass == 820) { xsec = 1.53; xsec_unc = 0.10326797385620914; return;}
    else if (glu_mass == 825) { xsec = 1.47; xsec_unc = 0.10340136054421768; return;}
    else if (glu_mass == 830) { xsec = 1.41; xsec_unc = 0.10354609929078014; return;}
    else if (glu_mass == 835) { xsec = 1.36; xsec_unc = 0.10441176470588234; return;}
    else if (glu_mass == 840) { xsec = 1.3; xsec_unc = 0.10461538461538462; return;}
    else if (glu_mass == 845) { xsec = 1.25; xsec_unc = 0.1048; return;}
    else if (glu_mass == 850) { xsec = 1.2; xsec_unc = 0.10500000000000001; return;}
    else if (glu_mass == 855) { xsec = 1.15; xsec_unc = 0.10608695652173913; return;}
    else if (glu_mass == 860) { xsec = 1.11; xsec_unc = 0.10630630630630629; return;}
    else if (glu_mass == 865) { xsec = 1.07; xsec_unc = 0.10654205607476636; return;}
    else if (glu_mass == 870) { xsec = 1.03; xsec_unc = 0.10679611650485436; return;}
    else if (glu_mass == 875) { xsec = 0.986; xsec_unc = 0.1075050709939148; return;}
    else if (glu_mass == 880) { xsec = 0.948; xsec_unc = 0.10759493670886076; return;}
    else if (glu_mass == 885) { xsec = 0.912; xsec_unc = 0.10789473684210527; return;}
    else if (glu_mass == 890) { xsec = 0.877; xsec_unc = 0.10820980615735462; return;}
    else if (glu_mass == 895) { xsec = 0.844; xsec_unc = 0.10864928909952608; return;}
    else if (glu_mass == 900) { xsec = 0.812; xsec_unc = 0.10886699507389162; return;}
    else if (glu_mass == 905) { xsec = 0.781; xsec_unc = 0.10934699103713189; return;}
    else if (glu_mass == 910) { xsec = 0.752; xsec_unc = 0.10970744680851065; return;}
    else if (glu_mass == 915) { xsec = 0.723; xsec_unc = 0.1099585062240664; return;}
    else if (glu_mass == 920) { xsec = 0.696; xsec_unc = 0.1103448275862069; return;}
    else if (glu_mass == 925) { xsec = 0.67; xsec_unc = 0.11074626865671641; return;}
    else if (glu_mass == 930) { xsec = 0.646; xsec_unc = 0.11114551083591331; return;}
    else if (glu_mass == 935) { xsec = 0.622; xsec_unc = 0.11141479099678457; return;}
    else if (glu_mass == 940) { xsec = 0.599; xsec_unc = 0.11185308848080135; return;}
    else if (glu_mass == 945) { xsec = 0.577; xsec_unc = 0.1121317157712305; return;}
    else if (glu_mass == 950) { xsec = 0.556; xsec_unc = 0.11258992805755395; return;}
    else if (glu_mass == 955) { xsec = 0.535; xsec_unc = 0.11271028037383177; return;}
    else if (glu_mass == 960) { xsec = 0.516; xsec_unc = 0.11317829457364341; return;}
    else if (glu_mass == 965) { xsec = 0.497; xsec_unc = 0.11348088531187123; return;}
    else if (glu_mass == 970) { xsec = 0.479; xsec_unc = 0.1139874739039666; return;}
    else if (glu_mass == 975) { xsec = 0.462; xsec_unc = 0.11428571428571428; return;}
    else if (glu_mass == 980) { xsec = 0.445; xsec_unc = 0.1146067415730337; return;}
    else if (glu_mass == 985) { xsec = 0.43; xsec_unc = 0.11488372093023255; return;}
    else if (glu_mass == 990) { xsec = 0.414; xsec_unc = 0.11521739130434783; return;}
    else if (glu_mass == 995) { xsec = 0.399; xsec_unc = 0.11553884711779448; return;}
    else if (glu_mass == 1000) { xsec = 0.385; xsec_unc = 0.1161038961038961; return;}
    else if (glu_mass == 1005) { xsec = 0.372; xsec_unc = 0.11639784946236559; return;}
    else if (glu_mass == 1010) { xsec = 0.359; xsec_unc = 0.11671309192200557; return;}
    else if (glu_mass == 1015) { xsec = 0.346; xsec_unc = 0.1170520231213873; return;}
    else if (glu_mass == 1020) { xsec = 0.334; xsec_unc = 0.11736526946107784; return;}
    else if (glu_mass == 1025) { xsec = 0.322; xsec_unc = 0.11770186335403728; return;}
    else if (glu_mass == 1030) { xsec = 0.311; xsec_unc = 0.11800643086816721; return;}
    else if (glu_mass == 1035) { xsec = 0.3; xsec_unc = 0.11866666666666667; return;}
    else if (glu_mass == 1040) { xsec = 0.29; xsec_unc = 0.11896551724137933; return;}
    else if (glu_mass == 1045) { xsec = 0.28; xsec_unc = 0.11928571428571427; return;}
    else if (glu_mass == 1050) { xsec = 0.27; xsec_unc = 0.11962962962962963; return;}
    else if (glu_mass == 1055) { xsec = 0.261; xsec_unc = 0.11992337164750957; return;}
    else if (glu_mass == 1060) { xsec = 0.252; xsec_unc = 0.12023809523809524; return;}
    else if (glu_mass == 1065) { xsec = 0.243; xsec_unc = 0.1205761316872428; return;}
    else if (glu_mass == 1070) { xsec = 0.235; xsec_unc = 0.12085106382978725; return;}
    else if (glu_mass == 1075) { xsec = 0.227; xsec_unc = 0.1211453744493392; return;}
    else if (glu_mass == 1080) { xsec = 0.219; xsec_unc = 0.1219178082191781; return;}
    else if (glu_mass == 1085) { xsec = 0.212; xsec_unc = 0.12216981132075472; return;}
    else if (glu_mass == 1090) { xsec = 0.205; xsec_unc = 0.12243902439024391; return;}
    else if (glu_mass == 1095) { xsec = 0.198; xsec_unc = 0.12272727272727271; return;}
    else if (glu_mass == 1100) { xsec = 0.191; xsec_unc = 0.12303664921465969; return;}
    else if (glu_mass == 1105) { xsec = 0.185; xsec_unc = 0.12324324324324325; return;}
    else if (glu_mass == 1110) { xsec = 0.179; xsec_unc = 0.12402234636871509; return;}
    else if (glu_mass == 1115) { xsec = 0.173; xsec_unc = 0.12427745664739884; return;}
    else if (glu_mass == 1120) { xsec = 0.167; xsec_unc = 0.1245508982035928; return;}
    else if (glu_mass == 1125) { xsec = 0.162; xsec_unc = 0.12469135802469135; return;}
    else if (glu_mass == 1130) { xsec = 0.156; xsec_unc = 0.125; return;}
    else if (glu_mass == 1135) { xsec = 0.151; xsec_unc = 0.12582781456953643; return;}
    else if (glu_mass == 1140) { xsec = 0.146; xsec_unc = 0.12602739726027398; return;}
    else if (glu_mass == 1145) { xsec = 0.141; xsec_unc = 0.12624113475177307; return;}
    else if (glu_mass == 1150) { xsec = 0.137; xsec_unc = 0.12700729927007298; return;}
    else if (glu_mass == 1155) { xsec = 0.132; xsec_unc = 0.12727272727272726; return;}
    else if (glu_mass == 1160) { xsec = 0.128; xsec_unc = 0.12734374999999998; return;}
    else if (glu_mass == 1165) { xsec = 0.124; xsec_unc = 0.1274193548387097; return;}
    else if (glu_mass == 1170) { xsec = 0.12; xsec_unc = 0.12833333333333335; return;}
    else if (glu_mass == 1175) { xsec = 0.116; xsec_unc = 0.12844827586206897; return;}
    else if (glu_mass == 1180) { xsec = 0.112; xsec_unc = 0.12857142857142856; return;}
    else if (glu_mass == 1185) { xsec = 0.109; xsec_unc = 0.12935779816513762; return;}
    else if (glu_mass == 1190) { xsec = 0.105; xsec_unc = 0.1295238095238095; return;}
    else if (glu_mass == 1195) { xsec = 0.102; xsec_unc = 0.1303921568627451; return;}
    else if (glu_mass == 1200) { xsec = 0.0985; xsec_unc = 0.1299492385786802; return;}
    else if (glu_mass == 1205) { xsec = 0.0953; xsec_unc = 0.1311647429171039; return;}
    else if (glu_mass == 1210) { xsec = 0.0923; xsec_unc = 0.13109425785482123; return;}
    else if (glu_mass == 1215) { xsec = 0.0894; xsec_unc = 0.13087248322147652; return;}
    else if (glu_mass == 1220) { xsec = 0.0866; xsec_unc = 0.13163972286374134; return;}
    else if (glu_mass == 1225) { xsec = 0.0838; xsec_unc = 0.1324582338902148; return;}
    else if (glu_mass == 1230) { xsec = 0.0812; xsec_unc = 0.1330049261083744; return;}
    else if (glu_mass == 1235) { xsec = 0.0786; xsec_unc = 0.13231552162849872; return;}
    else if (glu_mass == 1240) { xsec = 0.0762; xsec_unc = 0.13254593175853016; return;}
    else if (glu_mass == 1245) { xsec = 0.0738; xsec_unc = 0.13333333333333333; return;}
    else if (glu_mass == 1250) { xsec = 0.0715; xsec_unc = 0.13384615384615386; return;}
    else if (glu_mass == 1255) { xsec = 0.0692; xsec_unc = 0.13410404624277458; return;}
    else if (glu_mass == 1260) { xsec = 0.0671; xsec_unc = 0.13442622950819672; return;}
    else if (glu_mass == 1265) { xsec = 0.065; xsec_unc = 0.13476923076923078; return;}
    else if (glu_mass == 1270) { xsec = 0.063; xsec_unc = 0.13507936507936508; return;}
    else if (glu_mass == 1275) { xsec = 0.061; xsec_unc = 0.13557377049180327; return;}
    else if (glu_mass == 1280) { xsec = 0.0591; xsec_unc = 0.1358714043993232; return;}
    else if (glu_mass == 1285) { xsec = 0.0573; xsec_unc = 0.13612565445026178; return;}
    else if (glu_mass == 1290) { xsec = 0.0556; xsec_unc = 0.1365107913669065; return;}
    else if (glu_mass == 1295) { xsec = 0.0539; xsec_unc = 0.13692022263450834; return;}
    else if (glu_mass == 1300) { xsec = 0.0522; xsec_unc = 0.13716475095785438; return;}
    else if (glu_mass == 1305) { xsec = 0.0506; xsec_unc = 0.13754940711462452; return;}
    else if (glu_mass == 1310) { xsec = 0.0491; xsec_unc = 0.13788187372708757; return;}
    else if (glu_mass == 1315) { xsec = 0.0476; xsec_unc = 0.13823529411764704; return;}
    else if (glu_mass == 1320) { xsec = 0.0461; xsec_unc = 0.1386117136659436; return;}
    else if (glu_mass == 1325) { xsec = 0.0447; xsec_unc = 0.13892617449664432; return;}
    else if (glu_mass == 1330) { xsec = 0.0434; xsec_unc = 0.13940092165898615; return;}
    else if (glu_mass == 1335) { xsec = 0.0421; xsec_unc = 0.13966745843230405; return;}
    else if (glu_mass == 1340) { xsec = 0.0408; xsec_unc = 0.14019607843137255; return;}
    else if (glu_mass == 1345) { xsec = 0.0396; xsec_unc = 0.1404040404040404; return;}
    else if (glu_mass == 1350) { xsec = 0.0384; xsec_unc = 0.14088541666666668; return;}
    else if (glu_mass == 1355) { xsec = 0.0372; xsec_unc = 0.14112903225806453; return;}
    else if (glu_mass == 1360) { xsec = 0.0361; xsec_unc = 0.14155124653739612; return;}
    else if (glu_mass == 1365) { xsec = 0.035; xsec_unc = 0.142; return;}
    else if (glu_mass == 1370) { xsec = 0.034; xsec_unc = 0.14205882352941177; return;}
    else if (glu_mass == 1375) { xsec = 0.033; xsec_unc = 0.1427272727272727; return;}
    else if (glu_mass == 1380) { xsec = 0.032; xsec_unc = 0.143125; return;}
    else if (glu_mass == 1385) { xsec = 0.031; xsec_unc = 0.14322580645161292; return;}
    else if (glu_mass == 1390) { xsec = 0.0301; xsec_unc = 0.14385382059800664; return;}
    else if (glu_mass == 1395) { xsec = 0.0292; xsec_unc = 0.14383561643835616; return;}
    else if (glu_mass == 1400) { xsec = 0.0284; xsec_unc = 0.1443661971830986; return;}
    else if (glu_mass == 1405) { xsec = 0.0275; xsec_unc = 0.14472727272727273; return;}
    else if (glu_mass == 1410) { xsec = 0.0267; xsec_unc = 0.1449438202247191; return;}
    else if (glu_mass == 1415) { xsec = 0.0259; xsec_unc = 0.14555984555984555; return;}
    else if (glu_mass == 1420) { xsec = 0.0252; xsec_unc = 0.14603174603174604; return;}
    else if (glu_mass == 1425) { xsec = 0.0244; xsec_unc = 0.14631147540983605; return;}
    else if (glu_mass == 1430) { xsec = 0.0237; xsec_unc = 0.14641350210970464; return;}
    else if (glu_mass == 1435) { xsec = 0.023; xsec_unc = 0.14695652173913046; return;}
    else if (glu_mass == 1440) { xsec = 0.0224; xsec_unc = 0.14732142857142858; return;}
    else if (glu_mass == 1445) { xsec = 0.0217; xsec_unc = 0.147926267281106; return;}
    else if (glu_mass == 1450) { xsec = 0.0211; xsec_unc = 0.14786729857819905; return;}
    else if (glu_mass == 1455) { xsec = 0.0205; xsec_unc = 0.14829268292682926; return;}
    else if (glu_mass == 1460) { xsec = 0.0199; xsec_unc = 0.1487437185929648; return;}
    else if (glu_mass == 1465) { xsec = 0.0193; xsec_unc = 0.14922279792746115; return;}
    else if (glu_mass == 1470) { xsec = 0.0187; xsec_unc = 0.1497326203208556; return;}
    else if (glu_mass == 1475) { xsec = 0.0182; xsec_unc = 0.15; return;}
    else if (glu_mass == 1480) { xsec = 0.0177; xsec_unc = 0.15028248587570622; return;}
    else if (glu_mass == 1485) { xsec = 0.0172; xsec_unc = 0.1505813953488372; return;}
    else if (glu_mass == 1490) { xsec = 0.0167; xsec_unc = 0.1508982035928144; return;}
    else if (glu_mass == 1495) { xsec = 0.0162; xsec_unc = 0.15123456790123457; return;}
    else if (glu_mass == 1500) { xsec = 0.0157; xsec_unc = 0.1515923566878981; return;}
    else if (glu_mass == 1505) { xsec = 0.0153; xsec_unc = 0.15228758169934642; return;}
    else if (glu_mass == 1510) { xsec = 0.0148; xsec_unc = 0.1527027027027027; return;}
    else if (glu_mass == 1515) { xsec = 0.0144; xsec_unc = 0.1527777777777778; return;}
    else if (glu_mass == 1520) { xsec = 0.014; xsec_unc = 0.15357142857142858; return;}
    else if (glu_mass == 1525) { xsec = 0.0136; xsec_unc = 0.1536764705882353; return;}
    else if (glu_mass == 1530) { xsec = 0.0132; xsec_unc = 0.1537878787878788; return;}
    else if (glu_mass == 1535) { xsec = 0.0128; xsec_unc = 0.1546875; return;}
    else if (glu_mass == 1540) { xsec = 0.0125; xsec_unc = 0.1552; return;}
    else if (glu_mass == 1545) { xsec = 0.0121; xsec_unc = 0.15537190082644628; return;}
    else if (glu_mass == 1550) { xsec = 0.0118; xsec_unc = 0.15593220338983052; return;}
    else if (glu_mass == 1555) { xsec = 0.0115; xsec_unc = 0.15565217391304348; return;}
    else if (glu_mass == 1560) { xsec = 0.0111; xsec_unc = 0.15675675675675674; return;}
    else if (glu_mass == 1565) { xsec = 0.0108; xsec_unc = 0.15648148148148147; return;}
    else if (glu_mass == 1570) { xsec = 0.0105; xsec_unc = 0.15714285714285714; return;}
    else if (glu_mass == 1575) { xsec = 0.0102; xsec_unc = 0.15784313725490196; return;}
    else if (glu_mass == 1580) { xsec = 0.00993; xsec_unc = 0.1581067472306143; return;}
    else if (glu_mass == 1585) { xsec = 0.00966; xsec_unc = 0.15838509316770186; return;}
    else if (glu_mass == 1590) { xsec = 0.00939; xsec_unc = 0.15867944621938232; return;}
    else if (glu_mass == 1595) { xsec = 0.00912; xsec_unc = 0.15899122807017543; return;}
    else if (glu_mass == 1600) { xsec = 0.00887; xsec_unc = 0.15896279594137544; return;}
    else if (glu_mass == 1605) { xsec = 0.00862; xsec_unc = 0.16009280742459397; return;}
    else if (glu_mass == 1610) { xsec = 0.00838; xsec_unc = 0.15990453460620524; return;}
    else if (glu_mass == 1615) { xsec = 0.00815; xsec_unc = 0.16073619631901842; return;}
    else if (glu_mass == 1620) { xsec = 0.00792; xsec_unc = 0.16161616161616163; return;}
    else if (glu_mass == 1625) { xsec = 0.0077; xsec_unc = 0.16103896103896104; return;}
    else if (glu_mass == 1630) { xsec = 0.00749; xsec_unc = 0.16154873164218958; return;}
    else if (glu_mass == 1635) { xsec = 0.00728; xsec_unc = 0.1620879120879121; return;}
    else if (glu_mass == 1640) { xsec = 0.00708; xsec_unc = 0.16242937853107345; return;}
    else if (glu_mass == 1645) { xsec = 0.00689; xsec_unc = 0.16255442670537007; return;}
    else if (glu_mass == 1650) { xsec = 0.0067; xsec_unc = 0.1626865671641791; return;}
    else if (glu_mass == 1655) { xsec = 0.00651; xsec_unc = 0.16436251920122888; return;}
    else if (glu_mass == 1660) { xsec = 0.00633; xsec_unc = 0.16429699842022116; return;}
    else if (glu_mass == 1665) { xsec = 0.00616; xsec_unc = 0.16396103896103897; return;}
    else if (glu_mass == 1670) { xsec = 0.00599; xsec_unc = 0.1649415692821369; return;}
    else if (glu_mass == 1675) { xsec = 0.00583; xsec_unc = 0.1653516295025729; return;}
    else if (glu_mass == 1680) { xsec = 0.00567; xsec_unc = 0.16560846560846562; return;}
    else if (glu_mass == 1685) { xsec = 0.00551; xsec_unc = 0.16606170598911071; return;}
    else if (glu_mass == 1690) { xsec = 0.00536; xsec_unc = 0.1664179104477612; return;}
    else if (glu_mass == 1695) { xsec = 0.00521; xsec_unc = 0.16679462571976966; return;}
    else if (glu_mass == 1700) { xsec = 0.00507; xsec_unc = 0.16725838264299803; return;}
    else if (glu_mass == 1705) { xsec = 0.00493; xsec_unc = 0.16754563894523325; return;}
    else if (glu_mass == 1710) { xsec = 0.0048; xsec_unc = 0.16812500000000002; return;}
    else if (glu_mass == 1715) { xsec = 0.00467; xsec_unc = 0.16852248394004285; return;}
    else if (glu_mass == 1720) { xsec = 0.00454; xsec_unc = 0.16894273127753304; return;}
    else if (glu_mass == 1725) { xsec = 0.00442; xsec_unc = 0.1692307692307692; return;}
    else if (glu_mass == 1730) { xsec = 0.0043; xsec_unc = 0.1697674418604651; return;}
    else if (glu_mass == 1735) { xsec = 0.00418; xsec_unc = 0.17009569377990433; return;}
    else if (glu_mass == 1740) { xsec = 0.00407; xsec_unc = 0.1705159705159705; return;}
    else if (glu_mass == 1745) { xsec = 0.00396; xsec_unc = 0.17095959595959595; return;}
    else if (glu_mass == 1750) { xsec = 0.00385; xsec_unc = 0.17142857142857143; return;}
    else if (glu_mass == 1755) { xsec = 0.00375; xsec_unc = 0.17173333333333335; return;}
    else if (glu_mass == 1760) { xsec = 0.00365; xsec_unc = 0.17232876712328768; return;}
    else if (glu_mass == 1765) { xsec = 0.00355; xsec_unc = 0.17267605633802818; return;}
    else if (glu_mass == 1770) { xsec = 0.00345; xsec_unc = 0.17304347826086958; return;}
    else if (glu_mass == 1775) { xsec = 0.00336; xsec_unc = 0.17351190476190476; return;}
    else if (glu_mass == 1780) { xsec = 0.00327; xsec_unc = 0.17400611620795106; return;}
    else if (glu_mass == 1785) { xsec = 0.00318; xsec_unc = 0.1742138364779874; return;}
    else if (glu_mass == 1790) { xsec = 0.0031; xsec_unc = 0.17483870967741935; return;}
    else if (glu_mass == 1795) { xsec = 0.00301; xsec_unc = 0.1750830564784053; return;}
    else if (glu_mass == 1800) { xsec = 0.00293; xsec_unc = 0.17576791808873723; return;}
    else if (glu_mass == 1805) { xsec = 0.00286; xsec_unc = 0.17587412587412585; return;}
    else if (glu_mass == 1810) { xsec = 0.00278; xsec_unc = 0.17625899280575538; return;}
    else if (glu_mass == 1815) { xsec = 0.00271; xsec_unc = 0.17675276752767527; return;}
    else if (glu_mass == 1820) { xsec = 0.00263; xsec_unc = 0.17718631178707225; return;}
    else if (glu_mass == 1825) { xsec = 0.00256; xsec_unc = 0.17773437499999997; return;}
    else if (glu_mass == 1830) { xsec = 0.00249; xsec_unc = 0.1783132530120482; return;}
    else if (glu_mass == 1835) { xsec = 0.00243; xsec_unc = 0.1786008230452675; return;}
    else if (glu_mass == 1840) { xsec = 0.00236; xsec_unc = 0.17881355932203388; return;}
    else if (glu_mass == 1845) { xsec = 0.0023; xsec_unc = 0.17956521739130435; return;}
    else if (glu_mass == 1850) { xsec = 0.00224; xsec_unc = 0.17991071428571428; return;}
    else if (glu_mass == 1855) { xsec = 0.00218; xsec_unc = 0.18027522935779816; return;}
    else if (glu_mass == 1860) { xsec = 0.00212; xsec_unc = 0.18066037735849055; return;}
    else if (glu_mass == 1865) { xsec = 0.00207; xsec_unc = 0.1811594202898551; return;}
    else if (glu_mass == 1870) { xsec = 0.00201; xsec_unc = 0.181592039800995; return;}
    else if (glu_mass == 1875) { xsec = 0.00196; xsec_unc = 0.18214285714285716; return;}
    else if (glu_mass == 1880) { xsec = 0.00191; xsec_unc = 0.18272251308900525; return;}
    else if (glu_mass == 1885) { xsec = 0.00186; xsec_unc = 0.18333333333333332; return;}
    else if (glu_mass == 1890) { xsec = 0.00181; xsec_unc = 0.18342541436464088; return;}
    else if (glu_mass == 1895) { xsec = 0.00176; xsec_unc = 0.18409090909090908; return;}
    else if (glu_mass == 1900) { xsec = 0.00171; xsec_unc = 0.1842105263157895; return;}
    else if (glu_mass == 1905) { xsec = 0.00167; xsec_unc = 0.18502994011976046; return;}
    else if (glu_mass == 1910) { xsec = 0.00163; xsec_unc = 0.18527607361963191; return;}
    else if (glu_mass == 1915) { xsec = 0.00158; xsec_unc = 0.1860759493670886; return;}
    else if (glu_mass == 1920) { xsec = 0.00154; xsec_unc = 0.18636363636363637; return;}
    else if (glu_mass == 1925) { xsec = 0.0015; xsec_unc = 0.18666666666666665; return;}
    else if (glu_mass == 1930) { xsec = 0.00146; xsec_unc = 0.18698630136986305; return;}
    else if (glu_mass == 1935) { xsec = 0.00142; xsec_unc = 0.18802816901408448; return;}
    else if (glu_mass == 1940) { xsec = 0.00139; xsec_unc = 0.18848920863309354; return;}
    else if (glu_mass == 1945) { xsec = 0.00135; xsec_unc = 0.18888888888888888; return;}
    else if (glu_mass == 1950) { xsec = 0.00131; xsec_unc = 0.18931297709923667; return;}
    else if (glu_mass == 1955) { xsec = 0.00128; xsec_unc = 0.18984374999999998; return;}
    else if (glu_mass == 1960) { xsec = 0.00125; xsec_unc = 0.1904; return;}
    else if (glu_mass == 1965) { xsec = 0.00121; xsec_unc = 0.19090909090909092; return;}
    else if (glu_mass == 1970) { xsec = 0.00118; xsec_unc = 0.19152542372881354; return;}
    else if (glu_mass == 1975) { xsec = 0.00115; xsec_unc = 0.19130434782608696; return;}
    else if (glu_mass == 1980) { xsec = 0.00112; xsec_unc = 0.19196428571428573; return;}
    else if (glu_mass == 1985) { xsec = 0.00109; xsec_unc = 0.1926605504587156; return;}
    else if (glu_mass == 1990) { xsec = 0.00106; xsec_unc = 0.19339622641509435; return;}
    else if (glu_mass == 1995) { xsec = 0.00104; xsec_unc = 0.1932692307692308; return;}
    else if (glu_mass == 2000) { xsec = 0.00101; xsec_unc = 0.19405940594059404; return;}
    else if (glu_mass == 2005) { xsec = 0.000983; xsec_unc = 0.19430315361139372; return;}
    else if (glu_mass == 2010) { xsec = 0.000957; xsec_unc = 0.19540229885057472; return;}
    else if (glu_mass == 2015) { xsec = 0.000933; xsec_unc = 0.19614147909967847; return;}
    else if (glu_mass == 2020) { xsec = 0.000908; xsec_unc = 0.1960352422907489; return;}
    else if (glu_mass == 2025) { xsec = 0.000885; xsec_unc = 0.1966101694915254; return;}
    else if (glu_mass == 2030) { xsec = 0.000862; xsec_unc = 0.19721577726218098; return;}
    else if (glu_mass == 2035) { xsec = 0.00084; xsec_unc = 0.1976190476190476; return;}
    else if (glu_mass == 2040) { xsec = 0.000818; xsec_unc = 0.1980440097799511; return;}
    else if (glu_mass == 2045) { xsec = 0.000797; xsec_unc = 0.19949811794228356; return;}
    else if (glu_mass == 2050) { xsec = 0.000776; xsec_unc = 0.19974226804123713; return;}
    else if (glu_mass == 2055) { xsec = 0.000756; xsec_unc = 0.19973544973544974; return;}
    else if (glu_mass == 2060) { xsec = 0.000737; xsec_unc = 0.20081411126187243; return;}
    else if (glu_mass == 2065) { xsec = 0.000718; xsec_unc = 0.201949860724234; return;}
    else if (glu_mass == 2070) { xsec = 0.000699; xsec_unc = 0.20171673819742492; return;}
    else if (glu_mass == 2075) { xsec = 0.000681; xsec_unc = 0.2026431718061674; return;}
    else if (glu_mass == 2080) { xsec = 0.000664; xsec_unc = 0.2033132530120482; return;}
    else if (glu_mass == 2085) { xsec = 0.000647; xsec_unc = 0.20401854714064915; return;}
    else if (glu_mass == 2090) { xsec = 0.00063; xsec_unc = 0.20476190476190473; return;}
    else if (glu_mass == 2095) { xsec = 0.000614; xsec_unc = 0.20521172638436483; return;}
    else if (glu_mass == 2100) { xsec = 0.000598; xsec_unc = 0.205685618729097; return;}
    else if (glu_mass == 2105) { xsec = 0.000583; xsec_unc = 0.2058319039451115; return;}
    else if (glu_mass == 2110) { xsec = 0.000568; xsec_unc = 0.20598591549295772; return;}
    else if (glu_mass == 2115) { xsec = 0.000553; xsec_unc = 0.20795660036166366; return;}
    else if (glu_mass == 2120) { xsec = 0.000539; xsec_unc = 0.2077922077922078; return;}
    else if (glu_mass == 2125) { xsec = 0.000525; xsec_unc = 0.20761904761904765; return;}
    else if (glu_mass == 2130) { xsec = 0.000512; xsec_unc = 0.208984375; return;}
    else if (glu_mass == 2135) { xsec = 0.000499; xsec_unc = 0.21042084168336675; return;}
    else if (glu_mass == 2140) { xsec = 0.000486; xsec_unc = 0.20987654320987653; return;}
    else if (glu_mass == 2145) { xsec = 0.000473; xsec_unc = 0.2109936575052854; return;}
    else if (glu_mass == 2150) { xsec = 0.000461; xsec_unc = 0.21149674620390455; return;}
    else if (glu_mass == 2155) { xsec = 0.000449; xsec_unc = 0.21224944320712694; return;}
    else if (glu_mass == 2160) { xsec = 0.000438; xsec_unc = 0.213013698630137; return;}
    else if (glu_mass == 2165) { xsec = 0.000427; xsec_unc = 0.21358313817330207; return;}
    else if (glu_mass == 2170) { xsec = 0.000416; xsec_unc = 0.21418269230769232; return;}
    else if (glu_mass == 2175) { xsec = 0.000405; xsec_unc = 0.21481481481481482; return;}
    else if (glu_mass == 2180) { xsec = 0.000395; xsec_unc = 0.21544303797468353; return;}
    else if (glu_mass == 2185) { xsec = 0.000385; xsec_unc = 0.21610389610389613; return;}
    else if (glu_mass == 2190) { xsec = 0.000375; xsec_unc = 0.2168; return;}
    else if (glu_mass == 2195) { xsec = 0.000365; xsec_unc = 0.21753424657534248; return;}
    else if (glu_mass == 2200) { xsec = 0.000356; xsec_unc = 0.21825842696629216; return;}
    else if (glu_mass == 2205) { xsec = 0.000347; xsec_unc = 0.21902017291066286; return;}
    else if (glu_mass == 2210) { xsec = 0.000338; xsec_unc = 0.21982248520710063; return;}
    else if (glu_mass == 2215) { xsec = 0.00033; xsec_unc = 0.22030303030303033; return;}
    else if (glu_mass == 2220) { xsec = 0.000321; xsec_unc = 0.22118380062305298; return;}
    else if (glu_mass == 2225) { xsec = 0.000313; xsec_unc = 0.22172523961661342; return;}
    else if (glu_mass == 2230) { xsec = 0.000305; xsec_unc = 0.22262295081967212; return;}
    else if (glu_mass == 2235) { xsec = 0.000297; xsec_unc = 0.22323232323232323; return;}
    else if (glu_mass == 2240) { xsec = 0.00029; xsec_unc = 0.22413793103448273; return;}
    else if (glu_mass == 2245) { xsec = 0.000283; xsec_unc = 0.2247349823321555; return;}
    else if (glu_mass == 2250) { xsec = 0.000275; xsec_unc = 0.22545454545454546; return;}
    else if (glu_mass == 2255) { xsec = 0.000268; xsec_unc = 0.22611940298507463; return;}
    else if (glu_mass == 2260) { xsec = 0.000262; xsec_unc = 0.22709923664122136; return;}
    else if (glu_mass == 2265) { xsec = 0.000255; xsec_unc = 0.22784313725490196; return;}
    else if (glu_mass == 2270) { xsec = 0.000248; xsec_unc = 0.22862903225806452; return;}
    else if (glu_mass == 2275) { xsec = 0.000242; xsec_unc = 0.22933884297520662; return;}
    else if (glu_mass == 2280) { xsec = 0.000236; xsec_unc = 0.23008474576271187; return;}
    else if (glu_mass == 2285) { xsec = 0.00023; xsec_unc = 0.2308695652173913; return;}
    else if (glu_mass == 2290) { xsec = 0.000224; xsec_unc = 0.23169642857142858; return;}
    else if (glu_mass == 2295) { xsec = 0.000219; xsec_unc = 0.23242009132420088; return;}
    else if (glu_mass == 2300) { xsec = 0.000213; xsec_unc = 0.23333333333333334; return;}
    else if (glu_mass == 2305) { xsec = 0.000208; xsec_unc = 0.23413461538461539; return;}
    else if (glu_mass == 2310) { xsec = 0.000202; xsec_unc = 0.23514851485148516; return;}
    else if (glu_mass == 2315) { xsec = 0.000197; xsec_unc = 0.23604060913705585; return;}
    else if (glu_mass == 2320) { xsec = 0.000192; xsec_unc = 0.23697916666666666; return;}
    else if (glu_mass == 2325) { xsec = 0.000187; xsec_unc = 0.23743315508021393; return;}
    else if (glu_mass == 2330) { xsec = 0.000183; xsec_unc = 0.23879781420765026; return;}
    else if (glu_mass == 2335) { xsec = 0.000178; xsec_unc = 0.2393258426966292; return;}
    else if (glu_mass == 2340) { xsec = 0.000174; xsec_unc = 0.24022988505747125; return;}
    else if (glu_mass == 2345) { xsec = 0.000169; xsec_unc = 0.2414201183431953; return;}
    else if (glu_mass == 2350) { xsec = 0.000165; xsec_unc = 0.24242424242424246; return;}
    else if (glu_mass == 2355) { xsec = 0.000161; xsec_unc = 0.24285714285714285; return;}
    else if (glu_mass == 2360) { xsec = 0.000157; xsec_unc = 0.2439490445859873; return;}
    else if (glu_mass == 2365) { xsec = 0.000153; xsec_unc = 0.24509803921568624; return;}
    else if (glu_mass == 2370) { xsec = 0.000149; xsec_unc = 0.24630872483221478; return;}
    else if (glu_mass == 2375) { xsec = 0.000145; xsec_unc = 0.24689655172413794; return;}
    else if (glu_mass == 2380) { xsec = 0.000142; xsec_unc = 0.24788732394366197; return;}
    else if (glu_mass == 2385) { xsec = 0.000138; xsec_unc = 0.2485507246376812; return;}
    else if (glu_mass == 2390) { xsec = 0.000134; xsec_unc = 0.25; return;}
    else if (glu_mass == 2395) { xsec = 0.000131; xsec_unc = 0.2511450381679389; return;}
    else if (glu_mass == 2400) { xsec = 0.000128; xsec_unc = 0.25156249999999997; return;}
    else if (glu_mass == 2405) { xsec = 0.000125; xsec_unc = 0.2528; return;}
    else if (glu_mass == 2410) { xsec = 0.000121; xsec_unc = 0.2537190082644628; return;}
    else if (glu_mass == 2415) { xsec = 0.000118; xsec_unc = 0.2550847457627119; return;}
    else if (glu_mass == 2420) { xsec = 0.000115; xsec_unc = 0.25652173913043474; return;}
    else if (glu_mass == 2425) { xsec = 0.000113; xsec_unc = 0.2575221238938053; return;}
    else if (glu_mass == 2430) { xsec = 0.00011; xsec_unc = 0.2581818181818182; return;}
    else if (glu_mass == 2435) { xsec = 0.000107; xsec_unc = 0.2588785046728972; return;}
    else if (glu_mass == 2440) { xsec = 0.000104; xsec_unc = 0.2605769230769231; return;}
    else if (glu_mass == 2445) { xsec = 0.000102; xsec_unc = 0.26176470588235295; return;}
    else if (glu_mass == 2450) { xsec = 9.91e-05; xsec_unc = 0.2623612512613522; return;}
    else if (glu_mass == 2455) { xsec = 9.66e-05; xsec_unc = 0.2639751552795031; return;}
    else if (glu_mass == 2460) { xsec = 9.41e-05; xsec_unc = 0.26461211477151964; return;}
    else if (glu_mass == 2465) { xsec = 9.18e-05; xsec_unc = 0.2657952069716776; return;}
    else if (glu_mass == 2470) { xsec = 8.95e-05; xsec_unc = 0.26703910614525145; return;}
    else if (glu_mass == 2475) { xsec = 8.72e-05; xsec_unc = 0.268348623853211; return;}
    else if (glu_mass == 2480) { xsec = 8.5e-05; xsec_unc = 0.26941176470588235; return;}
    else if (glu_mass == 2485) { xsec = 8.29e-05; xsec_unc = 0.27020506634499397; return;}
    else if (glu_mass == 2490) { xsec = 8.08e-05; xsec_unc = 0.2722772277227723; return;}
    else if (glu_mass == 2495) { xsec = 7.88e-05; xsec_unc = 0.27284263959390864; return;}
    else if (glu_mass == 2500) { xsec = 7.68e-05; xsec_unc = 0.27473958333333337; return;}
    else if (glu_mass == 2505) { xsec = 7.49e-05; xsec_unc = 0.27636849132176233; return;}
    else if (glu_mass == 2510) { xsec = 7.3e-05; xsec_unc = 0.27671232876712326; return;}
    else if (glu_mass == 2515) { xsec = 7.12e-05; xsec_unc = 0.27808988764044945; return;}
    else if (glu_mass == 2520) { xsec = 6.94e-05; xsec_unc = 0.2795389048991354; return;}
    else if (glu_mass == 2525) { xsec = 6.77e-05; xsec_unc = 0.28064992614475626; return;}
    else if (glu_mass == 2530) { xsec = 6.6e-05; xsec_unc = 0.2818181818181818; return;}
    else if (glu_mass == 2535) { xsec = 6.43e-05; xsec_unc = 0.28304821150855364; return;}
    else if (glu_mass == 2540) { xsec = 6.27e-05; xsec_unc = 0.28548644338118023; return;}
    else if (glu_mass == 2545) { xsec = 6.11e-05; xsec_unc = 0.2864157119476268; return;}
    else if (glu_mass == 2550) { xsec = 5.96e-05; xsec_unc = 0.28859060402684567; return;}
    else if (glu_mass == 2555) { xsec = 5.81e-05; xsec_unc = 0.28915662650602403; return;}
    else if (glu_mass == 2560) { xsec = 5.66e-05; xsec_unc = 0.2915194346289753; return;}
    else if (glu_mass == 2565) { xsec = 5.52e-05; xsec_unc = 0.29166666666666663; return;}
    else if (glu_mass == 2570) { xsec = 5.38e-05; xsec_unc = 0.2936802973977695; return;}
    else if (glu_mass == 2575) { xsec = 5.25e-05; xsec_unc = 0.29523809523809524; return;}
    else if (glu_mass == 2580) { xsec = 5.12e-05; xsec_unc = 0.296875; return;}
    else if (glu_mass == 2585) { xsec = 4.99e-05; xsec_unc = 0.2985971943887776; return;}
    else if (glu_mass == 2590) { xsec = 4.86e-05; xsec_unc = 0.3004115226337449; return;}
    else if (glu_mass == 2595) { xsec = 4.74e-05; xsec_unc = 0.30168776371308015; return;}
    else if (glu_mass == 2600) { xsec = 4.62e-05; xsec_unc = 0.30303030303030304; return;}
    else if (glu_mass == 2605) { xsec = 4.51e-05; xsec_unc = 0.30376940133037694; return;}
    else if (glu_mass == 2610) { xsec = 4.39e-05; xsec_unc = 0.3052391799544419; return;}
    else if (glu_mass == 2615) { xsec = 4.28e-05; xsec_unc = 0.30841121495327106; return;}
    else if (glu_mass == 2620) { xsec = 4.18e-05; xsec_unc = 0.30861244019138756; return;}
    else if (glu_mass == 2625) { xsec = 4.07e-05; xsec_unc = 0.31203931203931207; return;}
    else if (glu_mass == 2630) { xsec = 3.97e-05; xsec_unc = 0.31234256926952136; return;}
    else if (glu_mass == 2635) { xsec = 3.87e-05; xsec_unc = 0.3152454780361757; return;}
    else if (glu_mass == 2640) { xsec = 3.77e-05; xsec_unc = 0.3156498673740053; return;}
    else if (glu_mass == 2645) { xsec = 3.68e-05; xsec_unc = 0.3179347826086956; return;}
    else if (glu_mass == 2650) { xsec = 3.59e-05; xsec_unc = 0.3203342618384401; return;}
    else if (glu_mass == 2655) { xsec = 3.5e-05; xsec_unc = 0.3228571428571429; return;}
    else if (glu_mass == 2660) { xsec = 3.41e-05; xsec_unc = 0.3225806451612903; return;}
    else if (glu_mass == 2665) { xsec = 3.32e-05; xsec_unc = 0.3253012048192771; return;}
    else if (glu_mass == 2670) { xsec = 3.24e-05; xsec_unc = 0.3271604938271605; return;}
    else if (glu_mass == 2675) { xsec = 3.16e-05; xsec_unc = 0.3291139240506329; return;}
    else if (glu_mass == 2680) { xsec = 3.08e-05; xsec_unc = 0.33116883116883117; return;}
    else if (glu_mass == 2685) { xsec = 3e-05; xsec_unc = 0.33333333333333337; return;}
    else if (glu_mass == 2690) { xsec = 2.93e-05; xsec_unc = 0.3351535836177474; return;}
    else if (glu_mass == 2695) { xsec = 2.85e-05; xsec_unc = 0.3371929824561403; return;}
    else if (glu_mass == 2700) { xsec = 2.78e-05; xsec_unc = 0.33920863309352517; return;}
    else if (glu_mass == 2705) { xsec = 2.71e-05; xsec_unc = 0.3413284132841328; return;}
    else if (glu_mass == 2710) { xsec = 2.65e-05; xsec_unc = 0.3433962264150943; return;}
    else if (glu_mass == 2715) { xsec = 2.58e-05; xsec_unc = 0.3457364341085271; return;}
    else if (glu_mass == 2720) { xsec = 2.51e-05; xsec_unc = 0.347808764940239; return;}
    else if (glu_mass == 2725) { xsec = 2.45e-05; xsec_unc = 0.3497959183673469; return;}
    else if (glu_mass == 2730) { xsec = 2.39e-05; xsec_unc = 0.3523012552301255; return;}
    else if (glu_mass == 2735) { xsec = 2.33e-05; xsec_unc = 0.35450643776824037; return;}
    else if (glu_mass == 2740) { xsec = 2.27e-05; xsec_unc = 0.3568281938325991; return;}
    else if (glu_mass == 2745) { xsec = 2.21e-05; xsec_unc = 0.35882352941176476; return;}
    else if (glu_mass == 2750) { xsec = 2.16e-05; xsec_unc = 0.3611111111111111; return;}
    else if (glu_mass == 2755) { xsec = 2.11e-05; xsec_unc = 0.36350710900473926; return;}
    else if (glu_mass == 2760) { xsec = 2.05e-05; xsec_unc = 0.36585365853658536; return;}
    else if (glu_mass == 2765) { xsec = 2e-05; xsec_unc = 0.36799999999999994; return;}
    else if (glu_mass == 2770) { xsec = 1.95e-05; xsec_unc = 0.37025641025641026; return;}
    else if (glu_mass == 2775) { xsec = 1.9e-05; xsec_unc = 0.3731578947368421; return;}
    else if (glu_mass == 2780) { xsec = 1.85e-05; xsec_unc = 0.37513513513513513; return;}
    else if (glu_mass == 2785) { xsec = 1.81e-05; xsec_unc = 0.37734806629834255; return;}
    else if (glu_mass == 2790) { xsec = 1.76e-05; xsec_unc = 0.3801136363636364; return;}
    else if (glu_mass == 2795) { xsec = 1.72e-05; xsec_unc = 0.38255813953488366; return;}
    else if (glu_mass == 2800) { xsec = 1.68e-05; xsec_unc = 0.38452380952380955; return;}
    else if (glu_mass == 2805) { xsec = 1.63e-05; xsec_unc = 0.3871165644171779; return;}
    else if (glu_mass == 2810) { xsec = 1.59e-05; xsec_unc = 0.3893081761006289; return;}
    else if (glu_mass == 2815) { xsec = 1.55e-05; xsec_unc = 0.39161290322580644; return;}
    else if (glu_mass == 2820) { xsec = 1.51e-05; xsec_unc = 0.39403973509933776; return;}
    else if (glu_mass == 2825) { xsec = 1.48e-05; xsec_unc = 0.39662162162162157; return;}
    else if (glu_mass == 2830) { xsec = 1.44e-05; xsec_unc = 0.3993055555555556; return;}
    else if (glu_mass == 2835) { xsec = 1.4e-05; xsec_unc = 0.40142857142857147; return;}
    else if (glu_mass == 2840) { xsec = 1.37e-05; xsec_unc = 0.4036496350364964; return;}
    else if (glu_mass == 2845) { xsec = 1.33e-05; xsec_unc = 0.406015037593985; return;}
    else if (glu_mass == 2850) { xsec = 1.3e-05; xsec_unc = 0.4084615384615385; return;}
    else if (glu_mass == 2855) { xsec = 1.27e-05; xsec_unc = 0.4110236220472441; return;}
    else if (glu_mass == 2860) { xsec = 1.24e-05; xsec_unc = 0.4137096774193548; return;}
    else if (glu_mass == 2865) { xsec = 1.21e-05; xsec_unc = 0.415702479338843; return;}
    else if (glu_mass == 2870) { xsec = 1.18e-05; xsec_unc = 0.4186440677966102; return;}
    else if (glu_mass == 2875) { xsec = 1.15e-05; xsec_unc = 0.42086956521739133; return;}
    else if (glu_mass == 2880) { xsec = 1.12e-05; xsec_unc = 0.42321428571428577; return;}
    else if (glu_mass == 2885) { xsec = 1.09e-05; xsec_unc = 0.4256880733944953; return;}
    else if (glu_mass == 2890) { xsec = 1.06e-05; xsec_unc = 0.4283018867924528; return;}
    else if (glu_mass == 2895) { xsec = 1.04e-05; xsec_unc = 0.4307692307692308; return;}
    else if (glu_mass == 2900) { xsec = 1.01e-05; xsec_unc = 0.43267326732673267; return;}
    else if (glu_mass == 2905) { xsec = 9.86e-06; xsec_unc = 0.4350912778904665; return;}
    else if (glu_mass == 2910) { xsec = 9.61e-06; xsec_unc = 0.43808532778355885; return;}
    else if (glu_mass == 2915) { xsec = 9.37e-06; xsec_unc = 0.4407684098185699; return;}
    else if (glu_mass == 2920) { xsec = 9.14e-06; xsec_unc = 0.4431072210065645; return;}
    else if (glu_mass == 2925) { xsec = 8.91e-06; xsec_unc = 0.44556677890011226; return;}
    else if (glu_mass == 2930) { xsec = 8.69e-06; xsec_unc = 0.44764096662830843; return;}
    else if (glu_mass == 2935) { xsec = 8.48e-06; xsec_unc = 0.45047169811320753; return;}
    else if (glu_mass == 2940) { xsec = 8.27e-06; xsec_unc = 0.45344619105199513; return;}
    else if (glu_mass == 2945) { xsec = 8.06e-06; xsec_unc = 0.4553349875930521; return;}
    else if (glu_mass == 2950) { xsec = 7.86e-06; xsec_unc = 0.45801526717557256; return;}
    else if (glu_mass == 2955) { xsec = 7.67e-06; xsec_unc = 0.4602346805736637; return;}
    else if (glu_mass == 2960) { xsec = 7.48e-06; xsec_unc = 0.46256684491978606; return;}
    else if (glu_mass == 2965) { xsec = 7.29e-06; xsec_unc = 0.4650205761316873; return;}
    else if (glu_mass == 2970) { xsec = 7.11e-06; xsec_unc = 0.46835443037974683; return;}
    else if (glu_mass == 2975) { xsec = 6.94e-06; xsec_unc = 0.47118155619596547; return;}
    else if (glu_mass == 2980) { xsec = 6.77e-06; xsec_unc = 0.4726735598227474; return;}
    else if (glu_mass == 2985) { xsec = 6.6e-06; xsec_unc = 0.4757575757575757; return;}
    else if (glu_mass == 2990) { xsec = 6.44e-06; xsec_unc = 0.4782608695652174; return;}
    else if (glu_mass == 2995) { xsec = 6.28e-06; xsec_unc = 0.48089171974522293; return;}
    else if (glu_mass == 3000) { xsec = 6.12e-06; xsec_unc = 0.48366013071895425; return;}
    else { xsec = 0; xsec_unc = 0; }
  }

  // xsec is pb, xsec_unc is relative uncertainty
  void higgsino2DCrossSection(int hig_mass, double &xsec, double &xsec_unc) {
    if(hig_mass == 127) { xsec = .5824*.5824*1.44725; xsec_unc = 0.0395277; return;}
    else if(hig_mass == 150) { xsec = .5824*.5824*0.71514; xsec_unc = 0.0421496; return;}
    else if(hig_mass == 175) { xsec = .5824*.5824*0.419059; xsec_unc = 0.0453279; return;}
    else if(hig_mass == 200) { xsec = .5824*.5824*0.244213; xsec_unc = 0.047925; return;}
    else if(hig_mass == 225) { xsec = .5824*.5824*0.156286; xsec_unc = 0.0502876; return;}
    else if(hig_mass == 250) { xsec = .5824*.5824*0.104252; xsec_unc = 0.0526169; return;}
    else if(hig_mass == 275) { xsec = .5824*.5824*0.0719125; xsec_unc = 0.0549666; return;}
    else if(hig_mass == 300) { xsec = .5824*.5824*0.0509994; xsec_unc = 0.0572762; return;}
    else if(hig_mass == 325) { xsec = .5824*.5824*0.0369715; xsec_unc = 0.0590317; return;}
    else if(hig_mass == 350) { xsec = .5824*.5824*0.0273286; xsec_unc = 0.0607766; return;}
    else if(hig_mass == 375) { xsec = .5824*.5824*0.0205429; xsec_unc = 0.0625031; return;}
    else if(hig_mass == 400) { xsec = .5824*.5824*0.0156691; xsec_unc = 0.0642085; return;}
    else if(hig_mass == 425) { xsec = .5824*.5824*0.0120965; xsec_unc = 0.0657801; return;}
    else if(hig_mass == 450) { xsec = .5824*.5824*0.00944017; xsec_unc = 0.0674544; return;}
    else if(hig_mass == 475) { xsec = .5824*.5824*0.00743587; xsec_unc = 0.0686033; return;}
    else if(hig_mass == 500) { xsec = .5824*.5824*0.00590757; xsec_unc = 0.0699909; return;}
    else if(hig_mass == 525) { xsec = .5824*.5824*0.00473235; xsec_unc = 0.0713166; return;}
    else if(hig_mass == 550) { xsec = .5824*.5824*0.0038167; xsec_unc = 0.0722834; return;}
    else if(hig_mass == 575) { xsec = .5824*.5824*0.00309847; xsec_unc = 0.0739435; return;}
    else if(hig_mass == 600) { xsec = .5824*.5824*0.00253015; xsec_unc = 0.0754291; return;}
    else if(hig_mass == 625) { xsec = .5824*.5824*0.00207755; xsec_unc = 0.0763142; return;}
    else if(hig_mass == 650) { xsec = .5824*.5824*0.00171418; xsec_unc = 0.0775695; return;}
    else if(hig_mass == 675) { xsec = .5824*.5824*0.0014199; xsec_unc = 0.0782907; return;}
    else if(hig_mass == 700) { xsec = .5824*.5824*0.00118113; xsec_unc = 0.0796388; return;}
    else if(hig_mass == 725) { xsec = .5824*.5824*0.00098639; xsec_unc = 0.0809291; return;}
    else if(hig_mass == 750) { xsec = .5824*.5824*0.000826366; xsec_unc = 0.081879; return;}
    else if(hig_mass == 775) { xsec = .5824*.5824*0.000694985; xsec_unc = 0.0841352; return;}
    else if(hig_mass == 800) { xsec = .5824*.5824*0.000586211; xsec_unc = 0.0862527; return;}
    else if(hig_mass == 825) { xsec = .5824*.5824*0.000495914; xsec_unc = 0.0863944; return;}
    else if(hig_mass == 850) { xsec = .5824*.5824*0.000420556; xsec_unc = 0.085742; return;}
    else if(hig_mass == 875) { xsec = .5824*.5824*0.000361029; xsec_unc = 0.0888678; return;}
    else if(hig_mass == 900) { xsec = .5824*.5824*0.000305935; xsec_unc = 0.0912439; return;}
    else if(hig_mass == 925) { xsec = .5824*.5824*0.000262621; xsec_unc = 0.0913227; return;}
    else if(hig_mass == 950) { xsec = .5824*.5824*0.00022285; xsec_unc = 0.0919538; return;}
    else if(hig_mass == 975) { xsec = .5824*.5824*0.0001909; xsec_unc = 0.0937616; return;}
    else if(hig_mass == 1000) { xsec = .5824*.5824*0.00016428; xsec_unc = 0.0954285; return;}
    else if(hig_mass == 1025) { xsec = .5824*.5824*0.00014139; xsec_unc = 0.095765; return;}
    else if(hig_mass == 1050) { xsec = .5824*.5824*0.000121865; xsec_unc = 0.0967595; return;}
    else if(hig_mass == 1075) { xsec = .5824*.5824*0.000105913; xsec_unc = 0.0978621; return;}
    else if(hig_mass == 1100) { xsec = .5824*.5824*9.12469e-05; xsec_unc = 0.0964142; return;}
    else if(hig_mass == 1125) { xsec = .5824*.5824*7.93058e-05; xsec_unc = 0.0999433; return;}
    else if(hig_mass == 1150) { xsec = .5824*.5824*6.84561e-05; xsec_unc = 0.103594; return;}
    else if(hig_mass == 1175) { xsec = .5824*.5824*5.93602e-05; xsec_unc = 0.10201; return;}
    else if(hig_mass == 1200) { xsec = .5824*.5824*5.16263e-05; xsec_unc = 0.102499; return;}
    else if(hig_mass == 1225) { xsec = .5824*.5824*4.4906e-05; xsec_unc = 0.104071; return;}
    else if(hig_mass == 1250) { xsec = .5824*.5824*3.91587e-05; xsec_unc = 0.104736; return;}
    else if(hig_mass == 1275) { xsec = .5824*.5824*3.43135e-05; xsec_unc = 0.10615; return;}
    else if(hig_mass == 1300) { xsec = .5824*.5824*2.99353e-05; xsec_unc = 0.10783; return;}
    else if(hig_mass == 1325) { xsec = .5824*.5824*2.62223e-05; xsec_unc = 0.108061; return;}
    else if(hig_mass == 1350) { xsec = .5824*.5824*2.28072e-05; xsec_unc = 0.109427; return;}
    else if(hig_mass == 1375) { xsec = .5824*.5824*2.00393e-05; xsec_unc = 0.109789; return;}
    else if(hig_mass == 1400) { xsec = .5824*.5824*1.75031e-05; xsec_unc = 0.111631; return;}
    else if(hig_mass == 1425) { xsec = .5824*.5824*1.53144e-05; xsec_unc = 0.11145; return;}
    else if(hig_mass == 1450) { xsec = .5824*.5824*1.34572e-05; xsec_unc = 0.11084; return;}
    else if(hig_mass == 1475) { xsec = .5824*.5824*1.17047e-05; xsec_unc = 0.113027; return;}
    else{ xsec = 0; xsec_unc = 0;}
  }

  // xsec is pb, xsec_unc is relative uncertainty
  void higgsinoCrossSection(int hig_mass, double &xsec, double &xsec_unc) {
    if(hig_mass ==127) { xsec = .5824*.5824*7.6022; xsec_unc = 0.0393921; return;}
    else if(hig_mass ==150) { xsec = .5824*.5824*3.83231; xsec_unc = 0.0413612; return;}
    else if(hig_mass ==175) { xsec = .5824*.5824*2.26794; xsec_unc = 0.044299; return;}
    else if(hig_mass ==200) { xsec = .5824*.5824*1.33562; xsec_unc = 0.0474362; return;}
    else if(hig_mass ==225) { xsec = .5824*.5824*0.860597; xsec_unc = 0.0504217; return;}
    else if(hig_mass ==250) { xsec = .5824*.5824*0.577314; xsec_unc = 0.0532731; return;}
    else if(hig_mass ==275) { xsec = .5824*.5824*0.400107; xsec_unc = 0.0560232; return;}
    else if(hig_mass ==300) { xsec = .5824*.5824*0.284855; xsec_unc = 0.0586867; return;}
    else if(hig_mass ==325) { xsec = .5824*.5824*0.20736; xsec_unc = 0.0613554; return;}
    else if(hig_mass ==350) { xsec = .5824*.5824*0.153841; xsec_unc = 0.0640598; return;}
    else if(hig_mass ==375) { xsec = .5824*.5824*0.116006; xsec_unc = 0.066892; return;}
    else if(hig_mass ==400) { xsec = .5824*.5824*0.0887325; xsec_unc = 0.0697517; return;}
    else if(hig_mass ==425) { xsec = .5824*.5824*0.0686963; xsec_unc = 0.0723531; return;}
    else if(hig_mass ==450) { xsec = .5824*.5824*0.0537702; xsec_unc = 0.0748325; return;}
    else if(hig_mass ==475) { xsec = .5824*.5824*0.0424699; xsec_unc = 0.0775146; return;}
    else if(hig_mass ==500) { xsec = .5824*.5824*0.0338387; xsec_unc = 0.0802572; return;}
    else if(hig_mass ==525) { xsec = .5824*.5824*0.0271867; xsec_unc = 0.0825803; return;}
    else if(hig_mass ==550) { xsec = .5824*.5824*0.0219868; xsec_unc = 0.0849278; return;}
    else if(hig_mass ==575) { xsec = .5824*.5824*0.0179062; xsec_unc = 0.087561; return;}
    else if(hig_mass ==600) { xsec = .5824*.5824*0.0146677; xsec_unc = 0.0900693; return;}
    else if(hig_mass ==625) { xsec = .5824*.5824*0.012062; xsec_unc = 0.091959; return;}
    else if(hig_mass ==650) { xsec = .5824*.5824*0.00996406; xsec_unc = 0.094065; return;}
    else if(hig_mass ==675) { xsec = .5824*.5824*0.00828246; xsec_unc = 0.0957436; return;}
    else if(hig_mass ==700) { xsec = .5824*.5824*0.00689981; xsec_unc = 0.0982894; return;}
    else if(hig_mass ==725) { xsec = .5824*.5824*0.00578355; xsec_unc = 0.0999915; return;}
    else if(hig_mass ==750) { xsec = .5824*.5824*0.0048731; xsec_unc = 0.101211; return;}
    else if(hig_mass ==775) { xsec = .5824*.5824*0.00409781; xsec_unc = 0.104646; return;}
    else if(hig_mass ==800) { xsec = .5824*.5824*0.00346143; xsec_unc = 0.107618; return;}
    else if(hig_mass ==825) { xsec = .5824*.5824*0.0029337; xsec_unc = 0.108353; return;}
    else if(hig_mass ==850) { xsec = .5824*.5824*0.0024923; xsec_unc = 0.110016; return;}
    else if(hig_mass ==875) { xsec = .5824*.5824*0.00213679; xsec_unc = 0.112636; return;}
    else if(hig_mass ==900) { xsec = .5824*.5824*0.00180616; xsec_unc = 0.1134; return;}
    else if(hig_mass ==925) { xsec = .5824*.5824*0.00155453; xsec_unc = 0.116949; return;}
    else if(hig_mass ==950) { xsec = .5824*.5824*0.00132692; xsec_unc = 0.117027; return;}
    else if(hig_mass ==975) { xsec = .5824*.5824*0.00112975; xsec_unc = 0.121244; return;}
    else if(hig_mass ==1000) { xsec = .5824*.5824*0.000968853; xsec_unc = 0.126209; return;}
    else if(hig_mass ==1025) { xsec = .5824*.5824*0.000840602; xsec_unc = 0.121654; return;}
    else if(hig_mass ==1050) { xsec = .5824*.5824*0.000731306; xsec_unc = 0.118502; return;}
    else if(hig_mass ==1075) { xsec = .5824*.5824*0.000627083; xsec_unc = 0.127723; return;}
    else if(hig_mass ==1100) { xsec = .5824*.5824*0.000538005; xsec_unc = 0.134099; return;}
    else if(hig_mass ==1125) { xsec = .5824*.5824*0.00046747; xsec_unc = 0.133755; return;}
    else if(hig_mass ==1150) { xsec = .5824*.5824*0.000405108; xsec_unc = 0.120607; return;}
    else if(hig_mass ==1175) { xsec = .5824*.5824*0.000348261; xsec_unc = 0.139744; return;}
    else if(hig_mass ==1200) { xsec = .5824*.5824*0.000299347; xsec_unc = 0.162604; return;}
    else if(hig_mass ==1225) { xsec = .5824*.5824*0.000265935; xsec_unc = 0.137575; return;}
    else if(hig_mass ==1250) { xsec = .5824*.5824*0.000240471; xsec_unc = 0.119271; return;}
    else if(hig_mass ==1275) { xsec = .5824*.5824*0.000190411; xsec_unc = 0.138061; return;}
    else if(hig_mass ==1300) { xsec = .5824*.5824*0.000160765; xsec_unc = 0.122224; return;}
    else if(hig_mass ==1325) { xsec = .5824*.5824*0.000136272; xsec_unc = 0.138533; return;}
    else if(hig_mass ==1350) { xsec = .5824*.5824*0.000111174; xsec_unc = 0.177681; return;}
    else if(hig_mass ==1375) { xsec = .5824*.5824*9.74728e-05; xsec_unc = 0.138992; return;}
    else if(hig_mass ==1400) { xsec = .5824*.5824*7.80263e-05; xsec_unc = 0.118718; return;}
    else if(hig_mass ==1425) { xsec = .5824*.5824*6.96843e-05; xsec_unc = 0.139439; return;}
    else if(hig_mass ==1450) { xsec = .5824*.5824*6.96962e-05; xsec_unc = 0.198887; return;}
    else if(hig_mass ==1475) { xsec = .5824*.5824*4.98006e-05; xsec_unc = 0.139874; return;}
    else{ xsec = 0; xsec_unc = 0;}
  }
}
