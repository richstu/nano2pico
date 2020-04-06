// Lists the MC cross sections

#include <iostream>
#include "cross_sections.hpp"

using namespace std;

namespace xsec{

  float crossSection(const TString &file, bool is2016){
    float xsec(-999999.), Htobb(0.5824), HToZG(0.001533), ZToLL(0.100974);

    if (is2016) {
        //  Cross-section taken from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
        // Alternative option: https://twiki.cern.ch/twiki/bin/view/Sandbox/FullNNLOcrossSections#Top_cross_section_for_13_TeV
        if(file.Contains("TTJets_Tune") || file.Contains("TT_"))  xsec = 815.96;
        if(file.Contains("TTJets_HT")){//LO cross sections with k-factor of 1.625 already applied
          if(file.Contains("2500toInf")) xsec = 0.0023234211;
          if(file.Contains("1200to2500")) xsec = 0.194972521;
          if(file.Contains("800to1200")) xsec = 1.07722318;
          if(file.Contains("600to800")) xsec = 2.61537118;

        }
        // The efficiency of the mtt>1000 cut is taken from sigma(mtt>1000)/sigma(inclusive) from mcm
        double mtt_1000_eff=(11.16/670.3);
        if(file.Contains("TTJets_Mtt-1000toInf")) xsec = 815.96*mtt_1000_eff;

        if(file.Contains("TTJets_DiLept") || file.Contains("TTTo2L2Nu")) xsec = 85.66; // (3*0.108)^2*815.96
        if(file.Contains("TTJets_SingleLept")) xsec = 178.7; //(1- ((1-3*0.108)^2+(3*0.108)^2))*815.96*0.5 per half
        if(file.Contains("TTToSemiLeptonic")) xsec = 357.4;
        
        if(file.Contains("TTJets_DiLept_genMET-150")) xsec = 0.06392*85.66; // filter_eff*(3*0.108)^2*815.96
        if(file.Contains("TTJets_SingleLept") && file.Contains("genMET-150")) xsec = 0.05217*178.7; //filter_eff*(1- ((1-3*0.108)^2+(3*0.108)^2))*815.96*0.5 per half

        // cross sections from mcm
        if(file.Contains("TTG")) xsec = 3.697;                
        if(file.Contains("TTTT_Tune")) xsec = 0.009103;
        // mcm cross section times the same kfactors used for leptonic samples
        if(file.Contains("WJetsToQQ_HT-600ToInf")) xsec = 95.14*1.21;
        if(file.Contains("ZJetsToQQ_HT600toInf")) xsec = 5.67*1.23;
        
        // From https://cms-pdmv.cern.ch/mcm
        // k-factors from https://mangano.web.cern.ch/mangano/public/MECCA/samples_50ns_miniaod.txt
        // k-factors are ratio of https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
        // NLO/NNLO cross-sections to that of an inclusive sample in mcm at lower order (LO/NLO)

        if(file.Contains("WJetsToLNu_Tune")) xsec=61526.7; //NNLO from Lesya's summary table 

        //cross-section per slice changed due to change in genHT definition
        if(file.Contains("WJetsToLNu_HT-70To100"))    xsec = 1372.*1.21; 
        if(file.Contains("WJetsToLNu_HT-100To200"))   xsec = 1347.*1.21; 
        if(file.Contains("WJetsToLNu_HT-200To400"))   xsec = 360.*1.21;
        if(file.Contains("WJetsToLNu_HT-400To600"))   xsec = 48.98*1.21;
        if(file.Contains("WJetsToLNu_HT-600ToInf"))   xsec = 18.77*1.21;
        if(file.Contains("WJetsToLNu_HT-600To800"))   xsec = 12.05*1.21;
        if(file.Contains("WJetsToLNu_HT-800To1200"))  xsec = 5.501*1.21;
        if(file.Contains("WJetsToLNu_HT-1200To2500")) xsec = 1.329*1.21;
        if(file.Contains("WJetsToLNu_HT-2500ToInf"))  xsec = 0.03216*1.21;

        //updated 02-2019 with XSDB
        if(file.Contains("QCD_HT100to200_Tune"))   xsec = 28060000;
        if(file.Contains("QCD_HT200to300_Tune"))   xsec = 1710000;
        if(file.Contains("QCD_HT300to500_Tune"))   xsec = 347500;
        if(file.Contains("QCD_HT500to700_Tune"))   xsec = 32060;
        if(file.Contains("QCD_HT700to1000_Tune"))  xsec = 6829;
        if(file.Contains("QCD_HT1000to1500_Tune")) xsec = 1207;
        if(file.Contains("QCD_HT1500to2000_Tune")) xsec = 120.0;
        if(file.Contains("QCD_HT2000toInf_Tune"))  xsec = 25.25;

        // Cross sections from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec
        // multiplied by BF(W->mu,e,tau) = 0.324
        if (file.Contains("ST_s-channel_4f_leptonDecays"))     xsec = 3.34;
        if (file.Contains("ST_t-channel_antitop_4f_InclusiveDecays")) xsec = 80.95;
        if (file.Contains("ST_t-channel_top_4f_InclusiveDecays")) xsec = 136.02;
        if (file.Contains("ST_tW_antitop_5f_NoFullyHadronicDecays"))     xsec = 35.85*0.543;
        if (file.Contains("ST_tW_top_5f_NoFullyHadronicDecays"))     xsec = 35.85*0.543; 

        if(file.Contains("DYJetsToLL_M-10to50_Tune")) xsec = 18610*1.23;
        if(file.Contains("DYJetsToLL_M-50_Tune"))     xsec = 4895*1.23;

        if(file.Contains("DYJetsToLL_M-50_HT-70to100"))    xsec = 175.3*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-100to200"))   xsec = 139.4*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-200to400"))   xsec = 42.75*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-400to600"))   xsec = 5.497*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-600to800"))   xsec = 1.363*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-800to1200"))  xsec = 0.6759*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-1200to2500")) xsec = 0.116*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-2500toInf"))  xsec = 0.002592*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-600toInf"))   xsec = 2.21*1.23;

        if(file.Contains("ZJetsToNuNu_HT-100To200"))       xsec = 280.35*1.27;
        if(file.Contains("ZJetsToNuNu_HT-200To400"))       xsec = 77.67*1.27;
        if(file.Contains("ZJetsToNuNu_HT-400To600"))       xsec = 10.73*1.27;
        if(file.Contains("ZJetsToNuNu_HT-600To800"))       xsec = 2.536*1.27;
        if(file.Contains("ZJetsToNuNu_HT-800To1200"))      xsec = 1.161*1.27;
        if(file.Contains("ZJetsToNuNu_HT-1200To2500"))     xsec = 0.2824*1.27;
        if(file.Contains("ZJetsToNuNu_HT-2500ToInf"))      xsec = 0.006459*1.27;
        if(file.Contains("ZJetsToNuNu_HT-600ToInf"))       xsec = 3.986*1.27;

        if(file.Contains("TTZToQQ"))                       xsec = 0.5297;
        if(file.Contains("TTZToLLNuNu_M-10"))              xsec = 0.2529;
        if(file.Contains("TTWJetsToQQ"))                   xsec = 0.4062;
        if(file.Contains("TTWJetsToLNu"))                  xsec = 0.2043;
       
        //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
        if(file.Contains("WWTo2L2Nu"))    xsec = 12.178; //NNLO
        if(file.Contains("WWToLNuQQ"))    xsec = 49.997; //NNLO
        if(file.Contains("ttHTobb_M125")) xsec = 0.2934;

        if(file.Contains("WZTo1L3Nu"))    xsec = 3.05;
        if(file.Contains("WZTo1L1Nu2Q"))  xsec = 10.96;
        if(file.Contains("WZTo2L2Q"))     xsec = 5.595;
        if(file.Contains("WZTo3LNu"))     xsec = 4.42965;
        if(file.Contains("VVTo2L2Nu"))    xsec = 11.95;
        if(file.Contains("ZZ_Tune"))      xsec = 16.523;

	//cross sections stolen from AN for TOP-18-009
	
        if(file.Contains("tZq"))   xsec = 0.09418;
	if(file.Contains("ttHToNonbb"))    xsec = 0.2151;

	//cross sections stolen from AN for TOP-18-009
	
        if(file.Contains("tZq"))   xsec = 0.09418;
	if(file.Contains("ttHToNonbb"))    xsec = 0.2151;

        // Calculated at 13 TeV in
        // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV
        // Higgs branching ratios from
        // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR
        if(file.Contains("ZH_HToBB_ZToLL_M-125"))   xsec = 0.883*Htobb*0.033658;
        if(file.Contains("ZH_HToBB_ZToNuNu_M-125")) xsec = 0.883*Htobb*0.2;
        if(file.Contains("WH_HToBB_WToLNu_M-125"))  xsec = 1.373*Htobb*(0.1071+0.1063+0.1138);
        if(file.Contains("ZH_HToBB_ZToNuNu_M125"))  xsec = 0.883*Htobb*0.2;
        if(file.Contains("WH_HToBB_WToLNu_M125"))   xsec = 1.373*Htobb*(0.1071+0.1063+0.1138);

        // Zgamma cross sections at 13 TeV
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
        if(file.Contains("DYJetsToLL") &&  
           file.Contains("amcatnlo"))      xsec = 6077.22;
        if(file.Contains("ZGTo2LG"))       xsec =  117.864;
        if(file.Contains("ZGToLLG"))       xsec =  117.864; // Other sample name
        if(file.Contains("ZZTo2L2Q"))      xsec =    3.22;
        if(file.Contains("ZZTo2L2Nu"))     xsec =    0.564;
        if(file.Contains("ZZTo4L"))        xsec =    1.256;
        if(file.Contains("TTTo2L2Nu"))     xsec =   87.31;
        if(file.Contains("TGJets"))        xsec =    2.967;
        if(file.Contains("LLAJJ"))         xsec =    0.1084; // from XSDB
        if(file.Contains("WZ_Tune"))       xsec =   47.13; //must come before WWZ for string matching
        if(file.Contains("WW_Tune"))       xsec =  116.7; // From arXiv:1408.5243
        if(file.Contains("WWW"))           xsec =    0.2086;
        if(file.Contains("WWZ"))           xsec =    0.1651;
        if(file.Contains("WZZ"))           xsec =    0.05565;
        if(file.Contains("ZZZ"))           xsec =    0.01398;
        if(file.Contains("WGGJets"))       xsec =    1.715;
        if(file.Contains("WWG"))           xsec =    0.2147; 
        if(file.Contains("ZGGJetsToLLGG")) xsec =    0.1699;
        // Zgamma signal
        if(file.Contains("GluGluHToZG"))   xsec = HToZG*ZToLL*48.58  ;// OG: 44.08;   J+M:48.58
        if(file.Contains("VBFHToZG"))      xsec = HToZG*ZToLL* 3.782 ;// OG:  3.779;  J+M: 3.782
        if(file.Contains("WplusH_HToZG"))  xsec = HToZG*ZToLL* 0.831 ;// OG:  0.8380  J+M: 0.831
        if(file.Contains("WminusH_HToZG")) xsec = HToZG*ZToLL* 0.527 ;// OG:  0.5313  J+M: 0.527
        if(file.Contains("ZH_HToZG"))      xsec = HToZG*ZToLL* 0.8839;// OG:  0.8824  J+M: 0.8839
        if(file.Contains("ttHToZG"))       xsec = HToZG*ZToLL* 0.5071;// OG:  0.5065  J+M: 0.5071
    } else {
        if(file.Contains("SMS-T1tttt_mGluino-1200_mLSP-800_Tune")) xsec = 0.0985;
        if(file.Contains("SMS-T1tttt_mGluino-2000_mLSP-100_Tune")) xsec = 0.00101;

        //  Cross-section taken from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
        // Alternative option: https://twiki.cern.ch/twiki/bin/view/Sandbox/FullNNLOcrossSections#Top_cross_section_for_13_TeV
        if(file.Contains("TTJets_Tune") || file.Contains("TT_"))  xsec = 815.96;

        if(file.Contains("TTJets_DiLept") || file.Contains("TTTo2L2Nu")) xsec = 85.66; // (3*0.108)^2*815.96
        if(file.Contains("TTJets_SingleLept")) xsec = 178.7; //(1- ((1-3*0.108)^2+(3*0.108)^2))*815.96*0.5 per half
        
        if(file.Contains("TTJets_DiLept_genMET-150")) xsec = 0.0676543*85.66; // filter_eff*(3*0.108)^2*815.96
        if(file.Contains("TTJets_SingleLept") && file.Contains("genMET-150")) xsec = 0.0568246*178.7; //filter_eff*(1- ((1-3*0.108)^2+(3*0.108)^2))*815.96*0.5 per half

        if(file.Contains("TTJets_DiLept_genMET-80")) xsec = 0.412882*85.66; // filter_eff*(3*0.108)^2*815.96
        if(file.Contains("TTJets_SingleLept") && file.Contains("genMET-80")) xsec = 0.293137*178.7; //filter_eff*(1- ((1-3*0.108)^2+(3*0.108)^2))*815.96*0.5 per half

        if(file.Contains("TTJets_HT")){//LO cross sections with k-factor of 1.625 already applied
          if(file.Contains("2500toInf")) xsec = 0.0015; // 0.0023234211;TTJets_HT-1200to2500
          if(file.Contains("1200to2500")) xsec = 0.1326; // 0.194972521;
          if(file.Contains("800to1200")) xsec = 0.7378; // 1.07722318;
          if(file.Contains("600to800")) xsec = 1.7966; // 2.61537118;
        }

        // from cross XSDB
        if(file.Contains("TTG")) xsec = 4.078;                
        if(file.Contains("TTTT_Tune")) xsec = 0.008213;
        
        // From https://cms-pdmv.cern.ch/mcm
        // k-factors from https://mangano.web.cern.ch/mangano/public/MECCA/samples_50ns_miniaod.txt
        // k-factors are ratio of https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
        // NLO/NNLO cross-sections to that of an inclusive sample in mcm at lower order (LO/NLO)

        if(file.Contains("WJetsToLNu_Tune")) xsec=20508.9*3; //NNLO from Lesya's summary table 

        //cross-section per slice based on inclusive sample, roughly 10% higher than 2016, less in extreme tail
        if(file.Contains("WJetsToLNu_HT-70To100"))  xsec = 0.0243795*20508.9*3;
        if(file.Contains("WJetsToLNu_HT-100To200"))  xsec = 0.0262096*20508.9*3;
        if(file.Contains("WJetsToLNu_HT-200To400"))  xsec = 0.00772818*20508.9*3;
        if(file.Contains("WJetsToLNu_HT-400To600"))  xsec = 0.00109366*20508.9*3;
        if(file.Contains("WJetsToLNu_HT-600To800"))  xsec = 0.000272388*20508.9*3;
        if(file.Contains("WJetsToLNu_HT-800To1200"))  xsec = 0.000122233*20508.9*3;
        if(file.Contains("WJetsToLNu_HT-1200To2500"))  xsec = 2.71060e-5*20508.9*3;
        if(file.Contains("WJetsToLNu_HT-2500ToInf"))  xsec = 3.94174e-07*20508.9*3;

        if(file.Contains("QCD_HT100to200_Tune")) xsec = 23700000;
        if(file.Contains("QCD_HT200to300_Tune"))   xsec = 1547000;
        if(file.Contains("QCD_HT300to500_Tune"))   xsec = 322600;
        if(file.Contains("QCD_HT500to700_Tune"))   xsec = 29980;
        if(file.Contains("QCD_HT700to1000_Tune"))  xsec = 6334;
        if(file.Contains("QCD_HT1000to1500_Tune")) xsec = 1088;
        if(file.Contains("QCD_HT1500to2000_Tune")) xsec = 99.11;
        if(file.Contains("QCD_HT2000toInf_Tune"))  xsec = 20.23;

        // Cross sections from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec
        // multiplied by BF(W->mu,e,tau) = 0.324
        if (file.Contains("ST_s-channel_4f_leptonDecays"))     xsec = 3.34;
        if (file.Contains("ST_t-channel_antitop_4f_inclusiveDecays") ||
            file.Contains("ST_t-channel_antitop_4f_InclusiveDecays")) xsec = 80.95;
        if (file.Contains("ST_t-channel_top_4f_inclusiveDecays") || 
            file.Contains("ST_t-channel_top_4f_InclusiveDecays")) xsec = 136.02;
        if (file.Contains("ST_tW_antitop_5f_NoFullyHadronicDecays"))     xsec = 35.85*0.543;
        if (file.Contains("ST_tW_top_5f_NoFullyHadronicDecays"))     xsec = 35.85*0.543; 

        if(file.Contains("DYJetsToLL_M-50_Tune"))     xsec = 2075.14*3;

        if(file.Contains("DYJetsToLL_M-50_HT-70to100"))    xsec = 0.0275032*2075.14*3;
        if(file.Contains("DYJetsToLL_M-50_HT-100to200"))    xsec = 0.0302083*2075.14*3;
        if(file.Contains("DYJetsToLL_M-50_HT-200to400"))    xsec = 0.00907651*2075.14*3;
        if(file.Contains("DYJetsToLL_M-50_HT-400to600"))    xsec = 0.00129238*2075.14*3;
        if(file.Contains("DYJetsToLL_M-50_HT-600to800"))    xsec = 0.000316039*2075.14*3;
        if(file.Contains("DYJetsToLL_M-50_HT-800to1200"))    xsec = 0.000137432*2075.14*3;
        if(file.Contains("DYJetsToLL_M-50_HT-1200to2500"))    xsec = 3.09368e-05*2075.14*3;
        if(file.Contains("DYJetsToLL_M-50_HT-2500toInf"))    xsec = 4.39860e-07*2075.14*3;

        // k-factor from DYJets 1.165
        if(file.Contains("ZJetsToNuNu_HT-100To200"))  xsec = 302.8*1.165;
        if(file.Contains("ZJetsToNuNu_HT-200To400"))  xsec = 92.59*1.165;
        if(file.Contains("ZJetsToNuNu_HT-400To600"))  xsec = 13.18*1.165;
        if(file.Contains("ZJetsToNuNu_HT-600To800"))  xsec = 3.257*1.165;
        if(file.Contains("ZJetsToNuNu_HT-800To1200"))  xsec = 1.49*1.165;
        if(file.Contains("ZJetsToNuNu_HT-1200To2500"))  xsec = 0.3419*1.165;
        if(file.Contains("ZJetsToNuNu_HT-2500ToInf"))  xsec = 0.005146*1.165;

        if(file.Contains("TTZToQQ"))                xsec = 0.5104;
        if(file.Contains("TTZToLLNuNu_M-10"))       xsec = 0.2432;
        if(file.Contains("TTWJetsToQQ"))            xsec = 0.4316;
        if(file.Contains("TTWJetsToLNu"))           xsec = 0.2149;
       
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
        if(file.Contains("WWTo2L2Nu"))   xsec = 12.178; //NNLO
        if(file.Contains("WWToLNuQQ"))   xsec = 49.997; //NNLO
        if(file.Contains("ttHTobb_M125"))   xsec = 0.2934;

        if(file.Contains("WZTo1L3Nu"))   xsec = 3.05;
        if(file.Contains("WZTo1L1Nu2Q"))   xsec = 10.73;
        if(file.Contains("WZTo2L2Q"))   xsec = 5.606;
        if(file.Contains("WZTo3LNu"))   xsec = 4.42965;
        if(file.Contains("VVTo2L2Nu"))   xsec = 11.95;
        if(file.Contains("ZZ_Tune"))   xsec = 16.523; //from twiki NLO

        // Calculated at 13 TeV in
        // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV
        // Higgs branching ratios from
        // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR
        if(file.Contains("ZH_HToBB_ZToLL_M-125"))      xsec = 0.883*Htobb*0.033658;
        if(file.Contains("ZH_HToBB_ZToNuNu_M-125"))    xsec = 0.883*Htobb*0.2;
        if(file.Contains("WH_HToBB_WToLNu_M-125"))     xsec = 1.373*Htobb*(0.1071+0.1063+0.1138);
        if(file.Contains("ZH_HToBB_ZToNuNu_M125"))    xsec = 0.883*Htobb*0.2;
        if(file.Contains("WHiggs0PMToBB"))     xsec = 1.373*Htobb*(0.1071+0.1063+0.1138); //using same as 2016

	//cross sections stolen from AN for TOP-18-009
	
        if(file.Contains("tZq"))   xsec = 0.09418;
	if(file.Contains("ttHToNonbb"))    xsec = 0.2151;
    }
    if(xsec<=0) std::cout<<"ERROR:: Cross section not found for "<<file<<std::endl;

    return xsec;
  }

  float fractionNegWeights(const TString &file){
    float fneg(0.);

    // ttH, ttZ, ttW, ttGamma, tZq
    if(file.Contains("ttHJetTobb_M125_13TeV_amcatnloFXFX"))     fneg = 0.3515;
    if(file.Contains("TTZToQQ"))                                fneg = 0.2657;
    if(file.Contains("TTZToLLNuNu_M-10"))                       fneg = 0.2676;
    if(file.Contains("TTWJetsToQQ"))                            fneg = 0.2412;
    if(file.Contains("TTWJetsToLNu"))                           fneg = 0.2433;
    if(file.Contains("TTG"))                                    fneg = 0.342; // from MCM

    if(file.Contains("TTTT_TuneCUETP8M1_13TeV-amcatnlo"))       fneg = 0.41; // from MCM
    if(file.Contains("VVTo2L2Nu_13TeV_amcatnloFXFX"))       fneg = 0.20; // from MCM
    if(file.Contains("TTJets_Mtt-1000toInf"))                   fneg = 0.376996;
    if(file.Contains("tZq"))					fneg = 0.3675;

    // Single top
    if (file.Contains("ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8")) fneg = 0.1884;

    return fneg;
  }


  void higgsino2DCrossSection(int hig_mass, double &xsec, double &xsec_unc) {
    if(hig_mass ==127) { xsec = .5824*.5824*1.44725; xsec_unc = 0.0395277; return;}
    else if(hig_mass ==150) { xsec = .5824*.5824*0.71514; xsec_unc = 0.0421496; return;}
    else if(hig_mass ==175) { xsec = .5824*.5824*0.419059; xsec_unc = 0.0453279; return;}
    else if(hig_mass ==200) { xsec = .5824*.5824*0.244213; xsec_unc = 0.047925; return;}
    else if(hig_mass ==225) { xsec = .5824*.5824*0.156286; xsec_unc = 0.0502876; return;}
    else if(hig_mass ==250) { xsec = .5824*.5824*0.104252; xsec_unc = 0.0526169; return;}
    else if(hig_mass ==275) { xsec = .5824*.5824*0.0719125; xsec_unc = 0.0549666; return;}
    else if(hig_mass ==300) { xsec = .5824*.5824*0.0509994; xsec_unc = 0.0572762; return;}
    else if(hig_mass ==325) { xsec = .5824*.5824*0.0369715; xsec_unc = 0.0590317; return;}
    else if(hig_mass ==350) { xsec = .5824*.5824*0.0273286; xsec_unc = 0.0607766; return;}
    else if(hig_mass ==375) { xsec = .5824*.5824*0.0205429; xsec_unc = 0.0625031; return;}
    else if(hig_mass ==400) { xsec = .5824*.5824*0.0156691; xsec_unc = 0.0642085; return;}
    else if(hig_mass ==425) { xsec = .5824*.5824*0.0120965; xsec_unc = 0.0657801; return;}
    else if(hig_mass ==450) { xsec = .5824*.5824*0.00944017; xsec_unc = 0.0674544; return;}
    else if(hig_mass ==475) { xsec = .5824*.5824*0.00743587; xsec_unc = 0.0686033; return;}
    else if(hig_mass ==500) { xsec = .5824*.5824*0.00590757; xsec_unc = 0.0699909; return;}
    else if(hig_mass ==526) { xsec = .5824*.5824*0.00469101; xsec_unc = 0.0713704; return;}
    else if(hig_mass ==550) { xsec = .5824*.5824*0.0038167; xsec_unc = 0.0722834; return;}
    else if(hig_mass ==576) { xsec = .5824*.5824*0.003073; xsec_unc = 0.0739957; return;}
    else if(hig_mass ==600) { xsec = .5824*.5824*0.00253015; xsec_unc = 0.0754291; return;}
    else if(hig_mass ==626) { xsec = .5824*.5824*0.00206136; xsec_unc = 0.0763466; return;}
    else if(hig_mass ==650) { xsec = .5824*.5824*0.00171418; xsec_unc = 0.0775695; return;}
    else if(hig_mass ==676) { xsec = .5824*.5824*0.00140934; xsec_unc = 0.0783375; return;}
    else if(hig_mass ==700) { xsec = .5824*.5824*0.00118113; xsec_unc = 0.0796388; return;}
    else if(hig_mass ==726) { xsec = .5824*.5824*0.000979349; xsec_unc = 0.0809883; return;}
    else if(hig_mass ==750) { xsec = .5824*.5824*0.000826366; xsec_unc = 0.081879; return;}
    else if(hig_mass ==776) { xsec = .5824*.5824*0.000690208; xsec_unc = 0.0842049; return;}
    else if(hig_mass ==800) { xsec = .5824*.5824*0.000586211; xsec_unc = 0.0862527; return;}
    else if(hig_mass ==826) { xsec = .5824*.5824*0.00049277; xsec_unc = 0.0864444; return;}
    else if(hig_mass ==850) { xsec = .5824*.5824*0.000420556; xsec_unc = 0.085742; return;}
    else if(hig_mass ==876) { xsec = .5824*.5824*0.000358734; xsec_unc = 0.0889174; return;}
    else if(hig_mass ==900) { xsec = .5824*.5824*0.000305935; xsec_unc = 0.0912439; return;}
    else if(hig_mass ==926) { xsec = .5824*.5824*0.000260948; xsec_unc = 0.091372; return;}
    else if(hig_mass ==950) { xsec = .5824*.5824*0.00022285; xsec_unc = 0.0919538; return;}
    else if(hig_mass ==976) { xsec = .5824*.5824*0.000189681; xsec_unc = 0.0938108; return;}
    else if(hig_mass ==1000) { xsec = .5824*.5824*0.00016428; xsec_unc = 0.0954285; return;}
    else if(hig_mass ==1024) { xsec = .5824*.5824*0.000142206; xsec_unc = 0.0957231; return;}
    else if(hig_mass ==1052) { xsec = .5824*.5824*0.000120971; xsec_unc = 0.0968997; return;}
    else if(hig_mass ==1076) { xsec = .5824*.5824*0.000105301; xsec_unc = 0.0979041; return;}
    else if(hig_mass ==1100) { xsec = .5824*.5824*9.12469e-05; xsec_unc = 0.0964142; return;}
    else if(hig_mass ==1124) { xsec = .5824*.5824*7.9765e-05; xsec_unc = 0.099902; return;}
    else if(hig_mass ==1152) { xsec = .5824*.5824*6.78234e-05; xsec_unc = 0.101061; return;}
    else if(hig_mass ==1176) { xsec = .5824*.5824*5.9016e-05; xsec_unc = 0.102051; return;}
    else if(hig_mass ==1200) { xsec = .5824*.5824*5.16263e-05; xsec_unc = 0.102499; return;}
    else if(hig_mass ==1224) { xsec = .5824*.5824*4.5147e-05; xsec_unc = 0.10403; return;}
    else if(hig_mass ==1252) { xsec = .5824*.5824*3.88343e-05; xsec_unc = 0.105206; return;}
    else if(hig_mass ==1276) { xsec = .5824*.5824*3.41304e-05; xsec_unc = 0.10619; return;}
    else if(hig_mass ==1300) { xsec = .5824*.5824*2.99353e-05; xsec_unc = 0.10783; return;}
    else if(hig_mass ==1324) { xsec = .5824*.5824*2.63637e-05; xsec_unc = 0.108024; return;}
    else if(hig_mass ==1352) { xsec = .5824*.5824*2.26779e-05; xsec_unc = 0.109016; return;}
    else if(hig_mass ==1376) { xsec = .5824*.5824*1.99318e-05; xsec_unc = 0.109822; return;}
    else if(hig_mass ==1400) { xsec = .5824*.5824*1.75031e-05; xsec_unc = 0.111631; return;}
    else if(hig_mass ==1424) { xsec = .5824*.5824*1.53974e-05; xsec_unc = 0.111417; return;}
    else if(hig_mass ==1452) { xsec = .5824*.5824*1.3245e-05; xsec_unc = 0.112313; return;}
    else if(hig_mass ==1476) { xsec = .5824*.5824*1.16416e-05; xsec_unc = 0.113058; return;}
    else{ xsec = 0; xsec_unc = 0;}
  }

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
