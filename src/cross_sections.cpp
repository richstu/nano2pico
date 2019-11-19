// Lists the MC cross sections

#include <iostream>
#include "cross_sections.hpp"

using namespace std;

namespace xsec{

  float crossSection(const TString &file, bool is2016){
    float xsec(-999999.), Htobb(0.5824), HToZG(0.001541);

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
        if(file.Contains("WJetsToLNu_HT-70To100"))  xsec = 1372.*1.21; 
        if(file.Contains("WJetsToLNu_HT-100To200"))  xsec = 1347.*1.21; 
        if(file.Contains("WJetsToLNu_HT-200To400"))  xsec = 360.*1.21;
        if(file.Contains("WJetsToLNu_HT-400To600"))  xsec = 48.98*1.21;
        if(file.Contains("WJetsToLNu_HT-600ToInf"))  xsec = 18.77*1.21;
        if(file.Contains("WJetsToLNu_HT-600To800"))  xsec = 12.05*1.21;
        if(file.Contains("WJetsToLNu_HT-800To1200"))  xsec = 5.501*1.21;
        if(file.Contains("WJetsToLNu_HT-1200To2500"))  xsec = 1.329*1.21;
        if(file.Contains("WJetsToLNu_HT-2500ToInf"))  xsec = 0.03216*1.21;

        //updated 02-2019 with XSDB
        if(file.Contains("QCD_HT100to200_Tune")) xsec = 28060000;
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
        if(file.Contains("DYJetsToLL_M-50_HT-100to200"))    xsec = 139.4*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-200to400"))    xsec = 42.75*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-400to600"))    xsec = 5.497*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-600to800"))    xsec = 1.363*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-800to1200"))    xsec = 0.6759*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-1200to2500"))    xsec = 0.116*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-2500toInf"))    xsec = 0.002592*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-600toInf"))    xsec = 2.21*1.23;

        if(file.Contains("ZJetsToNuNu_HT-100To200"))  xsec = 280.35*1.27;
        if(file.Contains("ZJetsToNuNu_HT-200To400"))  xsec = 77.67*1.27;
        if(file.Contains("ZJetsToNuNu_HT-400To600"))  xsec = 10.73*1.27;
        if(file.Contains("ZJetsToNuNu_HT-600To800"))  xsec = 2.536*1.27;
        if(file.Contains("ZJetsToNuNu_HT-800To1200"))  xsec = 1.161*1.27;
        if(file.Contains("ZJetsToNuNu_HT-1200To2500"))  xsec = 0.2824*1.27;
        if(file.Contains("ZJetsToNuNu_HT-2500ToInf"))  xsec = 0.006459*1.27;
        if(file.Contains("ZJetsToNuNu_HT-600ToInf"))  xsec = 3.986*1.27;

        if(file.Contains("TTZToQQ"))                xsec = 0.5297;
        if(file.Contains("TTZToLLNuNu_M-10"))       xsec = 0.2529;
        if(file.Contains("TTWJetsToQQ"))            xsec = 0.4062;
        if(file.Contains("TTWJetsToLNu"))           xsec = 0.2043;
       
        //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
        if(file.Contains("WWTo2L2Nu"))   xsec = 12.178; //NNLO
        if(file.Contains("WWToLNuQQ"))   xsec = 49.997; //NNLO
        if(file.Contains("ttHTobb_M125"))   xsec = 0.2934;

        if(file.Contains("WZTo1L3Nu"))   xsec = 3.05;
        if(file.Contains("WZTo1L1Nu2Q"))   xsec = 10.96;
        if(file.Contains("WZTo2L2Q"))   xsec = 5.595;
        if(file.Contains("WZTo3LNu"))   xsec = 4.42965;
        if(file.Contains("VVTo2L2Nu"))   xsec = 11.95;
        if(file.Contains("ZZ_Tune"))   xsec = 16.523;

	//cross sections stolen from AN for TOP-18-009
	
        if(file.Contains("tZq"))   xsec = 0.09418;
	if(file.Contains("ttHToNonbb"))    xsec = 0.2151;

        // Calculated at 13 TeV in
        // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV
        // Higgs branching ratios from
        // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR
        if(file.Contains("ZH_HToBB_ZToLL_M-125"))      xsec = 0.883*Htobb*0.033658;
        if(file.Contains("ZH_HToBB_ZToNuNu_M-125"))    xsec = 0.883*Htobb*0.2;
        if(file.Contains("WH_HToBB_WToLNu_M-125"))     xsec = 1.373*Htobb*(0.1071+0.1063+0.1138);
        if(file.Contains("ZH_HToBB_ZToNuNu_M125"))    xsec = 0.883*Htobb*0.2;
        if(file.Contains("WH_HToBB_WToLNu_M125"))     xsec = 1.373*Htobb*(0.1071+0.1063+0.1138);

        // Zgamma cross sections at 13 TeV
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
        if(file.Contains("DYJetsToLL") &&  
           file.Contains("amcatnlo"))      xsec = 6077.22;
        if(file.Contains("ZGTo2LG"))       xsec =  117.864;
        if(file.Contains("ZZTo2L2Q"))      xsec =    3.22;
        if(file.Contains("ZZTo2L2Nu"))     xsec =    0.564;
        if(file.Contains("ZZTo4L"))        xsec =    1.256;
        if(file.Contains("TTTo2L2Nu"))     xsec =   87.31;
        if(file.Contains("WWW"))           xsec =    0.2086;
        if(file.Contains("WWZ"))           xsec =    0.1651;
        if(file.Contains("WZZ"))           xsec =    0.05565;
        if(file.Contains("ZZZ"))           xsec =    0.01398;
        if(file.Contains("WGGJets"))       xsec =    1.715;
        if(file.Contains("WWG"))           xsec =    0.2147; 
        if(file.Contains("ZGGJetsToLLGG")) xsec =    0.1699;
        // Zgamma signal
        if(file.Contains("GluGluHToZG"))   xsec = HToZG*44.08;
        if(file.Contains("VBFHToZG"))      xsec = HToZG* 3.779;
        if(file.Contains("WplusH_HToZG"))  xsec = HToZG* 0.8380;
        if(file.Contains("WminusH_HToZG")) xsec = HToZG* 0.5313;
        if(file.Contains("ZH_HToZG"))      xsec = HToZG* 0.8824;
        if(file.Contains("ttHToZG"))       xsec = HToZG* 0.5065;
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

        // from cross XSDB
        if(file.Contains("TTG")) xsec = 4.078;                
        if(file.Contains("TTTT_Tune")) xsec = 0.008213;
        
        // From https://cms-pdmv.cern.ch/mcm
        // k-factors from https://mangano.web.cern.ch/mangano/public/MECCA/samples_50ns_miniaod.txt
        // k-factors are ratio of https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
        // NLO/NNLO cross-sections to that of an inclusive sample in mcm at lower order (LO/NLO)

        if(file.Contains("WJetsToLNu_Tune")) xsec=20508.9*3; //NNLO from Lesya's summary table 

        //cross-section per slice based on inclusive sample, roughly 10% higher than 2016, less in extreme tail
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
        if(file.Contains("WH_HToBB_WToLNu_M125"))     xsec = 1.373*Htobb*(0.1071+0.1063+0.1138);

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


  void higgsinoCrossSection(int hig_mass, double &xsec, double &xsec_unc) {
        if(hig_mass == 127) { xsec = .5824*.5824*7602.2/1000; xsec_unc = 0.0393921; return;}
        else if(hig_mass == 150) { xsec = .5824*.5824*3832.31/1000; xsec_unc = 0.0413612; return;}
        else if(hig_mass == 175) { xsec = .5824*.5824*2267.94/1000; xsec_unc = 0.044299;  return;}
        else if(hig_mass == 200) { xsec = .5824*.5824*1335.62/1000; xsec_unc = 0.0474362; return;}
        else if(hig_mass == 225) { xsec = .5824*.5824*860.597/1000; xsec_unc = 0.0504217; return;}
        else if(hig_mass == 250) { xsec = .5824*.5824*577.314/1000; xsec_unc = 0.0532731; return;}
        else if(hig_mass == 275) { xsec = .5824*.5824*400.107/1000; xsec_unc = 0.0560232; return;}
        else if(hig_mass == 300) { xsec = .5824*.5824*284.855/1000; xsec_unc = 0.0586867; return;}
        else if(hig_mass == 325) { xsec = .5824*.5824* 207.36/1000; xsec_unc = 0.0613554; return;}
        else if(hig_mass == 350) { xsec = .5824*.5824*153.841/1000; xsec_unc = 0.0640598; return;}
        else if(hig_mass == 375) { xsec = .5824*.5824*116.006/1000; xsec_unc = 0.066892;  return;}
        else if(hig_mass == 400) { xsec = .5824*.5824*88.7325/1000; xsec_unc = 0.0697517; return;}
        else if(hig_mass == 425) { xsec = .5824*.5824*68.6963/1000; xsec_unc = 0.0723531; return;}
        else if(hig_mass == 450) { xsec = .5824*.5824*53.7702/1000; xsec_unc = 0.0748325; return;}
        else if(hig_mass == 475) { xsec = .5824*.5824*42.4699/1000; xsec_unc = 0.0775146; return;}
        else if(hig_mass == 500) { xsec = .5824*.5824*33.8387/1000; xsec_unc = 0.0802572; return;}
        else if(hig_mass == 525) { xsec = .5824*.5824*27.1867/1000; xsec_unc = 0.0825803; return;}
        else if(hig_mass == 550) { xsec = .5824*.5824*21.9868/1000; xsec_unc = 0.0849278; return;}
        else if(hig_mass == 575) { xsec = .5824*.5824*17.9062/1000; xsec_unc = 0.087561;  return;}
        else if(hig_mass == 600) { xsec = .5824*.5824*14.6677/1000; xsec_unc = 0.0900693; return;}
        else if(hig_mass == 625) { xsec = .5824*.5824* 12.062/1000; xsec_unc = 0.091959;  return;}
        else if(hig_mass == 650) { xsec = .5824*.5824*9.96406/1000; xsec_unc = 0.094065;  return;}
        else if(hig_mass == 675) { xsec = .5824*.5824*8.28246/1000; xsec_unc = 0.0957436; return;}
        else if(hig_mass == 700) { xsec = .5824*.5824*6.89981/1000; xsec_unc = 0.0982894; return;}
        else if(hig_mass == 725) { xsec = .5824*.5824*5.78355/1000; xsec_unc = 0.0999915; return;}
        else if(hig_mass == 750) { xsec = .5824*.5824* 4.8731/1000; xsec_unc = 0.101211;  return;}
        else if(hig_mass == 775) { xsec = .5824*.5824*4.09781/1000; xsec_unc = 0.104646;  return;}
        else if(hig_mass == 800) { xsec = .5824*.5824*3.46143/1000; xsec_unc = 0.107618;  return;}
        else if(hig_mass == 825) { xsec = .5824*.5824* 2.9337/1000; xsec_unc = 0.108353;  return;}
        else if(hig_mass == 850) { xsec = .5824*.5824* 2.4923/1000; xsec_unc = 0.110016;  return;}
        else if(hig_mass == 875) { xsec = .5824*.5824*2.13679/1000; xsec_unc = 0.112636;  return;}
        else if(hig_mass == 900) { xsec = .5824*.5824*1.80616/1000; xsec_unc = 0.1134;    return;}
        else if(hig_mass == 925) { xsec = .5824*.5824*1.55453/1000; xsec_unc = 0.116949;  return;}
        else if(hig_mass == 950) { xsec = .5824*.5824*1.32692/1000; xsec_unc = 0.117027;  return;}
        else if(hig_mass == 975) { xsec = .5824*.5824*1.12975/1000; xsec_unc = 0.121244;  return;}
        else if(hig_mass ==1000) { xsec = .5824*.5824*0.968853/1000; xsec_unc = 0.126209; return;}
        else{ xsec = 0; xsec_unc = 0;}
  }
}
