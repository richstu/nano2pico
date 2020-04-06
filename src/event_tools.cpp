#include "event_tools.hpp"

#include "utilities.hpp"
#include "hig_trig_eff.hpp"

using namespace std;

EventTools::EventTools(const string &name_, int year_):
  name(name_),
  year(year_),
  isTTJets_LO_MET(false),
  isTTJets_LO_Incl(false),
  isWJets_LO(false),
  isDYJets_LO(false){

  if(Contains(name, "TTJets_") && Contains(name, "genMET-") && Contains(name, "madgraphMLM")) 
    isTTJets_LO_MET = true;

  if(Contains(name, "TTJets_") && !Contains(name, "TTJets_HT") && !Contains(name, "genMET-") &&Contains(name, "madgraphMLM")) 
    isTTJets_LO_Incl = true;

  if(Contains(name, "WJetsToLNu_Tune")  && Contains(name,"madgraphMLM"))
    isWJets_LO = true;
  
  if(Contains(name, "DYJetsToLL_M-50_Tune")  && Contains(name,"madgraphMLM"))
    isDYJets_LO = true;
}

EventTools::~EventTools(){
}

void EventTools::WriteStitch(nano_tree &nano, pico_tree &pico){
  pico.out_stitch_htmet() = pico.out_stitch() = pico.out_stitch_ht() = true;
  if(isTTJets_LO_Incl) {
    if (nano.LHE_HTIncoming()>600) 
      pico.out_stitch_htmet() = pico.out_stitch_ht() = false;
    if(year==2018 && nano.GenMET_pt()>80)
      pico.out_stitch_htmet() = pico.out_stitch() = false;
    else if (nano.GenMET_pt()>150) 
      pico.out_stitch_htmet() = pico.out_stitch() = false;
  }

  if (isTTJets_LO_MET && nano.LHE_HTIncoming()>600) 
      pico.out_stitch_htmet() = pico.out_stitch_ht() = false;

  if(isDYJets_LO  && nano.LHE_HT()>70) 
    pico.out_stitch_htmet() = pico.out_stitch_ht() = pico.out_stitch() = false;
  
  if(isWJets_LO  && nano.LHE_HT()>70) 
    pico.out_stitch_htmet() = pico.out_stitch_ht() = pico.out_stitch() = false;
  return;
}


void EventTools::WriteDataQualityFilters(nano_tree& nano, pico_tree& pico, vector<int> sig_jet_nano_idx,
                                         float min_jet_pt, float max_jet_eta, bool isData, bool isFastsim){
  // jet quality filter
  pico.out_pass_jets() = true;
  if (isFastsim) {
    // Fastsim: veto if certain central jets have no matching GenJet as per SUSY recommendation:
    // https://twiki.cern.ch/twiki/bin/view/CMS/SUSRecommendations18#Cleaning_up_of_fastsim_jets_from
    for(int ijet(0); ijet<nano.nJet(); ++ijet){
      if(nano.Jet_pt()[ijet] > 20 && fabs(nano.Jet_eta()[ijet])<=2.5 && nano.Jet_chHEF()[ijet] < 0.1) {
        bool found_match = false;
        for(int igenjet(0); igenjet<nano.nGenJet(); ++igenjet){
          if (dR(nano.Jet_eta()[ijet], nano.GenJet_eta()[igenjet], nano.Jet_phi()[ijet], nano.GenJet_phi()[igenjet])<=0.3) {
            found_match = true;
            break;
          }
        }
        if (!found_match) {
          pico.out_pass_jets() = false;
          break;
        }
      }
    }
  } else { // Fullsim: require just loosest possible ID for now (for all jets, not just central!)
    for(int ijet(0); ijet<nano.nJet(); ++ijet){
      if (nano.Jet_pt()[ijet] > min_jet_pt && nano.Jet_jetId()[ijet] < 1) 
        pico.out_pass_jets() = false;
    } 
  }

  // RA2b filters  
  pico.out_pass_muon_jet() = true; 
  for (auto &idx: sig_jet_nano_idx){
    // if (abs(nano.Jet_eta()[idx])>2.4) continue; -> already enforced in signal jet selection
    // if is overlapping with lepton -> already enforced in signal jet selection
    if (nano.Jet_pt()[idx]<=200.) continue;
    if (nano.Jet_muEF()[idx]<=0.5) continue;
    if (DeltaPhi(nano.Jet_phi()[idx],nano.MET_phi())<(3.14159-0.4)) continue;
    pico.out_pass_muon_jet() = false;
    break;
  }

  pico.out_pass_low_neutral_jet() = true;
  for(int ijet(0); ijet<nano.nJet(); ++ijet){  
    if (nano.Jet_pt()[ijet]<=min_jet_pt || fabs(nano.Jet_eta()[ijet])>max_jet_eta) continue;
    if (nano.Jet_neEmEF()[ijet] <0.03 && DeltaPhi(nano.Jet_phi()[ijet], pico.out_mht_phi())>(3.14159-0.4))
      pico.out_pass_low_neutral_jet() = false;
    break; //only apply to leading jet that passes pt and eta
  }

  pico.out_pass_htratio_dphi_tight() = true;
  float htratio = pico.out_ht5()/pico.out_ht();
  for(int ijet(0); ijet<nano.nJet(); ++ijet){  
    if (nano.Jet_pt()[ijet]<=min_jet_pt || fabs(nano.Jet_eta()[ijet])>max_jet_eta) continue;
    if (htratio > 1.2 && DeltaPhi(nano.Jet_phi()[ijet], pico.out_mht_phi()) > (5.3*htratio - 4.78)) 
      pico.out_pass_htratio_dphi_tight() = false;
    break; //only apply to leading jet that passes pt and eta
  }

  // filters directly from Nano
  pico.out_pass_hbhe() = nano.Flag_HBHENoiseFilter();
  pico.out_pass_hbheiso() = nano.Flag_HBHENoiseIsoFilter();
  pico.out_pass_goodv() = nano.Flag_goodVertices();
  pico.out_pass_cschalo_tight() = nano.Flag_globalSuperTightHalo2016Filter();
  pico.out_pass_eebadsc() = nano.Flag_eeBadScFilter();
  pico.out_pass_ecaldeadcell() = nano.Flag_EcalDeadCellTriggerPrimitiveFilter();
  if (year==2016) {
    pico.out_pass_badcalib() = true;
    pico.out_pass_badpfmu() = nano.Flag_BadPFMuonSummer16Filter();
    pico.out_pass_badchhad() = nano.Flag_BadChargedCandidateSummer16Filter();
  } else {
    pico.out_pass_badcalib() = nano.Flag_ecalBadCalibFilterV2();
    pico.out_pass_badpfmu() = nano.Flag_BadPFMuonFilter();
    pico.out_pass_badchhad() = nano.Flag_BadChargedCandidateFilter();
  }
  pico.out_pass_mubadtrk() = nano.Flag_muonBadTrackFilter();

  // Combined pass variable, as recommended here:
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Analysis_Recommendations_for_ana
  pico.out_pass() = pico.out_pass_muon_jet() && pico.out_pass_badpfmu() && 
                    pico.out_met()/pico.out_met_calo()<5 &&
                    pico.out_pass_goodv() &&
                    pico.out_pass_hbhe() && pico.out_pass_hbheiso() && 
                    pico.out_pass_ecaldeadcell() && pico.out_pass_badcalib() &&
                    pico.out_pass_jets();

  if (!isFastsim) {
    pico.out_pass() = pico.out_pass() && pico.out_pass_cschalo_tight();
    if (isData) 
      pico.out_pass() = pico.out_pass() && pico.out_pass_eebadsc();
  }
  
  // Combined pass RA2b-like variable
  // https://github.com/rpatelCERN/boostedHiggsPlusMET/blob/Higgsino/src/definitions.cc#L1137
  pico.out_pass_boosted() = pico.out_pass_hbhe() && 
                         pico.out_pass_hbheiso() && 
                         pico.out_pass_eebadsc() && 
                         // pico.out_pass_ecaldeadcell() && 
                         pico.out_pass_goodv() &&
                         pico.out_met()/pico.out_met_calo() < 5. &&
                         pico.out_pass_badpfmu() && 
                         pico.out_pass_cschalo_tight() && 
                         pico.out_pass_low_neutral_jet() && 
                         pico.out_pass_htratio_dphi_tight() && 
                         pico.out_pass_jets();

  return;
}

void EventTools::CopyTriggerDecisions(nano_tree& nano, pico_tree& pico){
  //single lepton triggers
  pico.out_HLT_IsoMu24() = nano.HLT_IsoMu24();
  pico.out_HLT_IsoMu27() = nano.HLT_IsoMu27();
  pico.out_HLT_IsoTkMu24() = nano.HLT_IsoTkMu24();
  pico.out_HLT_Mu50() = nano.HLT_Mu50();
  pico.out_HLT_Ele27_WPTight_Gsf() = nano.HLT_Ele27_WPTight_Gsf();
  pico.out_HLT_Ele35_WPTight_Gsf() = nano.HLT_Ele35_WPTight_Gsf();
  pico.out_HLT_Ele115_CaloIdVT_GsfTrkIdT() = nano.HLT_Ele115_CaloIdVT_GsfTrkIdT();
  pico.out_HLT_Mu8() = nano.HLT_Mu8();
  pico.out_HLT_Mu17() = nano.HLT_Mu17();
  pico.out_HLT_Ele8_CaloIdM_TrackIdM_PFJet30() = nano.HLT_Ele8_CaloIdM_TrackIdM_PFJet30();
  pico.out_HLT_Ele17_CaloIdM_TrackIdM_PFJet30() = nano.HLT_Ele17_CaloIdM_TrackIdM_PFJet30();

  // MET triggers
  pico.out_HLT_PFMET110_PFMHT110_IDTight() = nano.HLT_PFMET110_PFMHT110_IDTight();
  pico.out_HLT_PFMET120_PFMHT120_IDTight() = nano.HLT_PFMET120_PFMHT120_IDTight();
  pico.out_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight() = nano.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight();
  pico.out_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight() = nano.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight();
  pico.out_HLT_PFMET120_PFMHT120_IDTight_PFHT60() = nano.HLT_PFMET120_PFMHT120_IDTight_PFHT60();
  pico.out_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60() = nano.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60();

  // Jet trigger
  pico.out_HLT_PFJet500() = nano.HLT_PFJet500();

  // ZGamma triggers
  pico.out_HLT_Mu17_Photon30_IsoCaloId()               = nano.HLT_Mu17_Photon30_IsoCaloId();
  pico.out_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() = nano.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ();
  pico.out_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL()    = nano.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL();
  pico.out_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL()          = nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL();
  pico.out_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL()        = nano.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL();
  pico.out_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ()       = nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ();
  pico.out_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ()     = nano.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ();
  pico.out_HLT_Ele27_eta2p1_WPTight_Gsf()              = nano.HLT_Ele27_eta2p1_WPTight_Gsf();
  pico.out_HLT_Photon175()                             = nano.HLT_Photon175();
  return;
}

int EventTools::GetEventType(){
  int sample = -999, category = -9, bin = -99;
  if(Contains(name, "Run201")){ sample = 0;
    if(Contains(name, "SingleElectron") || Contains(name, "EGamma")){ category = 0;
    }else if(Contains(name, "SingleMuon")){ category = 1;
    }else if(Contains(name, "DoubleEG")){   category = 2;
    }else if(Contains(name, "DoubleMuon")){ category = 3;
    }else if(Contains(name, "MET")){        category = 4;
    }else if(Contains(name, "HTMHT")){      category = 5;
    }else if(Contains(name, "JetHT")){      category = 6;
    }
    auto pos = name.find("Run201")+7;
    if(pos < name.size()
       && isalpha(name.at(pos))
       && pos-1 < name.size()
       && isdigit(name.at(pos-1))){
      int run = toupper(name.at(pos))-65;
      int year_ = name.at(pos-1)-48;
      bin = (year_-5)*26+run+1;
      //2015A=1, ..., 2015D=4, ..., 2016A=27, 2016B=28, 2016C=29, ...
      if(bin > 99) bin = 0;//Sanity check
    }
  }else if((Contains(name, "TTJets") || Contains(name, "TT_")) && !Contains(name, "TTTT_")){ sample = 1;
    if(Contains(name, "TTJets_Tune")){ category = 0; bin = 0;
    }else if(Contains(name, "SingleLept")){ category = 1; bin = 0;
      if(Contains(name, "genMET-150")) bin = 1;
    }else if(Contains(name, "DiLept")){ category = 2; bin = 0;
      if(Contains(name, "genMET-150")) bin = 1;
    }else if(Contains(name, "TTJets_HT")){ category = 3;
      if(Contains(name, "HT-600to800")){ bin = 0;
      }else if(Contains(name, "HT-800to1200")){ bin = 1;
      }else if(Contains(name, "HT-1200to2500")){ bin = 2;
      }else if(Contains(name, "HT-2500toInf")){ bin = 3;
      }
    }else if(Contains(name, "TT_")){ category = 4; bin = 0;
    }else if(Contains(name, "TTJets_Mtt")){ category = 5; bin = 0;
    }
  }else if(Contains(name, "WJets") && !Contains(name, "TTWJets")){ sample = 2;
    if(Contains(name, "WJetsToLNu_Tune")){ category = 0; bin = 0;
    }else if(Contains(name, "WJetsToLNu_HT")){ category = 1;
      if(Contains(name, "HT-70To100")){ bin = 0;
      }else if(Contains(name, "HT-100To200")){ bin = 1;
      }else if(Contains(name, "HT-200To400")){ bin = 2;
      }else if(Contains(name, "HT-400To600")){ bin = 3;
      }else if(Contains(name, "HT-600To800")){ bin = 4;
      }else if(Contains(name, "HT-800To1200")){ bin = 5;
      }else if(Contains(name, "HT-1200To2500")){ bin = 6;
      }else if(Contains(name, "HT-2500ToInf")){ bin = 7;
      }else if(Contains(name, "HT-600ToInf")){ bin = 10;
      }
    }else if(Contains(name, "WJetsToQQ_HT")){ category = 2;
      if(Contains(name, "HT-600ToInf")){ bin = 0;
      }
    }
  }else if(Contains(name, "ST_")){ sample = 3;
    if(Contains(name, "ST_s-channel")){ category = 0; bin = 0;
    }else if(Contains(name, "ST_t-channel_top")){ category = 1; bin = 0;
    }else if(Contains(name, "ST_t-channel_antitop")){ category = 2; bin = 0;
    }else if(Contains(name, "ST_tW_top")){ category = 3; bin = 0;
    }else if(Contains(name, "ST_tW_antitop")){ category = 4; bin = 0;
    }
  }else if(Contains(name, "TTWJets")){ sample = 4;
    if(Contains(name, "TTWJetsToLNu")){ category = 0; bin = 0;
    }else if(Contains(name, "TTWJetsToQQ")){ category = 1; bin = 0;
    }
  }else if(Contains(name, "TTZ")){ sample = 5;
    if(Contains(name, "TTZToLLNuNu")){ category = 0; bin = 0;
    }else if(Contains(name, "TTZToQQ")){ category = 1; bin = 0;
    }
  }else if(Contains(name, "DYJetsToLL")){ sample = 6;
    if(Contains(name, "DYJetsToLL_M-50_Tune") && Contains(name, "madgraphMLM")){ category = 0; bin = 0;
    }else if(Contains(name, "DYJetsToLL_M-50_HT")){ category = 1;
      if(Contains(name, "HT-70to100")){ bin = 0;
      }else if(Contains(name, "HT-100to200")){ bin = 1;
      }else if(Contains(name, "HT-200to400")){ bin = 2;
      }else if(Contains(name, "HT-400to600")){ bin = 3;
      }else if(Contains(name, "HT-600to800")){ bin = 4;
      }else if(Contains(name, "HT-800to1200")){ bin = 5;
      }else if(Contains(name, "HT-1200to2500")){ bin = 6;
      }else if(Contains(name, "HT-2500toInf")){ bin = 7;
      }else if(Contains(name, "HT-600toInf")){ bin = 10;
      }
    } else if(Contains(name, "DYJetsToLL") && Contains(name, "amcatnloFXFX")){ category = 2; 
      if(Contains(name, "M-50_Tune")){ bin = 0;
      }else if(Contains(name, "Pt-50To100")){ bin = 1;
      }else if(Contains(name, "Pt-100To250")){ bin = 2;
      }else if(Contains(name, "Pt-250To400")){ bin = 3;
      }else if(Contains(name, "Pt-400To650")){ bin = 4;
      }else if(Contains(name, "Pt-650ToInf")){ bin = 5;
      }
    }

  }else if(Contains(name, "QCD")){ sample = 7;
    if(Contains(name, "QCD_HT")){ category = 0;
      if(Contains(name, "HT50to100")){ bin = 0;
      }else if(Contains(name, "HT100to200")){ bin = 1;
      }else if(Contains(name, "HT200to300")){ bin = 2;
      }else if(Contains(name, "HT300to500")){ bin = 3;
      }else if(Contains(name, "HT500to700")){ bin = 4;
      }else if(Contains(name, "HT700to1000")){ bin = 5;
      }else if(Contains(name, "HT1000to1500")){ bin = 6;
      }else if(Contains(name, "HT1500to2000")){ bin = 7;
      }else if(Contains(name, "HT2000toInf")){ bin = 8;
      }
    }else if(Contains(name, "QCD_Pt")){ category = 1;
      if(Contains(name, "5to10")){ bin = 0;
      }else if(Contains(name, "10to15")){ bin = 1;
      }else if(Contains(name, "15to30")){ bin = 2;
      }else if(Contains(name, "30to50")){ bin = 3;
      }else if(Contains(name, "50to80")){ bin = 4;
      }else if(Contains(name, "80to120")){ bin = 5;
      }else if(Contains(name, "120to170")){ bin = 6;
      }else if(Contains(name, "170to300")){ bin = 7;
      }else if(Contains(name, "300to470")){ bin = 8;
      }else if(Contains(name, "470to600")){ bin = 9;
      }else if(Contains(name, "600to800")){ bin = 10;
      }else if(Contains(name, "800to1000")){ bin = 11;
      }else if(Contains(name, "1000to1400")){ bin = 12;
      }else if(Contains(name, "1400to1800")){ bin = 13;
      }else if(Contains(name, "1800to2400")){ bin = 14;
      }else if(Contains(name, "2400to3200")){ bin = 15;
      }else if(Contains(name, "3200toInf")){ bin = 16;
      }
    }
  }else if(Contains(name, "ZJets")){ sample = 8;
    if(Contains(name, "ZJetsToNuNu")){ category = 0;
      if(Contains(name, "HT-70To100")){ bin = 0;
      }else if(Contains(name, "HT-100To200")){ bin = 1;
      }else if(Contains(name, "HT-200To400")){ bin = 2;
      }else if(Contains(name, "HT-400To600")){ bin = 3;
      }else if(Contains(name, "HT-600To800")){ bin = 4;
      }else if(Contains(name, "HT-800To1200")){ bin = 5;
      }else if(Contains(name, "HT-1200To2500")){ bin = 6;
      }else if(Contains(name, "HT-2500ToInf")){ bin = 7;
      }else if(Contains(name, "HT-600ToInf")){ bin = 10;
      }
    }else if(Contains(name, "ZJetsToQQ")){ category = 1;
      if(Contains(name, "HT600toInf")){ bin = 0;
      }
    }
  }else if(Contains(name, "ttH")){ sample = 9;
    if(Contains(name, "ttHJetTobb")){ category = 0; bin = 0;
    }
  }else if(Contains(name, "TTGJets")){ sample = 10;
    if(Contains(name, "TTGJets_Tune")){ category = 0; bin = 0;
    }
  }else if(Contains(name, "TTTT")){ sample = 11;
    if(Contains(name, "TTTT_Tune")){ category = 0; bin = 0;
    }
  }else if(Contains(name, "WH_") && !Contains(name,"TChiWH")){ sample = 12;
    if(Contains(name, "WH_HToBB_WToLNu")){ category = 0; bin = 0;
    }
  }else if(Contains(name, "ZH_")){ sample = 13;
    if(Contains(name, "ZH_HToBB_ZToNuNu")){ category = 0; bin = 0;
    }
  }else if(Contains(name, "WW") && !Contains(name,"TChiHH")){ sample = 14;
    if(Contains(name, "WWToLNuQQ")){ category = 0; bin = 0;
    }else if(Contains(name, "WWTo2L2Nu")){ category = 1; bin = 0;
    }else if(Contains(name, "WW_Tune"))  { category = 2; bin = 0;
    }
  }else if(Contains(name, "WZ") && !Contains(name,"TChiHH")){ sample = 15;
    if(Contains(name, "WZTo1L3Nu")){ category = 0; bin = 0;
    }else if(Contains(name, "WZTo1L1Nu2Q")){ category = 1; bin = 0;
    }else if(Contains(name, "WZTo2L2Q")){    category = 2; bin = 0;
    }else if(Contains(name, "WZTo3LNu")){    category = 3; bin = 0;
    }else if(Contains(name, "WZ_Tune")){     category = 4; bin = 0;
    }
  }else if(Contains(name, "ZZ") && !Contains(name,"TChiHH")){ sample = 16;
    if(Contains(name, "ZZ_Tune")){ category = 0; bin = 0;
    }
  }else if(Contains(name, "ZGTo")) { sample = 17; bin = 0;
    if(Contains(name,"ZGTo2LG_Tune")) category = 0;
    else if(Contains(name,"ZGToLLG_01J")) category = 1;
  }else if(Contains(name, "TGJets")) { sample = 18; bin = 0;
    if(Contains(name, "TGJets_Tune")) category = 0;
  }else if(Contains(name, "LLAJJ")) { sample = 19; bin = 0;
    if(Contains(name,"EWK_MLL-50")) category = 0;
  }else if(Contains(name, "T1tttt")){ sample = 100; category = 0; bin = 0;
  }else if(Contains(name, "T2tt")){   sample = 101; category = 0; bin = 0;
  }else if(Contains(name, "T1bbbb")){ sample = 102; category = 0; bin = 0;
  }else if(Contains(name, "T2bb")){   sample = 103; category = 0; bin = 0;
  }else if(Contains(name, "T1qqqq")){ sample = 104; category = 0; bin = 0;
  }else if(Contains(name, "RPV")){    sample = 105; category = 0; bin = 0;
  }else if(Contains(name, "TChiHH")){ sample = 106; category = 0; bin = 0;
  }else if(Contains(name, "HToZG")) { sample = 200; bin = 0;
    if(Contains(name,"GluGluH"))      category = 0; 
    else if(Contains(name,"VBF"))     category = 1; 
    else if(Contains(name,"WPlusH"))  category = 2; 
    else if(Contains(name,"WMinusH")) category = 3; 
    else if(Contains(name,"ZH"))      category = 4; 
    else if(Contains(name,"ttH"))     category = 5; 
  }

  if(sample < 0 || category < 0 || bin < 0
     || category > 9 || bin > 99){
    // DBG("Could not find type code for " << name << ": sample=" << sample << ", category=" << category << ", bin=" << bin);
    int code = -(1000*abs(sample)+100*abs(category)+abs(bin));
    if(code >= 0) code = -99999;
    return code;
  }else{
    int code = 1000*sample+100*category+bin;
    if(code < 0 || code > 107000){
      cout<<"ERROR:: Type code out of range for " << name << ": sample=" << sample << ", category=" << category << ", bin=" << bin;
    }
    return code;
  }
}

void EventTools::WriteTriggerEfficiency(pico_tree &pico) {
  // trigger efficiency and uncertainty - @todo, needs to be updated to full Run 2 trig eff. measurement
  pico.out_eff_trig() = hig_trig_eff::eff(pico);
  float effunc = hig_trig_eff::eff_unc(pico);
  pico.out_sys_trig().resize(2,0);
  pico.out_sys_trig()[0] = 1+effunc;
  pico.out_sys_trig()[1] = 1-effunc;

}
