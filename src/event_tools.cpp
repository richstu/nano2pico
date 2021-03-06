#include "event_tools.hpp"

#include "utilities.hpp"
#include "hig_trig_eff.hpp"
#include "TMath.h"

using namespace std;

EventTools::EventTools(const string &name_, int year_):
  name(name_),
  year(year_),
  isTTJets_LO_Incl(false),
  isTTJets_LO_MET(false),
  isTTJets_LO_HT(false),
  isWJets_LO(false),
  isDYJets_LO(false),
  dataset(-1){

  if(Contains(name, "TTJets_") && Contains(name, "genMET-") && Contains(name, "madgraphMLM")) 
    isTTJets_LO_MET = true;

  if(Contains(name, "TTJets_") && !Contains(name, "TTJets_HT") && !Contains(name, "genMET-") &&Contains(name, "madgraphMLM")) 
    isTTJets_LO_Incl = true;

  if(Contains(name, "TTJets_HT") && Contains(name, "madgraphMLM")) 
    isTTJets_LO_HT = true;

  if(Contains(name, "WJetsToLNu_Tune")  && Contains(name,"madgraphMLM"))
    isWJets_LO = true;
  
  if(Contains(name, "DYJetsToLL_M-50_Tune")  && Contains(name,"madgraphMLM"))
    isDYJets_LO = true;

  if(Contains(name, "EGamma")) // looks like this replaced SingleElectron and DoubleEG starting in 2018
    dataset = Dataset::EGamma;
  else if(Contains(name, "SingleElectron")) 
    dataset = Dataset::SingleElectron;
  else if(Contains(name, "SingleMuon")) 
    dataset = Dataset::SingleMuon;
  else if(Contains(name, "DoubleEG")) 
    dataset = Dataset::DoubleEG;
  else if(Contains(name, "DoubleMuon")) 
    dataset = Dataset::DoubleMuon;
  else if(Contains(name, "MET")) 
    dataset = Dataset::MET;
  else if(Contains(name, "JetHT")) 
    dataset = Dataset::JetHT;
}

EventTools::~EventTools(){
}

void EventTools::WriteStitch(nano_tree &nano, pico_tree &pico){
  pico.out_stitch_photon() = pico.out_stitch_htmet() = pico.out_stitch() = pico.out_stitch_ht() = true;
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

  if (isTTJets_LO_Incl || isTTJets_LO_MET || isTTJets_LO_HT) {
    //remove events covered by TTG, see AN-17-197 (TOP-18-010)
    //stitch if prompt photon w pt>13 |eta|<3 deltaR(genPart[pt>5])>0.2
    for (unsigned int mc_idx = 0; mc_idx < nano.GenPart_pdgId().size(); mc_idx++) {
      if (nano.GenPart_pdgId().at(mc_idx) == 22) {
	float ph_pt = nano.GenPart_pt().at(mc_idx);
	float ph_eta = nano.GenPart_eta().at(mc_idx);
	float ph_phi = nano.GenPart_phi().at(mc_idx);
        if (ph_pt > 13 && fabs(ph_eta)<3.0 && (nano.GenPart_statusFlags().at(mc_idx) & 0x1) == 1) {
	  //check if another genparticle nearby
	  bool deltar_fail = false;
          for (unsigned int mc_idx_2 = 0; mc_idx_2 < nano.GenPart_pdgId().size(); mc_idx_2++) {
            if (nano.GenPart_pt().at(mc_idx_2)>5 && dR(ph_eta,nano.GenPart_eta().at(mc_idx_2),ph_phi,nano.GenPart_phi().at(mc_idx_2))<0.2) {
              deltar_fail = true;
	      break;
	    }
          }
          if (!deltar_fail) {
            pico.out_stitch_photon() = pico.out_stitch() = false;
          }
	}
      }
    }
  }

  if(isDYJets_LO  && nano.LHE_HT()>70) 
    pico.out_stitch_htmet() = pico.out_stitch_ht() = pico.out_stitch() = false;
  
  if(isWJets_LO  && nano.LHE_HT()>70) 
    pico.out_stitch_htmet() = pico.out_stitch_ht() = pico.out_stitch() = false;
  return;
}


void EventTools::WriteDataQualityFilters(nano_tree& nano, pico_tree& pico, vector<int> sig_jet_nano_idx,
                                         float min_jet_pt, bool isData, bool isFastsim){
  float MET_pt, MET_phi;
  getMETWithJEC(nano, year, isFastsim, MET_pt, MET_phi);
  vector<float> Jet_pt, Jet_mass;
  getJetWithJEC(nano, isFastsim, Jet_pt, Jet_mass);

  // jet quality filter
  pico.out_pass_jets() = true;
  if (isFastsim) {
    // Fastsim: veto if certain central jets have no matching GenJet as per SUSY recommendation:
    // https://twiki.cern.ch/twiki/bin/view/CMS/SUSRecommendations18#Cleaning_up_of_fastsim_jets_from
    for(int ijet(0); ijet<nano.nJet(); ++ijet){
      if(Jet_pt[ijet] > 20 && fabs(nano.Jet_eta()[ijet])<=2.5 && nano.Jet_chHEF()[ijet] < 0.1) {
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
      if (Jet_pt[ijet] > min_jet_pt && nano.Jet_jetId()[ijet] < 1) 
        pico.out_pass_jets() = false;
    } 
  }

  // RA2b filters  
  pico.out_pass_muon_jet() = true; 
  for (auto &idx: sig_jet_nano_idx){
    // if (abs(nano.Jet_eta()[idx])>2.4) continue; -> already enforced in signal jet selection
    // if is overlapping with lepton -> already enforced in signal jet selection
    if (Jet_pt[idx]<=200.) continue;
    if (nano.Jet_muEF()[idx]<=0.5) continue;
    if (DeltaPhi(nano.Jet_phi()[idx],MET_phi)<(TMath::Pi()-0.4)) continue;
    pico.out_pass_muon_jet() = false;
    break;
  }

  pico.out_pass_low_neutral_jet() = true;
  for(int ijet(0); ijet<nano.nJet(); ++ijet){  
    if (nano.Jet_neEmEF()[ijet] <0.03 && DeltaPhi(nano.Jet_phi()[ijet], pico.out_met_phi())>(TMath::Pi()-0.4))
      pico.out_pass_low_neutral_jet() = false;
    break; //only apply to leading jet
  }

  pico.out_pass_htratio_dphi_tight() = true;
  float htratio = pico.out_ht5()/pico.out_ht();
  for(int ijet(0); ijet<nano.nJet(); ++ijet){  
    if (htratio >= 1.2 && DeltaPhi(nano.Jet_phi()[ijet], pico.out_met_phi()) < (5.3*htratio - 4.78)) 
      pico.out_pass_htratio_dphi_tight() = false;
    break; //only apply to leading jet
  }

  pico.out_pass_ecalnoisejet() = true;
  if (year!=2016) {
    int counter = 0;
    bool goodjet[2] = {true, true};
    double dphi = 0.;
    for (int ijet(0); ijet < nano.nJet(); ijet++) {
      if (counter >= 2) break;
      float jet_pt = nano.Jet_pt()[ijet];
      if (isFastsim) jet_pt = nano.Jet_pt_nom()[ijet];
      if (jet_pt>30 && fabs(nano.Jet_eta()[ijet])>2.4 && fabs(nano.Jet_eta()[ijet])<5.0) {
        dphi = DeltaPhi(nano.Jet_phi()[ijet], pico.out_met_phi());
        if (nano.Jet_pt()[ijet]>250 && (dphi > 2.6 || dphi < 0.1)) goodjet[counter] = false;
        ++counter;
      }
    }
    pico.out_pass_ecalnoisejet() = goodjet[0] && goodjet[1];
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
  } else {
    pico.out_pass_badcalib() = nano.Flag_ecalBadCalibFilterV2();
  }
  pico.out_pass_badchhad() = nano.Flag_BadChargedCandidateFilter();
  pico.out_pass_badpfmu() = nano.Flag_BadPFMuonFilter();
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

bool EventTools::SaveTriggerDecisions(nano_tree& nano, pico_tree& pico, bool isZgamma){

  bool egamma_trigs = nano.HLT_Ele25_WPTight_Gsf() || nano.HLT_Ele27_WPTight_Gsf() || 
                      nano.HLT_Ele28_WPTight_Gsf() || nano.HLT_Ele32_WPTight_Gsf() ||
	              nano.HLT_Ele32_WPTight_Gsf_L1DoubleEG() || nano.HLT_Ele35_WPTight_Gsf() || 
                      nano.HLT_Ele20_WPLoose_Gsf() || nano.HLT_Ele45_WPLoose_Gsf() ||
                      nano.HLT_Ele105_CaloIdVT_GsfTrkIdT() || nano.HLT_Ele115_CaloIdVT_GsfTrkIdT() ||
                      nano.HLT_Ele135_CaloIdVT_GsfTrkIdT() || nano.HLT_Ele145_CaloIdVT_GsfTrkIdT() ||
                      nano.HLT_Ele25_eta2p1_WPTight_Gsf() || nano.HLT_Ele27_eta2p1_WPTight_Gsf() || 
                      nano.HLT_Ele20_eta2p1_WPLoose_Gsf() || nano.HLT_Ele25_eta2p1_WPLoose_Gsf() ||
                      nano.HLT_Ele27_eta2p1_WPLoose_Gsf() || nano.HLT_Ele15_IsoVVVL_PFHT350() ||
		      nano.HLT_Ele15_IsoVVVL_PFHT400() || nano.HLT_Ele15_IsoVVVL_PFHT450() ||
		      nano.HLT_Ele15_IsoVVVL_PFHT600() || nano.HLT_Ele50_IsoVVVL_PFHT450();

  pico.out_HLT_Ele25_WPTight_Gsf() = nano.HLT_Ele25_WPTight_Gsf();
  pico.out_HLT_Ele27_WPTight_Gsf() = nano.HLT_Ele27_WPTight_Gsf();
  pico.out_HLT_Ele28_WPTight_Gsf() = nano.HLT_Ele28_WPTight_Gsf();
  pico.out_HLT_Ele32_WPTight_Gsf() = nano.HLT_Ele32_WPTight_Gsf();
  pico.out_HLT_Ele32_WPTight_Gsf_L1DoubleEG() = nano.HLT_Ele32_WPTight_Gsf_L1DoubleEG();
  pico.out_HLT_Ele35_WPTight_Gsf() = nano.HLT_Ele35_WPTight_Gsf();
  pico.out_HLT_Ele20_WPLoose_Gsf() = nano.HLT_Ele20_WPLoose_Gsf();
  pico.out_HLT_Ele45_WPLoose_Gsf() = nano.HLT_Ele45_WPLoose_Gsf();
  pico.out_HLT_Ele105_CaloIdVT_GsfTrkIdT() = nano.HLT_Ele105_CaloIdVT_GsfTrkIdT();
  pico.out_HLT_Ele115_CaloIdVT_GsfTrkIdT() = nano.HLT_Ele115_CaloIdVT_GsfTrkIdT();
  pico.out_HLT_Ele135_CaloIdVT_GsfTrkIdT() = nano.HLT_Ele135_CaloIdVT_GsfTrkIdT();
  pico.out_HLT_Ele145_CaloIdVT_GsfTrkIdT() = nano.HLT_Ele145_CaloIdVT_GsfTrkIdT();
  pico.out_HLT_Ele25_eta2p1_WPTight_Gsf() = nano.HLT_Ele25_eta2p1_WPTight_Gsf();
  pico.out_HLT_Ele27_eta2p1_WPTight_Gsf() = nano.HLT_Ele27_eta2p1_WPTight_Gsf();
  pico.out_HLT_Ele20_eta2p1_WPLoose_Gsf() = nano.HLT_Ele20_eta2p1_WPLoose_Gsf();
  pico.out_HLT_Ele25_eta2p1_WPLoose_Gsf() = nano.HLT_Ele25_eta2p1_WPLoose_Gsf();
  pico.out_HLT_Ele27_eta2p1_WPLoose_Gsf() = nano.HLT_Ele27_eta2p1_WPLoose_Gsf();
  pico.out_HLT_Ele15_IsoVVVL_PFHT350() = nano.HLT_Ele15_IsoVVVL_PFHT350();
  pico.out_HLT_Ele15_IsoVVVL_PFHT400() = nano.HLT_Ele15_IsoVVVL_PFHT400();
  pico.out_HLT_Ele15_IsoVVVL_PFHT450() = nano.HLT_Ele15_IsoVVVL_PFHT450();
  pico.out_HLT_Ele15_IsoVVVL_PFHT600() = nano.HLT_Ele15_IsoVVVL_PFHT600();
  pico.out_HLT_Ele50_IsoVVVL_PFHT450() = nano.HLT_Ele50_IsoVVVL_PFHT450();

  bool muon_trigs = nano.HLT_IsoMu20() || nano.HLT_IsoMu22() || nano.HLT_IsoMu24() ||
                    nano.HLT_IsoMu27() || nano.HLT_IsoTkMu20() || nano.HLT_IsoTkMu22() ||
                    nano.HLT_IsoTkMu24() || nano.HLT_Mu50() || nano.HLT_Mu55() ||
                    nano.HLT_TkMu50() || nano.HLT_IsoMu22_eta2p1() || nano.HLT_IsoMu24_eta2p1() ||
                    nano.HLT_Mu45_eta2p1() || nano.HLT_Mu15_IsoVVVL_PFHT350() || nano.HLT_Mu15_IsoVVVL_PFHT400() ||
                    nano.HLT_Mu15_IsoVVVL_PFHT450() || nano.HLT_Mu15_IsoVVVL_PFHT600() || nano.HLT_Mu50_IsoVVVL_PFHT400() ||
                    nano.HLT_Mu50_IsoVVVL_PFHT450();

  pico.out_HLT_IsoMu20() = nano.HLT_IsoMu20();
  pico.out_HLT_IsoMu22() = nano.HLT_IsoMu22();
  pico.out_HLT_IsoMu24() = nano.HLT_IsoMu24();
  pico.out_HLT_IsoMu27() = nano.HLT_IsoMu27();
  pico.out_HLT_IsoTkMu20() = nano.HLT_IsoTkMu20();
  pico.out_HLT_IsoTkMu22() = nano.HLT_IsoTkMu22();
  pico.out_HLT_IsoTkMu24() = nano.HLT_IsoTkMu24();
  pico.out_HLT_Mu50() = nano.HLT_Mu50();
  pico.out_HLT_Mu55() = nano.HLT_Mu55();
  pico.out_HLT_TkMu50() = nano.HLT_TkMu50();
  pico.out_HLT_IsoMu22_eta2p1() = nano.HLT_IsoMu22_eta2p1();
  pico.out_HLT_IsoMu24_eta2p1() = nano.HLT_IsoMu24_eta2p1();
  pico.out_HLT_Mu45_eta2p1() = nano.HLT_Mu45_eta2p1();
  pico.out_HLT_Mu15_IsoVVVL_PFHT350() = nano.HLT_Mu15_IsoVVVL_PFHT350();
  pico.out_HLT_Mu15_IsoVVVL_PFHT400() = nano.HLT_Mu15_IsoVVVL_PFHT400();
  pico.out_HLT_Mu15_IsoVVVL_PFHT450() = nano.HLT_Mu15_IsoVVVL_PFHT450();
  pico.out_HLT_Mu15_IsoVVVL_PFHT600() = nano.HLT_Mu15_IsoVVVL_PFHT600();
  pico.out_HLT_Mu50_IsoVVVL_PFHT400() = nano.HLT_Mu50_IsoVVVL_PFHT400();
  pico.out_HLT_Mu50_IsoVVVL_PFHT450() = nano.HLT_Mu50_IsoVVVL_PFHT450();

  bool met_trigs = nano.HLT_PFMET90_PFMHT90_IDTight() || nano.HLT_PFMETNoMu90_PFMHTNoMu90_IDTight() || 
                   nano.HLT_PFMET100_PFMHT100_IDTight() || nano.HLT_PFMETNoMu100_PFMHTNoMu100_IDTight() || 
		   nano.HLT_PFMET110_PFMHT110_IDTight() || nano.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight() || 
                   nano.HLT_PFMET120_PFMHT120_IDTight() || nano.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight() || 
                   nano.HLT_PFMET130_PFMHT130_IDTight() || nano.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight() || 
                   nano.HLT_PFMET140_PFMHT140_IDTight() || nano.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight() || 
                   nano.HLT_PFMET100_PFMHT100_IDTight_PFHT60() || nano.HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60() ||
                   nano.HLT_PFMET110_PFMHT110_IDTight_PFHT60() || nano.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60() ||
                   nano.HLT_PFMET120_PFMHT120_IDTight_PFHT60() || nano.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60() ||
                   nano.HLT_PFMET130_PFMHT130_IDTight_PFHT60() || nano.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60() ||
                   nano.HLT_PFMET140_PFMHT140_IDTight_PFHT60() || nano.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_PFHT60() ||
                   nano.HLT_PFMET120_PFMHT120_IDTight_HFCleaned() || nano.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned() || nano.HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned() ||
                   nano.HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1() || nano.HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1() || nano.HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1();
  //MET_MHT[NoMu]
  pico.out_HLT_PFMET90_PFMHT90_IDTight() = nano.HLT_PFMET90_PFMHT90_IDTight();
  pico.out_HLT_PFMETNoMu90_PFMHTNoMu90_IDTight() = nano.HLT_PFMETNoMu90_PFMHTNoMu90_IDTight();
  pico.out_HLT_PFMET100_PFMHT100_IDTight() = nano.HLT_PFMET100_PFMHT100_IDTight();
  pico.out_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight() = nano.HLT_PFMETNoMu100_PFMHTNoMu100_IDTight();
  pico.out_HLT_PFMET110_PFMHT110_IDTight() = nano.HLT_PFMET110_PFMHT110_IDTight();
  pico.out_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight() = nano.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight();
  pico.out_HLT_PFMET120_PFMHT120_IDTight() = nano.HLT_PFMET120_PFMHT120_IDTight();
  pico.out_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight() = nano.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight();
  pico.out_HLT_PFMET130_PFMHT130_IDTight() = nano.HLT_PFMET130_PFMHT130_IDTight();
  pico.out_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight() = nano.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight();
  pico.out_HLT_PFMET140_PFMHT140_IDTight() = nano.HLT_PFMET140_PFMHT140_IDTight();
  pico.out_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight() = nano.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight();
  //MET_MHT[NoMu]_HT60
  pico.out_HLT_PFMET100_PFMHT100_IDTight_PFHT60() = nano.HLT_PFMET100_PFMHT100_IDTight_PFHT60();
  pico.out_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60() = nano.HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60();
  pico.out_HLT_PFMET110_PFMHT110_IDTight_PFHT60() = nano.HLT_PFMET110_PFMHT110_IDTight_PFHT60();
  pico.out_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60() = nano.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60();
  pico.out_HLT_PFMET120_PFMHT120_IDTight_PFHT60() = nano.HLT_PFMET120_PFMHT120_IDTight_PFHT60();
  pico.out_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60() = nano.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60();
  pico.out_HLT_PFMET130_PFMHT130_IDTight_PFHT60() = nano.HLT_PFMET130_PFMHT130_IDTight_PFHT60();
  pico.out_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60() = nano.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60();
  pico.out_HLT_PFMET140_PFMHT140_IDTight_PFHT60() = nano.HLT_PFMET140_PFMHT140_IDTight_PFHT60();
  pico.out_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_PFHT60() = nano.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_PFHT60();
  //HF cleaned
  pico.out_HLT_PFMET120_PFMHT120_IDTight_HFCleaned() = nano.HLT_PFMET120_PFMHT120_IDTight_HFCleaned();
  pico.out_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned() = nano.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned();
  pico.out_HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned() = nano.HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned();
  //experimental b-tag cross triggers
  pico.out_HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1() = nano.HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1();
  pico.out_HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1() = nano.HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1();
  pico.out_HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1() = nano.HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1();

  // Jet/HT trigger
  bool jetht_trigs = nano.HLT_PFJet500() || nano.HLT_PFHT125() || nano.HLT_PFHT200() || nano.HLT_PFHT300() || 
                                            nano.HLT_PFHT400() || nano.HLT_PFHT475() || nano.HLT_PFHT600() || 
					    nano.HLT_PFHT650() || nano.HLT_PFHT800() || nano.HLT_PFHT900() || 
					    nano.HLT_PFHT180() || nano.HLT_PFHT370() || nano.HLT_PFHT430() || 
					    nano.HLT_PFHT510() || nano.HLT_PFHT590() || nano.HLT_PFHT680() || 
					    nano.HLT_PFHT780() || nano.HLT_PFHT890() || nano.HLT_PFHT1050() || 
					    nano.HLT_PFHT250() || nano.HLT_PFHT350();
  pico.out_HLT_PFJet500() = nano.HLT_PFJet500();
  pico.out_HLT_PFHT125() = nano.HLT_PFHT125();
  pico.out_HLT_PFHT200() = nano.HLT_PFHT125();
  pico.out_HLT_PFHT300() = nano.HLT_PFHT300();
  pico.out_HLT_PFHT400() = nano.HLT_PFHT400();
  pico.out_HLT_PFHT475() = nano.HLT_PFHT475();
  pico.out_HLT_PFHT600() = nano.HLT_PFHT600();
  pico.out_HLT_PFHT650() = nano.HLT_PFHT650();
  pico.out_HLT_PFHT800() = nano.HLT_PFHT800();
  pico.out_HLT_PFHT900() = nano.HLT_PFHT900();
  pico.out_HLT_PFHT180() = nano.HLT_PFHT180();
  pico.out_HLT_PFHT370() = nano.HLT_PFHT370();
  pico.out_HLT_PFHT430() = nano.HLT_PFHT430();
  pico.out_HLT_PFHT510() = nano.HLT_PFHT510();
  pico.out_HLT_PFHT590() = nano.HLT_PFHT590();
  pico.out_HLT_PFHT680() = nano.HLT_PFHT680();
  pico.out_HLT_PFHT780() = nano.HLT_PFHT780();
  pico.out_HLT_PFHT890() = nano.HLT_PFHT890();
  pico.out_HLT_PFHT1050() = nano.HLT_PFHT1050();
  pico.out_HLT_PFHT250() = nano.HLT_PFHT250();
  pico.out_HLT_PFHT350() = nano.HLT_PFHT350();

  // ZGamma triggers
  pico.out_HLT_Mu17_Photon30_IsoCaloId()               = nano.HLT_Mu17_Photon30_IsoCaloId();
  pico.out_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() = nano.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ();
  pico.out_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL()    = nano.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL();
  pico.out_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL()          = nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL();
  pico.out_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL()        = nano.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL();
  pico.out_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ()       = nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ();
  pico.out_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ()     = nano.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ();
  pico.out_HLT_Photon175()                             = nano.HLT_Photon175();


  if (isZgamma)
    return true;
  else {
    // this assumes that we process either all the datasets or at least an ordered subset starting with the MET 
    // e.g. after measuring trigger efficiency, it would no longer be necessary to run on JetHT
    if (dataset==Dataset::MET                                                                && met_trigs) return true;
    else if (year==2018 && dataset==Dataset::EGamma                         && egamma_trigs && !met_trigs) return true;
    else if ((year==2016 || year==2017) && dataset==Dataset::SingleElectron && egamma_trigs && !met_trigs) return true;
    else if (dataset==Dataset::SingleMuon                    && muon_trigs && !egamma_trigs && !met_trigs) return true;
    else if (dataset==Dataset::JetHT         && jetht_trigs && !muon_trigs && !egamma_trigs && !met_trigs) return true;
    else return false;
  }

  return false;
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
  }else if(Contains(name, "ttHJet") && !Contains(name,"HToZG")){ sample = 9;
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
  }else if(Contains(name, "ZH") && !Contains(name,"HToZG") && !Contains(name, "T5qqqqZH")){ sample = 13;
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
  }else if(Contains(name, "T5qqqqZH")){ sample = 107; category = 0; bin = 0;
  }else if(Contains(name, "HToZG")) { sample = 200; bin = 0;
    if(Contains(name,"GluGluH"))      category = 0; 
    else if(Contains(name,"VBF"))     category = 1; 
    else if(Contains(name,"WplusH"))  category = 2; 
    else if(Contains(name,"WminusH")) category = 3; 
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
