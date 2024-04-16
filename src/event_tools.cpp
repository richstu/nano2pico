#include "event_tools.hpp"
#include "utilities.hpp"
#include "hig_trig_eff.hpp"
#include "TMath.h"
#include <bitset>

using namespace std;

EventTools::EventTools(const string &name_, int year_, bool isData_, float nanoaod_version_):
  name(name_),
  year(year_),
  isTTJets_LO_Incl(false),
  isTTJets_LO_MET(false),
  isTTJets_LO_HT(false),
  isWJets_LO(false),
  isDYJets_LO(false),
  isEWKZ(false),
  isWJ(false),
  isWW(false),
  isWZ(false),
  isZZ(false),
  isFastSim(false),
  isData(isData_),
  nanoaod_version(nanoaod_version_),
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

  if(Contains(name, "Fast"))
    isFastSim = true;


  //These four variables control the generator settings of the overlap removal variable in MC
  if(Contains(name, "WJetsToLNu") || Contains(name,"WGToLNuG_01J"))
    isWJ = true;

  if(Contains(name,"EWKZ2Jets") || Contains(name, "ZGamma2JToGamma2L2J_EWK"))
    isEWKZ = true;

  if(Contains(name, "WW") && !Contains(name,"WWW") && !Contains(name,"WWZ") && !Contains(name,"HToWW") && !Contains(name,"TChiHH"))
    isWW = true; //WW or WWG

  if(Contains(name, "WZ") && !Contains(name,"WWZ") && !Contains(name,"WZZ") && !Contains(name,"TChiHH"))
    isWZ = true; //WZ or WZG

  if(Contains(name, "ZZ") && !Contains(name,"WZZ") && !Contains(name,"ZZZ") && !Contains(name,"HToZZ") && !Contains(name,"TChiHH"))
    isZZ = true; //ZZ or ZZG

  if(Contains(name, "EGamma")) // replaced SingleElectron and DoubleEG starting in 2018
    dataset = Dataset::EGamma;
  else if(Contains(name, "SingleElectron")) 
    dataset = Dataset::SingleElectron;
  else if(Contains(name, "SingleMuon")) 
    dataset = Dataset::SingleMuon;
  else if(Contains(name, "DoubleEG")) 
    dataset = Dataset::DoubleEG;
  else if(Contains(name, "DoubleMuon")) 
    dataset = Dataset::DoubleMuon;
  else if(Contains(name, "MuonEG"))
    dataset = Dataset::MuonEG;
  else if(Contains(name, "Muon") && !Contains(name,"DoubleMuon") && !Contains(name,"SingleMuon") && !Contains(name,"MuonEG"))  //replaced SingleMuon and DoubleMuon starting in 2022
    dataset = Dataset::Muon;
  else if(Contains(name, "MET")) 
    dataset = Dataset::MET;
  else if(Contains(name, "JetHT")) 
    dataset = Dataset::JetHT;
  else if(Contains(name, "JetMET")) //replaced JetHT and MET starting in 2022
    dataset = Dataset::JetMET;
}

EventTools::~EventTools(){
}

void EventTools::WriteStitch(nano_tree &nano, pico_tree &pico){
  pico.out_stitch_photon() = true;
  pico.out_stitch_htmet() = true;
  pico.out_stitch() = true;
  pico.out_stitch_ht() = true;
  pico.out_is_overlap_old() = true;
  pico.out_is_overlap() = false; 
  pico.out_old_stitch_dy() = true;

  if (isData){
    pico.out_is_overlap() = false; pico.out_use_event() = true; 
    return;
  }

  if(isTTJets_LO_Incl && !isFastSim) {
    if (nano.LHE_HTIncoming()>600.f) 
      pico.out_stitch_htmet() = pico.out_stitch_ht() = false;
    if(year==2018 && nano.GenMET_pt()>80.f)
      pico.out_stitch_htmet() = pico.out_stitch() = false;
    else if (nano.GenMET_pt()>150.f) 
      pico.out_stitch_htmet() = pico.out_stitch() = false;
  }

  if (isTTJets_LO_MET && nano.LHE_HTIncoming()>600.f && !isFastSim) 
      pico.out_stitch_htmet() = pico.out_stitch_ht() = false;

  if ((isTTJets_LO_Incl || isTTJets_LO_MET || isTTJets_LO_HT) && !isFastSim) {
    //remove events covered by TTG, see AN-17-197 (TOP-18-010)
    //stitch if prompt photon w pt>13 |eta|<3 deltaR(genPart[pt>5])>0.2
    for (unsigned int mc_idx = 0; mc_idx < nano.GenPart_pdgId().size(); mc_idx++) {
      if (nano.GenPart_pdgId().at(mc_idx) == 22) {
        float ph_pt = nano.GenPart_pt().at(mc_idx);
        float ph_eta = nano.GenPart_eta().at(mc_idx);
        float ph_phi = nano.GenPart_phi().at(mc_idx);
        if (ph_pt > 13.f && fabs(ph_eta)<3.0.f && (nano.GenPart_statusFlags().at(mc_idx) & 0x1) == 1) {
          //check if another genparticle nearby
          bool deltar_fail = false;
          for (unsigned int mc_idx_2 = 0; mc_idx_2 < nano.GenPart_pdgId().size(); mc_idx_2++) {
            if (nano.GenPart_pt().at(mc_idx_2)>5.f && dR(ph_eta,nano.GenPart_eta().at(mc_idx_2),ph_phi,nano.GenPart_phi().at(mc_idx_2))<0.2f) {
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

  //Need to include the new overlap removal by including the newer values for ZGtoLLG_lowMll_lowGPt
  double ptmin = 9.0;
  double isocone = 0.05;
  if(isZZ || isTTJets_LO_Incl || Contains(name,"TTGJets") || isEWKZ){
    ptmin = 10.0;
  }
  if(isWZ || isWW){
    ptmin = 20.0;
  }
  if(isWJ){
    ptmin = 15.0;
  }

  double ptmin_old = 15.0;
  double etamax_old = 2.6;
  double isocone_old = 0.05;
  if(isWW || isZZ || isTTJets_LO_Incl || Contains(name,"TTGJets")){
    ptmin_old = 10.0;
    etamax_old= 99.0;
  }

  if(isWW || isWZ){
    ptmin_old = 20;
  }

  bool found_hadronic_w = false;
  bool found_higgs = false;
  int ntrulep = 0; //includes taus, unlike pico branch

  for(int mc_idx(0); mc_idx<nano.nGenPart(); mc_idx++) {

    bitset<15> mc_statusFlags(nano.GenPart_statusFlags().at(mc_idx));

    //check if event has a hadronic w decay
    int mom_id = -1;
    int mc_id = nano.GenPart_pdgId().at(mc_idx);
    if (nano.GenPart_genPartIdxMother().at(mc_idx) != -1)
      mom_id = nano.GenPart_pdgId().at(nano.GenPart_genPartIdxMother().at(mc_idx));
    if (abs(mom_id)==24 && abs(mc_id)>=1 && abs(mc_id)<=4)
      found_hadronic_w = true;

    //check if event has Higgs
    if (mc_id==25)
      found_higgs = true;

    //calcaulate number of leptons in event - require is FirstCopy and not isDirectPromptTauDecayProduct and isPrompt
    if ((abs(mc_id)==11 || abs(mc_id)==13 || abs(mc_id)==15) && !mc_statusFlags[5] && mc_statusFlags[12] && mc_statusFlags[0])
      ntrulep++;

    if( nano.GenPart_pdgId().at(mc_idx) == 22 ){ // photons 
      TVector3 compPart,genPhoton;

      if( (mc_statusFlags[0] || mc_statusFlags[8]) ){  // Which are isPrompt or fromHardProcess
        genPhoton.SetPtEtaPhi(nano.GenPart_pt().at(mc_idx), 
                            nano.GenPart_eta().at(mc_idx), 
                            nano.GenPart_phi().at(mc_idx));

        if( genPhoton.Pt() > ptmin){
          //check if another generator particle nearby
          bool found_other_particles = false;
          for (int mc_idx_2 = 0; mc_idx_2 < nano.nGenPart(); mc_idx_2++) {
            bitset<15> mc_statusFlags2(nano.GenPart_statusFlags().at(mc_idx_2));
            compPart.SetPtEtaPhi(nano.GenPart_pt().at(mc_idx_2), 
                                 nano.GenPart_eta().at(mc_idx_2), 
                                 nano.GenPart_phi().at(mc_idx_2)); 
            

            //isPrompt and fromHardProcess already applied
            if ( (compPart.Pt() > 5.0f) && (genPhoton.DeltaR(compPart) < isocone) && (mc_idx != mc_idx_2) && (mc_statusFlags2[8]) && (nano.GenPart_pdgId().at(mc_idx_2) != 22 )  ) {
              found_other_particles = true;
              //Basically saying that a photon is not isolated so this is not SM Zgamma sample!
            }
          }

          if(!found_other_particles){
            pico.out_is_overlap() = true; 
          }       

        }
      }


      if( (mc_statusFlags[0] || mc_statusFlags[8]) ){  // Which are isPrompt or fromHardProcess
        genPhoton.SetPtEtaPhi(nano.GenPart_pt().at(mc_idx), 
                            nano.GenPart_eta().at(mc_idx), 
                            nano.GenPart_phi().at(mc_idx));

        if( genPhoton.Pt() > ptmin_old && fabs(genPhoton.Eta())< etamax_old ){
          //check if another generator particle nearby
          bool found_other_particles = false;
          for (int mc_idx_2 = 0; mc_idx_2 < nano.nGenPart(); mc_idx_2++) {
            bitset<15> mc_statusFlags2(nano.GenPart_statusFlags().at(mc_idx_2));
            compPart.SetPtEtaPhi(nano.GenPart_pt().at(mc_idx_2), 
                                 nano.GenPart_eta().at(mc_idx_2), 
                                 nano.GenPart_phi().at(mc_idx_2)); 
            

            //isPrompt and fromHardProcess already applied
            if ( (compPart.Pt() > 5.0f) && (genPhoton.DeltaR(compPart) < isocone_old) && (mc_idx != mc_idx_2) && (mc_statusFlags2[8]) && (nano.GenPart_pdgId().at(mc_idx_2) != 22 )  ) {
              found_other_particles = true;
              //Basically saying that a photon is not isolated so this is not SM Zgamma sample!
            }
          }

          if(!found_other_particles){
            pico.out_is_overlap_old() = false;
          }       

        }
      }

      if( (mc_statusFlags[0] || mc_statusFlags[8]) && nano.GenPart_status().at(mc_idx) == 1 ){  // Which are isPrompt or fromHardProcess and stable
        for(size_t reco_idx(0); reco_idx < pico.out_photon_pt().size(); reco_idx++){
          if(pico.out_photon_sig()[reco_idx]){
             compPart.SetPtEtaPhi(pico.out_photon_pt().at(reco_idx), 
                                  pico.out_photon_eta().at(reco_idx), 
                                  pico.out_photon_phi().at(reco_idx));
            if(genPhoton.DeltaR(compPart) < 0.1f){
              pico.out_old_stitch_dy() = false;
            }
          }
        }
      }

    } //GenPart_pdgId==22
  } //loop over GenParts

  //This bit of code uses the overlap removal variable to then select whether an event should be kept or not. 
  //If the event contains should and does (does not) contain a generator photon then is_overlap_old = false, is_overlap = false, use_event = true (is_overlap_old=true, is_overlap=true, use_event=false)
  //If the event contains should not and does not (does) contain a generator photon then is_overlap_old = true, is_overlap = false, use_event = true (is_overlap_old=false, is_overlap=true, use_event=false)
  if( Contains(name,"ZZG") && found_higgs ){ //remove Higgs decays from ZZG sample
    pico.out_is_overlap() = true; pico.out_use_event() = false;
  } else if( (Contains(name, "TTGJets") ) || (isWZ && Contains(name, "WZG")) 
             || (isZZ && Contains(name, "ZZG")) || Contains(name,"ZGToLLG") 
             || (isWJ && Contains(name,"WGToLNuG_01J")) 
             || (isEWKZ && Contains(name,"ZGamma2JToGamma2L2J_EWK")) ){
    pico.out_use_event() = pico.out_is_overlap();
  } else if( isTTJets_LO_Incl || Contains(name, "TTTo2L2Nu") || isWJ || isEWKZ
             || (Contains(name, "DY") && !Contains(name, "DYG")) ){ //no WW due to bug in WWG
    pico.out_use_event() = !pico.out_is_overlap();
  } else if( isWZ && !found_hadronic_w ){ //WZG only generated with W->lnu
    pico.out_use_event() = !pico.out_is_overlap();
  } else if( isZZ && ntrulep==4 ){ //ZZ only generated with ZZ->4l
    pico.out_use_event() = !pico.out_is_overlap();
  } else {
    pico.out_is_overlap() = false; pico.out_use_event() = true;
  }
  //Note: removed overlap removal for WW and WWG since WWG sample may have bug

  if(isDYJets_LO  && nano.LHE_HT()>70f) 
    pico.out_stitch_htmet() = pico.out_stitch_ht() = pico.out_stitch() = false;
  
  if(isWJets_LO  && nano.LHE_HT()>70f) 
    pico.out_stitch_htmet() = pico.out_stitch_ht() = pico.out_stitch() = false;
  return;
}


void EventTools::WriteDataQualityFilters(nano_tree& nano, pico_tree& pico, vector<int> sig_jet_nano_idx,
                                         float min_jet_pt, bool isFastsim, bool is_preUL){
  float MET_pt, MET_phi;
  getMETWithJEC(nano, year, isFastsim, MET_pt, MET_phi, is_preUL);
  vector<float> Jet_pt, Jet_mass;
  getJetWithJEC(nano, isFastsim, Jet_pt, Jet_mass);
  vector<int> Jet_jetId;
  getJetId(nano, nanoaod_version, Jet_jetId);

  // jet quality filter
  pico.out_pass_jets() = true;
  if (isFastsim) {
    // Fastsim: veto if certain central jets have no matching GenJet as per SUSY recommendation:
    // https://twiki.cern.ch/twiki/bin/view/CMS/SUSRecommendations18#Cleaning_up_of_fastsim_jets_from
    for(int ijet(0); ijet<nano.nJet(); ++ijet){
      if(Jet_pt[ijet] > 20f && fabs(nano.Jet_eta()[ijet])<=2.5f && nano.Jet_chHEF()[ijet] < 0.1f) {
        bool found_match = false;
        for(int igenjet(0); igenjet<nano.nGenJet(); ++igenjet){
          if (dR(nano.Jet_eta()[ijet], nano.GenJet_eta()[igenjet], nano.Jet_phi()[ijet], nano.GenJet_phi()[igenjet])<=0.3f) {
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
      if (Jet_pt[ijet] > min_jet_pt && Jet_jetId[ijet] < 1) 
        pico.out_pass_jets() = false;
    } 
  }

  // RA2b filters  
  pico.out_pass_muon_jet() = true; 
  for (auto &idx: sig_jet_nano_idx){
    // if (abs(nano.Jet_eta()[idx])>2.4) continue; -> already enforced in signal jet selection
    // if is overlapping with lepton -> already enforced in signal jet selection
    if (Jet_pt[idx]<=200.f) continue;
    if (nano.Jet_muEF()[idx]<=0.5f) continue;
    if (DeltaPhi(nano.Jet_phi()[idx],MET_phi)<(TMath::Pi()-0.4f)) continue;
    pico.out_pass_muon_jet() = false;
    break;
  }

  pico.out_pass_low_neutral_jet() = true;
  for(int ijet(0); ijet<nano.nJet();){  
    if (nano.Jet_neEmEF()[ijet] <0.03f && DeltaPhi(nano.Jet_phi()[ijet], pico.out_met_phi())>(TMath::Pi()-0.4f))
      pico.out_pass_low_neutral_jet() = false;
    break; //only apply to leading jet
  }

  pico.out_pass_htratio_dphi_tight() = true;
  float htratio = pico.out_ht5()/pico.out_ht();
  for(int ijet(0); ijet<nano.nJet();){  
    if (htratio >= 1.2f && DeltaPhi(nano.Jet_phi()[ijet], pico.out_met_phi()) < (5.3f*htratio - 4.78f)) 
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
      if (jet_pt>30 && fabs(nano.Jet_eta()[ijet])>2.4f && fabs(nano.Jet_eta()[ijet])<5.0f) {
        dphi = DeltaPhi(nano.Jet_phi()[ijet], pico.out_met_phi());
        if (nano.Jet_pt()[ijet]>250f && (dphi > 2.6f || dphi < 0.1f)) goodjet[counter] = false;
        ++counter;
      }
    }
    pico.out_pass_ecalnoisejet() = goodjet[0] && goodjet[1];
  }

  // filters directly from Nano
  pico.out_pass_goodv() = nano.Flag_goodVertices();
  pico.out_pass_cschalo_tight() = nano.Flag_globalSuperTightHalo2016Filter();
  pico.out_pass_hbhe() = nano.Flag_HBHENoiseFilter();
  pico.out_pass_hbheiso() = nano.Flag_HBHENoiseIsoFilter();
  pico.out_pass_ecaldeadcell() = nano.Flag_EcalDeadCellTriggerPrimitiveFilter();
  pico.out_pass_badpfmu() = nano.Flag_BadPFMuonFilter();
  if (nanoaod_version+0.01 > 9) {
    pico.out_pass_badpfmudz() = nano.Flag_BadPFMuonDzFilter();
    pico.out_pass_hfnoisyhits() = nano.Flag_hfNoisyHitsFilter();
  }
  pico.out_pass_eebadsc() = nano.Flag_eeBadScFilter();
  if (year==2016) {
    pico.out_pass_badcalib() = true;
  } else {
    if (nanoaod_version+0.01 < 9) {
      pico.out_pass_badcalib() = nano.Flag_ecalBadCalibFilterV2();
    } else {
      pico.out_pass_badcalib() = nano.Flag_ecalBadCalibFilter();
    }
  }
  pico.out_pass_badchhad() = nano.Flag_BadChargedCandidateFilter();
  pico.out_pass_mubadtrk() = nano.Flag_muonBadTrackFilter();

  if (nanoaod_version+0.01 > 9) { //H->Zy/UL
    //combine pass variable currently consists of recommended JME POG filters
    //as well as optional hfnoisyhits filter
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2 UL section
    pico.out_pass() = pico.out_pass_goodv() && pico.out_pass_cschalo_tight() &&
                      pico.out_pass_hbhe() && pico.out_pass_hbheiso() &&
                      pico.out_pass_ecaldeadcell() && pico.out_pass_badpfmu() &&
                      pico.out_pass_badpfmudz() && pico.out_pass_hfnoisyhits() && 
                      pico.out_pass_eebadsc() && pico.out_pass_badcalib();
    //muon_jet, met/met_calo, met/mht, low_neutral_jet - possibly useful RA2b filters to study
    //htratio_dphi_tight, ecalnoisejet - RA2b filters that may be superseded by hfnoisyhits
    //HEM to study
  }
  else { //HH+MET Run 2/pre-UL
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
                      nano.HLT_Ele28_WPTight_Gsf() || nano.HLT_Ele30_WPTight_Gsf() || nano.HLT_Ele32_WPTight_Gsf() ||
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
  pico.out_HLT_Ele30_WPTight_Gsf() = nano.HLT_Ele30_WPTight_Gsf();
  pico.out_HLT_Ele32_WPTight_Gsf() = nano.HLT_Ele32_WPTight_Gsf();
  pico.out_HLT_Ele32_WPTight_Gsf_L1DoubleEG() = nano.HLT_Ele32_WPTight_Gsf_L1DoubleEG();
  pico.out_HLT_Ele35_WPTight_Gsf() = nano.HLT_Ele35_WPTight_Gsf();
  pico.out_HLT_Ele20_WPLoose_Gsf() = nano.HLT_Ele20_WPLoose_Gsf();
  pico.out_HLT_Ele45_WPLoose_Gsf() = nano.HLT_Ele45_WPLoose_Gsf();
  pico.out_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30() = nano.HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30();
  pico.out_HLT_Ele8_CaloIdM_TrackIdM_PFJet30() = nano.HLT_Ele8_CaloIdM_TrackIdM_PFJet30();
  pico.out_HLT_Ele17_CaloIdM_TrackIdM_PFJet30() = nano.HLT_Ele17_CaloIdM_TrackIdM_PFJet30();
  pico.out_HLT_Ele23_CaloIdM_TrackIdM_PFJet30() = nano.HLT_Ele23_CaloIdM_TrackIdM_PFJet30();
  pico.out_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30() = nano.HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30();
  pico.out_HLT_Ele105_CaloIdVT_GsfTrkIdT() = nano.HLT_Ele105_CaloIdVT_GsfTrkIdT();
  pico.out_HLT_Ele115_CaloIdVT_GsfTrkIdT() = nano.HLT_Ele115_CaloIdVT_GsfTrkIdT();
  pico.out_HLT_Ele135_CaloIdVT_GsfTrkIdT() = nano.HLT_Ele135_CaloIdVT_GsfTrkIdT();
  pico.out_HLT_Ele145_CaloIdVT_GsfTrkIdT() = nano.HLT_Ele145_CaloIdVT_GsfTrkIdT();
  pico.out_HLT_Ele25_eta2p1_WPTight_Gsf() = nano.HLT_Ele25_eta2p1_WPTight_Gsf();
  pico.out_HLT_Ele27_eta2p1_WPTight_Gsf() = nano.HLT_Ele27_eta2p1_WPTight_Gsf();
  pico.out_HLT_Ele28_eta2p1_WPTight_Gsf_HT150() = nano.HLT_Ele28_eta2p1_WPTight_Gsf_HT150();
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

  // Dilepton and diphoton triggers
  bool doubleeg_trigs = nano.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() ||
      nano.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL() || nano.HLT_DoubleEle25_CaloIdL_MW() ||
      nano.HLT_DoublePhoton70() || nano.HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId() ||
      nano.HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55() ||
      nano.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90() ||
      nano.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95();

  bool doublemuon_trigs = nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL() ||
      nano.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL() || nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ() ||
      nano.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ() || nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8() || 
      nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8() || nano.HLT_Mu37_TkMu27();

  bool muoneg_trigs = nano.HLT_Mu17_Photon30_IsoCaloId();

  //Multilepton triggers
  pico.out_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL() = nano.HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL();
  pico.out_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() = nano.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ();
  pico.out_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL()    = nano.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL();
  pico.out_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL()          = nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL();
  pico.out_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL()        = nano.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL();
  pico.out_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ()       = nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ();
  pico.out_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ()     = nano.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ();
  pico.out_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8() = nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8();
  pico.out_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8()   = nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8();
  pico.out_HLT_DoubleEle25_CaloIdL_MW()                = nano.HLT_DoubleEle25_CaloIdL_MW();
  pico.out_HLT_Mu37_TkMu27()                           = nano.HLT_Mu37_TkMu27();
  // Photon triggers
  pico.out_HLT_Mu17_Photon30_IsoCaloId()               = nano.HLT_Mu17_Photon30_IsoCaloId();
  pico.out_HLT_Photon175()                             = nano.HLT_Photon175();
  // Double and di photon trigger
  pico.out_HLT_DoublePhoton70()                        = nano.HLT_DoublePhoton70();
  pico.out_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId() = nano.HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId();
  pico.out_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55() = nano.HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55();
  pico.out_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90() = nano.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90();
  pico.out_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95() = nano.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95();


  //Ele32_WPTight_Gsf was not in menu for most of 2017, but can be emulated by
  //AND of Ele32_WPTight_Gsf_L1DoubleEG and single EG L1 seeds
  if (year==2017) {
    bool pass_l1_singleeg = false;
    for (unsigned itrig = 0; itrig < nano.TrigObj_id().size(); itrig++) {
      if (nano.TrigObj_id()[itrig]==11) {
        if ((nano.TrigObj_filterBits()[itrig] & 0x400)!=0) {
          pass_l1_singleeg = true;
        }
      }
    }
    pico.out_HLT_Ele32_WPTight_Gsf_Emu() = nano.HLT_Ele32_WPTight_Gsf_L1DoubleEG() 
                                           && pass_l1_singleeg;
  }
  else {
    pico.out_HLT_Ele32_WPTight_Gsf_Emu() = false;
  }

  //trigger summary branches for H->Zgamma
  pico.out_trig_single_el() = false;
  pico.out_trig_double_el() = false;
  pico.out_trig_single_mu() = false;
  pico.out_trig_double_mu() = false;

  if (year==2016) {
    pico.out_trig_single_el() = nano.HLT_Ele27_WPTight_Gsf();
    pico.out_trig_double_el() = nano.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ();
    pico.out_trig_single_mu() = nano.HLT_IsoMu24() || nano.HLT_IsoTkMu24();
    pico.out_trig_double_mu() = nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL() ||
                                nano.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL() || 
                                nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ() ||
                                nano.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ();
  }
  if (year==2017) {
    pico.out_trig_single_el() = pico.out_HLT_Ele32_WPTight_Gsf_Emu();
    pico.out_trig_double_el() = nano.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL();
    pico.out_trig_single_mu() = nano.HLT_IsoMu27();
    pico.out_trig_double_mu() = nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8() || 
                                nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8();
  }
  if (year==2018) {
    pico.out_trig_single_el() = nano.HLT_Ele32_WPTight_Gsf();
    pico.out_trig_double_el() = nano.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL();
    pico.out_trig_single_mu() = nano.HLT_IsoMu24();
    pico.out_trig_double_mu() = nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8();
  }
  if (year==2022) {
    pico.out_trig_single_el() = nano.HLT_Ele30_WPTight_Gsf(); 
    pico.out_trig_double_el() = nano.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL();
    pico.out_trig_single_mu() = nano.HLT_IsoMu24();
    pico.out_trig_double_mu() = nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8();
  }
  if (year==2023) {
    //Ele30 is enabled & unprescaled in some of 2023, but it is unclear from
    //POG presentations whether it was enabled until the LHC RQX.L8 incident
    pico.out_trig_single_el() = nano.HLT_Ele30_WPTight_Gsf(); 
    pico.out_trig_double_el() = nano.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL();
    pico.out_trig_single_mu() = nano.HLT_IsoMu24();
    pico.out_trig_double_mu() = nano.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8();
  }

  //overlap removal performed based on trigger decisions
  if (isZgamma) {
    // this assumes that we process all the single and dileptondatasets
    if (dataset==Dataset::DoubleMuon                                                                      && doublemuon_trigs) return true;
    else if (dataset==Dataset::SingleMuon                                                  && muon_trigs && !doublemuon_trigs) return true;
    else if (dataset==Dataset::Muon                                                       && (muon_trigs || doublemuon_trigs)) return true;
    else if (dataset==Dataset::DoubleEG                                 && doubleeg_trigs && !muon_trigs && !doublemuon_trigs) return true;
    else if (dataset==Dataset::SingleElectron          && egamma_trigs && !doubleeg_trigs && !muon_trigs && !doublemuon_trigs) return true;
    else if (dataset==Dataset::EGamma                 && (doubleeg_trigs || egamma_trigs) && !muon_trigs && !doublemuon_trigs) return true;
    else if (dataset==Dataset::MuonEG && muoneg_trigs && !egamma_trigs && !doubleeg_trigs && !muon_trigs && !doublemuon_trigs) return true;
    else return false;
  }
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
  if(Contains(name, "Run20")){ sample = 0;
    if(Contains(name, "SingleElectron") || Contains(name, "EGamma")){ category = 0;
    }else if(Contains(name, "SingleMuon")){ category = 1;
    }else if(Contains(name, "DoubleEG")){   category = 2;
    }else if(Contains(name, "DoubleMuon")){ category = 3;
    }else if(Contains(name, "MET")){        category = 4;
    }else if(Contains(name, "HTMHT")){      category = 5;
    }else if(Contains(name, "JetHT")){      category = 6;
    }
    auto pos = name.find("Run20")+7;
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
  }else if((Contains(name, "TTJets") || Contains(name, "TT_") || Contains(name,"TTTo2L2Nu")) && !Contains(name, "TTTT_")){ sample = 1;
    if(Contains(name, "TTJets_Tune") && !Contains(name,"amcatnlo")){ category = 0; bin = 0;
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
    }else if(Contains(name, "TTTo2L2Nu") || Contains(name,"TTto2L2Nu")){ category = 6; bin = 0;
    }else if(Contains(name, "TTJets_Tune") && Contains(name,"amcatnlo")){ category = 7; bin = 0;
    }else if(Contains(name, "TTtoLNu2Q")){ category = 8; bin = 0;
    }
  }else if(Contains(name, "WJets") && !Contains(name, "TTWJets") && !Contains(name,"ttWJets")){ sample = 2;
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
    }else if(Contains(name, "WtoLNu-2Jets")){ category = 3; bin = 0;
    }
  }else if(Contains(name, "ST_") || Contains(name,"TWminusto") || Contains(name,"TbarWplusto")){ sample = 3;
    if(Contains(name, "ST_s-channel")){ category = 0; bin = 0;
    }else if(Contains(name, "ST_t-channel_top")){ category = 1; bin = 0;
    }else if(Contains(name, "ST_t-channel_antitop")){ category = 2; bin = 0;
    }else if(Contains(name, "ST_tW_top")){ category = 3; bin = 0;
    }else if(Contains(name, "ST_tW_antitop")){ category = 4; bin = 0;
    }else if(Contains(name, "ST_tWA")){ category = 5; bin = 0;
    }else if(Contains(name, "TWminusto2L2Nu")){ category = 6; bin = 0;
    }else if(Contains(name, "TWminustoLNu2Q")){ category = 7; bin = 0;
    }else if(Contains(name, "TbarWplusto2L2Nu")){ category = 8; bin = 0;
    }else if(Contains(name, "TbarWplustoLNu2Q")){ category = 9; bin = 0;
    }
  }else if(Contains(name, "TTWJets") || Contains(name,"ttW")){ sample = 4;
    if(Contains(name, "TTWJetsToLNu")){ category = 0; bin = 0;
    }else if(Contains(name, "TTWJetsToQQ")){ category = 1; bin = 0;
    }else if(Contains(name, "ttWJets_Tune")){ category = 2; bin = 0;
    }else if(Contains(name, "ttWJetsToLNu")){ category = 3; bin = 0;
    }
  }else if(Contains(name, "TTZ") || Contains(name,"ttZ")){ sample = 5;
    if(Contains(name, "TTZToLLNuNu")){ category = 0; bin = 0;
    }else if(Contains(name, "TTZToQQ")){ category = 1; bin = 0;
    }else if(Contains(name, "ttZJets_Tune")){ category = 2; bin = 0;
    }
  }else if(Contains(name, "DYJetsToLL") || Contains(name,"DYto2L")){ sample = 6;
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
    } else if(Contains(name, "DYto2L-2Jets_MLL-50")){ category = 3; bin = 0;
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
  }else if(Contains(name, "ZJets") && !Contains(name,"WZ")){ sample = 8;
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
  }else if(Contains(name, "ttH") && !Contains(name,"HToZG")){ sample = 9;
    if(Contains(name, "ttHJetToGG")){       category = 0; bin = 0;
    }else if(Contains(name, "HToZZ")){      category = 1; bin = 0;
    }else if(Contains(name, "HToWW")){      category = 2; bin = 0;
    }else if(Contains(name, "HToTauTau")){  category = 3; bin = 0;
    }else if(Contains(name, "ttHJetTobb")){ category = 4; bin = 0; //previously category 0
    }else if(Contains(name, "HToMuMu")){    category = 5; bin = 0;
    }
  }else if(Contains(name, "TTGJets")){ sample = 10;
    if(Contains(name, "TTGJets_Tune")){ category = 0; bin = 0;
    }
  }else if(Contains(name, "TTTT")){ sample = 11;
    if(Contains(name, "TTTT_Tune")){ category = 0; bin = 0;
    }
  }else if((Contains(name, "WH_") || Contains(name,"WplusH") || Contains(name,"WminusH") ) && !Contains(name,"TChiWH") && !Contains(name,"HToZG")){ sample = 12;
    if(Contains(name, "HToGG")){           category = 0; bin = 0;
    }else if(Contains(name, "HToZZ")){     category = 1; bin = 0;
    }else if(Contains(name, "HToWW")){     category = 2; bin = 0;
    }else if(Contains(name, "HToTauTau")){ category = 3; bin = 0;
    }else if(Contains(name, "HToBB")){     category = 4; bin = 0; //previously category 0
    }else if(Contains(name, "HToMuMu")){   category = 5; bin = 0;
    }
  }else if(Contains(name, "ZH") && !Contains(name,"HToZG") && !Contains(name, "T5qqqqZH")){ sample = 13;
    if(Contains(name, "HToGG")){           category = 0; bin = 0;
    }else if(Contains(name, "HToZZ")){     category = 1; bin = 0;
    }else if(Contains(name, "HToWW")){     category = 2; bin = 0;
    }else if(Contains(name, "HToTauTau")){ category = 3; bin = 0;
    }else if(Contains(name, "HToBB")){     category = 4; bin = 0; //previously category 0
    }else if(Contains(name, "HToMuMu")){   category = 5; bin = 0;
    }
  }else if(Contains(name, "WW") && !Contains(name,"WWG") && !Contains(name,"WWW") && !Contains(name,"WWZ") && !Contains(name,"HToWW") && !Contains(name,"TChiHH")){ sample = 14;
    if(Contains(name, "WWToLNuQQ")){ category = 0; bin = 0;
    }else if(Contains(name, "WWTo2L2Nu")){ category = 1; bin = 0;
    }else if(Contains(name, "WW_Tune"))  { category = 2; bin = 0;
    }
  }else if(Contains(name, "WZ") && !Contains(name,"WZG") && !Contains(name,"WWZ") && !Contains(name,"WZZ") && !Contains(name,"TChiHH")){ sample = 15;
    if(Contains(name, "WZTo1L3Nu")){ category = 0; bin = 0;
    }else if(Contains(name, "WZTo1L1Nu2Q")){ category = 1; bin = 0;
    }else if(Contains(name, "WZTo2L2Q")){    category = 2; bin = 0;
    }else if(Contains(name, "WZTo3LNu")){    category = 3; bin = 0;
    }else if(Contains(name, "WZ_Tune")){     category = 4; bin = 0;
    }
  }else if(Contains(name, "ZZ") && !Contains(name,"ZZG") && !Contains(name,"WZZ") && !Contains(name,"ZZZ") && !Contains(name,"HToZZ") && !Contains(name,"TChiHH")){ sample = 16;
    if(Contains(name, "ZZ_Tune")){ category = 0; bin = 0;
    }
  }else if((Contains(name, "ZGTo") || Contains(name,"DYG")) && !Contains(name,"ZZG")) { sample = 17; bin = 0;
    if(Contains(name,"ZGTo2LG_Tune")){ category = 0;
    }else if(Contains(name,"ZGToLLG_01J") && !Contains(name,"lowGPt")){ category = 1;
    }else if(Contains(name,"ZGToLLG_01J") && Contains(name,"lowGPt")){ category = 2;
    }else if(Contains(name,"DYGto2LG")) { category = 3;
      if(Contains(name,"PTG-10to50")){         bin = 0;
      }else if(Contains(name,"PTG-50to100")){  bin = 1;
      }else if(Contains(name,"PTG-100to200")){ bin = 2;
      }else if(Contains(name,"PTG-200")){      bin = 3;
      }
    }
  }else if(Contains(name, "TGJets")) { sample = 18; bin = 0;
    if(Contains(name, "TGJets_Tune")) category = 0;
  }else if(Contains(name, "LLAJJ") || Contains(name,"ZGamma2JToGamma2L2J") || Contains(name,"ZG2JtoG2L2J")) { sample = 19; bin = 0;
    if(Contains(name,"EWK_MLL-50")){ category = 0;
    }else if(Contains(name, "ZGamma2JToGamma2L2J") || Contains(name,"ZG2JtoG2L2J")){ category = 1;
    }
  }else if(Contains(name, "EWKZ2Jets")){ sample = 20; category = 0; bin = 0;
  }else if(Contains(name, "WWG")){ sample = 21; category = 0; bin = 0;
  }else if(Contains(name, "WZG")){ sample = 22; category = 0; bin = 0;
  }else if(Contains(name, "ZZG")){ sample = 23; category = 0; bin = 0;
  }else if(Contains(name, "WWW")){ sample = 24;
    if(Contains(name, "WWW_4F_Tune") && !Contains(name,"madspin")){ category = 0; bin = 0;
    }else if(Contains(name, "WWW_4F_DiLeptonFilter")){ category = 1; bin = 0;
    }else if(Contains(name, "WWW_4F_Tune") && Contains(name,"madspin")){ category = 2; bin = 0;
    }
  }else if(Contains(name, "WWZ")){ sample = 25;
    if(Contains(name, "WWZ_4F_Tune")){ category = 0; bin = 0;
    }else if(Contains(name, "WWZJetsTo4L2Nu")){ category = 1; bin = 0;
    }
  }else if(Contains(name, "WZZ")){ sample = 26; category = 0; bin = 0;
  }else if(Contains(name, "ZZZ")){ sample = 27; category = 0; bin = 0;
  }else if(Contains(name, "GluGluH") && !Contains(name,"HToZG")){ sample = 28;
    if(Contains(name, "HToGG")){           category = 0; bin = 0;
    }else if(Contains(name, "HToZZ")){     category = 1; bin = 0;
    }else if(Contains(name, "HToWW")){     category = 2; bin = 0;
    }else if(Contains(name, "HToTauTau")){ category = 3; bin = 0;
    }else if(Contains(name, "HToBB")){     category = 4; bin = 0;
    }else if(Contains(name, "HToMuMu")){   category = 5; bin = 0;
    }
  }else if(Contains(name, "VBFH") && !Contains(name,"HToZG")){ sample = 29;
    if(Contains(name, "HToGG")){           category = 0; bin = 0;
    }else if(Contains(name, "HToZZ")){     category = 1; bin = 0;
    }else if(Contains(name, "HToWW")){     category = 2; bin = 0;
    }else if(Contains(name, "HToTauTau")){ category = 3; bin = 0;
    }else if(Contains(name, "HToBB")){     category = 4; bin = 0;
    }else if(Contains(name, "HToMuMu")){   category = 5; bin = 0;
    }
  }else if(Contains(name, "WGTo") || !Contains(name,"WWG")){ sample = 30; category = 0; bin = 0;
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
    if(code < 0 || code > 207000){
      cout<<"ERROR:: Type code out of range for " << name << ": sample=" << sample << ", category=" << category << ", bin=" << bin;
    }
    return code;
  }
}

void EventTools::WriteTriggerEfficiency(pico_tree &pico) {
  // trigger efficiency and uncertainty - @todo, needs to be updated to full Run 2 trig eff. measurement
  pico.out_w_trig() = hig_trig_eff::eff(pico);
  float effunc = hig_trig_eff::eff_unc(pico);
  pico.out_sys_trig().resize(2,0);
  pico.out_sys_trig()[0] = 1+effunc;
  pico.out_sys_trig()[1] = 1-effunc;

}
