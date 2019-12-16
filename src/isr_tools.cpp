#include "isr_tools.hpp"

#include <bitset>

#include "TLorentzVector.h"

#include "utilities.hpp"

using namespace std;

ISRTools::ISRTools(const string &name_, int year_):
  name(name_),
  year(year_),
  isTTJets_LO(false),
  isGluino(false),
  isTChi(false){

  if(Contains(name, "TTJets_") && Contains(name, "madgraphMLM")) 
    isTTJets_LO = true;

  if(Contains(name, "SMS-T1") || Contains(name, "SMS-T5")) 
    isGluino = true;

  if(Contains(name, "SMS-TChi")) 
    isTChi = true;
}

ISRTools::~ISRTools(){
}

bool ISRTools::IsLastCopyBeforeFSR_or_LastCopy(nano_tree & nano, int imc){
  bitset<15> mc_statusFlags(nano.GenPart_statusFlags().at(imc));
  int mc_mom_index = nano.GenPart_genPartIdxMother().at(imc);
  int mc_id = nano.GenPart_pdgId().at(imc);

  // 14: isLastCopyBeforeFSR
  // 13: isLastCopy
  if (mc_statusFlags[13] == 1) {
    if (mc_mom_index == -1) return true;
    bitset<15> mom_statusFlags(nano.GenPart_statusFlags().at(mc_mom_index));
    int mom_id = nano.GenPart_pdgId().at(mc_mom_index);
    // A lastCopyBeforeFSR exists
    if (mom_id == mc_id && mom_statusFlags[14] == 1) return false;
    else return true;
  }
  if (mc_statusFlags[14] == 1) return true;
  return false;
}

void ISRTools::WriteISRSystemPt(nano_tree &nano, pico_tree &pico) {
  TLorentzVector isr_p4;
  float mprod(-999), mlsp(-999);
  for(int imc(0); imc<nano.nGenPart(); ++imc) {
    if (IsLastCopyBeforeFSR_or_LastCopy(nano, imc)) {
      int mc_absid = abs(nano.GenPart_pdgId().at(imc));
      //types defined in event tools
      TLorentzVector mc_v4; 
      mc_v4.SetPtEtaPhiM(nano.GenPart_pt()[imc], nano.GenPart_eta()[imc], nano.GenPart_phi()[imc], nano.GenPart_mass()[imc]);
      if (mc_absid==6 && pico.out_type()>=1000 && pico.out_type()<2000) isr_p4 -= mc_v4;
      else if (mc_absid==23 && pico.out_type()>=6000 && pico.out_type()<7000) isr_p4 -= mc_v4;
      else if (pico.out_type()==100e3 || pico.out_type()==102e3 || pico.out_type()==104e3) {
        if (mc_absid==1000021) {
          isr_p4 -= mc_v4;
          mprod = nano.GenPart_mass()[imc];
        } else if (mc_absid==1000022) {
          mlsp = nano.GenPart_mass()[imc];
        }
      } else if ((mc_absid==1000023 || mc_absid==1000025) && pico.out_type()==106e3) {
        isr_p4 -= mc_v4;
        mprod = nano.GenPart_mass()[imc];
        mlsp = 1.;
      }
    }
  }
  pico.out_isr_tru_pt() = isr_p4.Pt();
  pico.out_isr_tru_eta() = isr_p4.Eta();
  pico.out_isr_tru_phi() = isr_p4.Phi();
  pico.out_mprod() = mprod;
  pico.out_mlsp() = mlsp;

  return;
}

std::map<int, std::vector<int> > ISRTools::GetChildMap(nano_tree & nano)
{
  map<int, vector<int> > child_map;
  for(int imc(0); imc<nano.nGenPart(); ++imc) {
    int mc_mom_index = nano.GenPart_genPartIdxMother().at(imc);
    child_map[mc_mom_index].push_back(imc);
  }
  return child_map;
}

void ISRTools::WriteISRJetMultiplicity(nano_tree &nano, pico_tree &pico) {
  bool verbose = false;
  map<int, vector<int> > child_map = GetChildMap(nano);

  pico.out_nisr() = 0;
  for (size_t ijet(0); ijet<pico.out_jet_pt().size(); ijet++){
    if (!pico.out_jet_isgood()[ijet]) continue;

    bool matched = false;
    for (int imc(0); imc < nano.nGenPart(); imc++) {
      if (matched) break;
      if (nano.GenPart_status()[imc]!=23 || abs(nano.GenPart_pdgId()[imc])>5) continue;
      int momid = abs(nano.GenPart_pdgId()[nano.GenPart_genPartIdxMother()[imc]]);
      if(!(momid==6 || momid==23 || momid==24 || momid==25 || momid>1e6)) continue; 
      //check against daughter in case of hard initial splitting
      for (size_t idau(0); idau < child_map[imc].size(); idau++) {
        float dr_ = dR(pico.out_jet_eta()[ijet], nano.GenPart_eta()[child_map[imc][idau]], 
                      pico.out_jet_phi()[ijet], nano.GenPart_phi()[child_map[imc][idau]]);
        if(dr_<0.3){
          if (verbose) cout<<"Jet: ("<<pico.out_jet_pt()[ijet]<<", "<<pico.out_jet_eta()[ijet]<<", "
                           <<pico.out_jet_phi()[ijet]<<"), MC: ("<<nano.GenPart_pt()[child_map[imc][idau]]
                           <<", "<<nano.GenPart_eta()[child_map[imc][idau]]<<", "
                           <<nano.GenPart_phi()[child_map[imc][idau]]<<"), ID "
                           <<nano.GenPart_pdgId()[child_map[imc][idau]]<<". dr "<<dr_ <<endl;
            matched = true;
            break;
        }
      }
    } // Loop over MC particles
    if(!matched) {
      pico.out_nisr()++;
    }
  } // Loop over jets
}

void ISRTools::WriteISRWeights(pico_tree &pico) {
  pico.out_w_isr() = 1.; pico.out_sys_isr().resize(2,1.);
  if (isGluino || (year==2016 && isTTJets_LO)){
    const float isr_norm_tt = 1.117;
    float isr_wgt = -999.;
    if (pico.out_nisr()==0)      isr_wgt = 1.; 
    else if (pico.out_nisr()==1) isr_wgt = 0.920; 
    else if (pico.out_nisr()==2) isr_wgt = 0.821; 
    else if (pico.out_nisr()==3) isr_wgt = 0.715; 
    else if (pico.out_nisr()==4) isr_wgt = 0.662; 
    else if (pico.out_nisr()==5) isr_wgt = 0.561; 
    else if (pico.out_nisr()>=6) isr_wgt = 0.511; 
    pico.out_w_isr() = isr_wgt*isr_norm_tt;
    //assign relative unc = 50% of the deviation from flat
    float absolute_unc = (1-isr_wgt)/2.;
    pico.out_sys_isr()[0] = isr_norm_tt*(isr_wgt+absolute_unc); 
    pico.out_sys_isr()[1] = isr_norm_tt*(isr_wgt-absolute_unc); 
  } else if (isTChi) {
    float isr_wgt = 1.;
    if      (pico.out_isr_tru_pt()<=50)  isr_wgt = 1.;
    else if (pico.out_isr_tru_pt()<=100) isr_wgt = 1.052;
    else if (pico.out_isr_tru_pt()<=150) isr_wgt = 1.179;
    else if (pico.out_isr_tru_pt()<=200) isr_wgt = 1.150;
    else if (pico.out_isr_tru_pt()<=300) isr_wgt = 1.057;
    else if (pico.out_isr_tru_pt()<=400) isr_wgt = 1.000;
    else if (pico.out_isr_tru_pt()<=600) isr_wgt = 0.912;
    else                             isr_wgt = 0.783; 
    pico.out_w_isr() = isr_wgt;
    //assign relative unc = 100% of the deviation from flat
    if (isr_wgt>1) pico.out_sys_isr()[0] = 1+2*(isr_wgt-1);
    else pico.out_sys_isr()[0] = 1-2*(1-isr_wgt);
    pico.out_sys_isr()[1] = 1.;
  }
  return;
}