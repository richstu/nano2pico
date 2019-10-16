#include "mc_producer.hpp"
#define BACKWARD_HAS_BFD 0
#define BACKWARD_HAS_DWARF 0
#define BACKWARD_HAS_BACKTRACE_SYMBOL 0
#define BACKWARD_HAS_UNWIND 0
#define BACKWARD_HAS_BACKTRACE 0
#define BACKWARD_HAS_DW 0
// #include "backward.h"
#include "utilities.hpp"
#include <bitset>
#include <algorithm>

// namespace backward {
//   backward::SignalHandling sh;
// } // namespace backward

using namespace std;

GenParticleProducer::GenParticleProducer(int year_){
  year = year_;
}

GenParticleProducer::~GenParticleProducer(){
}

bool GenParticleProducer::IsInteresting(vector<int> const & interested_mc_ids, vector<pair<int, int> > const & interested_mc_ids_range, int mc_id){
  if (find(interested_mc_ids.begin(), interested_mc_ids.end(), mc_id)!=interested_mc_ids.end()) return true;
  if (find(interested_mc_ids.begin(), interested_mc_ids.end(), -1*mc_id)!=interested_mc_ids.end()) return true;
  for (auto id_range : interested_mc_ids_range) {
    int id_range_start = id_range.first;
    int id_range_end = id_range.second;
    if ((fabs(mc_id) >= fabs(id_range_start)) && (fabs(mc_id) <= fabs(id_range_end))) return true;
  }
  return false;
}

bool GenParticleProducer::IsLastCopyBeforeFSR_or_LastCopy(nano_tree & nano, int imc){
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

std::map<int, std::vector<int> > GenParticleProducer::GetChildMap(nano_tree & nano)
{
  map<int, vector<int> > child_map;
  for(int imc(0); imc<nano.nGenPart(); ++imc) {
    int mc_mom_index = nano.GenPart_genPartIdxMother().at(imc);
    child_map[mc_mom_index].push_back(imc);
  }
  return child_map;
}

int GenParticleProducer::GetFirstCopyIdx(nano_tree & nano, int imc)
{
  if (imc <= 0) {
    return imc;
  }
  int mc_id = nano.GenPart_pdgId().at(imc);
  int mc_mom_index = nano.GenPart_genPartIdxMother().at(imc);
  if (mc_mom_index == -1) {
    // Found that initial particles isFirstCopy is 0. Do not use below check.
    //if (bitset<15>(nano.GenPart_statusFlags().at(imc))[12] != 1) 
    //{
    //  cout<<"[Error] 1 GenParticleProducer::GetFirstCopyIdx is not firstCopy."<<endl;
    //  cout<<"  imc:"<<imc<<" mc_id:"<<mc_id<<endl;
    //}
    return imc;
  }
  int mc_mom_id = nano.GenPart_pdgId().at(mc_mom_index);
  //cout<<mc_id<<":"<<mc_mom_index<<" mc_id:"<<mc_id<<" mc_mom_id: "<<mc_mom_id<<endl;
  if (mc_id == mc_mom_id) return GetFirstCopyIdx(nano, mc_mom_index);
  if (bitset<15>(nano.GenPart_statusFlags().at(imc))[12] != 1) 
  {
    cout<<"[Error] GenParticleProducer::GetFirstCopyIdx is not firstCopy."<<endl;
    cout<<"  imc:"<<imc<<" mc_id:"<<mc_id<<endl;
  }
  return imc;
}

int GenParticleProducer::GetMotherIdx(nano_tree & nano, int imc)
{
  if (imc == 0) return -1;
  int mc_mom_index = nano.GenPart_genPartIdxMother().at(imc);
  return GetFirstCopyIdx(nano, mc_mom_index);
}

void GenParticleProducer::WriteGenParticles(nano_tree &nano, pico_tree &pico){

  // Saves if isPrompt and isFirstCopy
  // H, Z, W, b, t, 
  vector<int> interested_mc_ids = {25, 23, 24, 5, 6};
  // SUSY
  vector<pair<int, int> > interested_mc_ids_range = {{1000001, 2000015}};

  // Saves if isPrompt and isFirstCopy
  // Saves if isTauDecayProduct and isFirstCopy
  // e, mu, tau
  vector<int> interested_lepton_ids = {11, 13, 15};
  int ntrulep = 0;
  int ntrumu = 0;
  int ntruel = 0;
  int ntrutau = 0;
  int ntrutaul = 0;
  int ntrutauh = 0;

  // Collect interesting particle indices
  vector<int> interested_mc_indices;
  for(int imc(0); imc<nano.nGenPart(); ++imc) {
    int mc_id = nano.GenPart_pdgId().at(imc);
    bitset<15> mc_statusFlags(nano.GenPart_statusFlags().at(imc));
    bool is_interesting = IsInteresting(interested_mc_ids, interested_mc_ids_range, mc_id);
    bool lepton_interesting = IsInteresting(interested_lepton_ids, {}, mc_id);
    bool save_index = false;
    bool is_tauDecayProduct = false;
    if (is_interesting) {
      // 0: isPrompt, 12: isFirstCopy, 
      if (mc_statusFlags[0]==1 && mc_statusFlags[12]==1) save_index = true;
    }
    if (lepton_interesting) {
      // 0: isPrompt, 12: isFirstCopy, 
      if (mc_statusFlags[0]==1 && mc_statusFlags[12]==1) save_index = true;
      // 2: isTauDecayProduct, 12: isFirstCopy, 
      if (mc_statusFlags[2]==1 && mc_statusFlags[12]==1) {
        save_index = true;
        is_tauDecayProduct = true;
      }
    }
    // store information
    if (save_index) {
      interested_mc_indices.push_back(imc);
      if (abs(mc_id) == 11 && !is_tauDecayProduct) ntruel++;
      if (abs(mc_id) == 13 && !is_tauDecayProduct) ntrumu++;
      if (abs(mc_id) == 11 && is_tauDecayProduct) ntrutaul++;
      if (abs(mc_id) == 13 && is_tauDecayProduct) ntrutaul++;
      if (abs(mc_id) == 15) ntrutau++;
    }
  }
  ntrulep = ntrumu + ntruel;
  ntrutauh = ntrutau - ntrutaul;

  // Find relation between indices
  // mc_index_to_interested_index[imc] = interested_index
  map<int, int> mc_index_to_interested_index;
  for (unsigned interested_index=0; interested_index<interested_mc_indices.size() ; ++interested_index) {
    int imc = interested_mc_indices[interested_index];
    mc_index_to_interested_index[imc] = interested_index;
    //cout<<"imc: "<<imc<<" idx: "<<interested_index<<endl;
  }

  // Save interesting particles
  for (auto imc : interested_mc_indices) {
    // Parse information
    float mc_pt = nano.GenPart_pt().at(imc);
    float mc_eta = nano.GenPart_eta().at(imc);
    float mc_phi = nano.GenPart_phi().at(imc);
    float mc_mass = nano.GenPart_mass().at(imc);
    int mc_id = nano.GenPart_pdgId().at(imc);
    // Find index of mother that is not itself.
    //int mc_mom_index = nano.GenPart_genPartIdxMother().at(imc);
    int mc_mom_index = GetMotherIdx(nano, imc);
    //cout<<"imc: "<<imc<<" mc_mom_index: "<<mc_mom_index<<endl;
    int mc_mom = mc_mom_index>=0 ? nano.GenPart_pdgId().at(mc_mom_index): -1;
    int mc_mom_idx = -1;
    if (mc_index_to_interested_index.find(mc_mom_index) != mc_index_to_interested_index.end()) {
      mc_mom_idx = mc_index_to_interested_index[mc_mom_index];
    }
    int mc_status = nano.GenPart_status().at(imc);

    // Save information
    pico.out_mc_pt().push_back(mc_pt);
    pico.out_mc_eta().push_back(mc_eta);
    pico.out_mc_phi().push_back(mc_phi);
    pico.out_mc_mass().push_back(mc_mass);
    pico.out_mc_id().push_back(mc_id);
    pico.out_mc_mom().push_back(mc_mom);
    pico.out_mc_momidx().push_back(mc_mom_idx);
    pico.out_mc_status().push_back(mc_status);
    pico.out_ntrulep() = ntrulep;
    pico.out_ntruel() = ntruel;
    pico.out_ntrumu() = ntrumu;
    pico.out_ntrutaul() = ntrutaul;
    pico.out_ntrutauh() = ntrutauh;
  }

  return;
}
