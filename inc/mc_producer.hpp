#ifndef H_MC_PRODUCER
#define H_MC_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class GenParticleProducer{
public:

  explicit GenParticleProducer(int year);
  ~GenParticleProducer();

  void WriteGenParticles(nano_tree &nano, pico_tree &pico);
  bool IsInteresting(std::vector<int> const & interested_mc_ids, std::vector<std::pair<int, int> > const & interested_mc_ids_range, int mc_id);
  bool IsLastCopyBeforeFSR_or_LastCopy(nano_tree & nano, int mc_index);
  // Makes child_map
  // child_map[mom_index] = [child_index]
  std::map<int, std::vector<int> > GetChildMap(nano_tree & nano);
  int GetFirstCopyIdx(nano_tree & nano, int imc);
  // Returns id to mother not being itself.
  int GetMotherIdx(nano_tree & nano, int imc);
private:
  int year;
  
};

#endif
