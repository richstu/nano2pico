#ifndef H_MC_PRODUCER
#define H_MC_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class GenParticleProducer{
public:

  explicit GenParticleProducer(int year,float nanoaod_version);
  ~GenParticleProducer();

  void WriteGenParticles(nano_tree &nano, pico_tree &pico, bool isDY);
  bool IsInteresting(std::vector<int> const & interested_mc_ids, std::vector<std::pair<int, int> > const & interested_mc_ids_range, int mc_id);

  int GetFirstCopyIdxOrInterestingIdx(nano_tree & nano, int imc, std::map<int, int> & mc_index_to_interested_index);
  // Returns id to mother not being itself.
  int GetMotherIdx(nano_tree & nano, int imc, std::map<int, int> & mc_index_to_interested_index);
private:
  int year;
  float nanoaod_version;
};

#endif
