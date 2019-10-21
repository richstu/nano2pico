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

  int GetFirstCopyIdx(nano_tree & nano, int imc);
  // Returns id to mother not being itself.
  int GetMotherIdx(nano_tree & nano, int imc);
private:
  int year;
  
};

#endif
