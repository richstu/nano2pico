#ifndef H_MC_PRODUCER
#define H_MC_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class GenParticleProducer{
public:

  explicit GenParticleProducer(int year);
  ~GenParticleProducer();

  void WriteGenParticles(nano_tree &nano, pico_tree &pico);

private:
  int year;
  
};

#endif