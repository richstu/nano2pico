#ifndef H_DILEP_PRODUCER
#define H_DILEP_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class DileptonProducer{
public:

  explicit DileptonProducer(int year);
  ~DileptonProducer();

  void WriteDileptons();

private:
  int year;
  
};

#endif