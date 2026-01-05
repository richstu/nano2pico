#ifndef H_DILEP_PRODUCER
#define H_DILEP_PRODUCER

#include "pico_tree.hpp"

class DileptonProducer{
public:

  explicit DileptonProducer(int year);
  ~DileptonProducer();

  void WriteDileptons(pico_tree &pico, bool is_signal);

private:
  int year;
};

#endif
