#ifndef H_DILEP_PRODUCER
#define H_DILEP_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class DileptonProducer{
public:

  explicit DileptonProducer(int year);
  ~DileptonProducer();

  void WriteDielectrons(nano_tree &nano, pico_tree &pico, std::vector<int> sig_el_nano_idx);
  void WriteDimuons(nano_tree &nano, pico_tree &pico, std::vector<int> sig_mu_nano_idx);
  void WriteDileptons(nano_tree &nano, pico_tree &pico, 
                      std::vector<int> sig_el_nano_idx, std::vector<int> sig_mu_nano_idx);

private:
  int year;
  
};

#endif
