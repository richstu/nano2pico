#ifndef H_DILEP_PRODUCER
#define H_DILEP_PRODUCER

#include "pico_tree.hpp"
#include "KinZfitter.hpp"

class DileptonProducer{
public:

  explicit DileptonProducer(int year);
  ~DileptonProducer();

  void WriteDileptons(pico_tree &pico, 
                      std::vector<int> sig_el_pico_idx, std::vector<int> sig_mu_pico_idx);

private:
  int year;
  //KinZfitter *kinZfitter;

  //Below value is used for debug purposes
  //int cnt_refit;
};

#endif
