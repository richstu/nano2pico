#ifndef H_HIG_PRODUCER
#define H_HIG_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class HigVarProducer{
public:

  explicit HigVarProducer(int year);
  ~HigVarProducer();

  void WriteHigVars(pico_tree& pico, bool doDeepFlav);
  void WriteDPhiVars();

private:
  int year;
  
};

#endif