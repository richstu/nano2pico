#ifndef TTZ_PRODUCER
#define TTZ_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class TTZVarProducer{
public:

  explicit TTZVarProducer(int year);
  ~TTZVarProducer();

  void WriteTTZVars(pico_tree &pico);

private:
  int year;
  
};

#endif
