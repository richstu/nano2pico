#ifndef H_MET_PRODUCER
#define H_MET_PRODUCER


#include "nano_tree.hpp"
#include "pico_tree.hpp"

class MetProducer{
public:

  explicit MetProducer(int year, bool isData, bool verbose=false);
  ~MetProducer();

  void SetVerbose(bool verbose_){ verbose = verbose_; };

  void WriteMet(nano_tree &nano, pico_tree &pico);
private:
  int year;
  bool verbose;
  bool isData;
};

#endif
