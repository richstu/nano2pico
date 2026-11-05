/*!
  \brief Prints the number of entries and branches of a TTree (chained over the
  given path, wildcards allowed), optionally with a selection applied.
  Credit: Claude
  Lightweight cpp implementation to avoid pyROOT import. Not used in batch methods,
  but present as a backup in case batching fails.

  Usage: count_entries.exe [-t tree_name] [-c selection] path [path...]
  Output, one line per path, to stdout:  COUNT <entries> <branches> <path>
  ROOT's own messages (e.g. file not found) go to stderr; a missing file gives 0 0.
*/

#include <iostream>
#include <string>
#include <vector>

#include "TChain.h"

int main(int argc, char *argv[]) {
  std::string tree_name = "tree";
  std::string selection;
  std::vector<std::string> paths;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "-t" && i + 1 < argc) {
      tree_name = argv[++i];
    } else if (arg == "-c" && i + 1 < argc) {
      selection = argv[++i];
    } else {
      paths.push_back(arg);
    }
  }
  if (paths.empty()) {
    std::cerr << "Usage: " << argv[0] << " [-t tree_name] [-c selection] path [path...]" << std::endl;
    return 1;
  }

  for (const std::string &path : paths) {
    TChain chain(tree_name.c_str());
    chain.Add(path.c_str());
    const Long64_t nentries = selection.empty() ? chain.GetEntries() : chain.GetEntries(selection.c_str());
    const Int_t nbranches = chain.GetNbranches();
    std::cout << "COUNT " << nentries << ' ' << nbranches << ' ' << path << std::endl;
  }
  return 0;
}
 

