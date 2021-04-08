#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <regex>

#include "TH1D.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TEntryList.h"
#include "TFile.h"
#include "TError.h"
#include "TChainElement.h"

#include "JTreeReaderHelper.h"

string getEntryListName(vector<float> const & POI_masses) {
  string entryListName = "POI_mass";
  for (unsigned iPOI=0; iPOI < POI_masses.size(); ++iPOI) {
    entryListName += "_"+to_string(int(POI_masses[iPOI]));
  }
  return entryListName;
}

vector<string> splitString(string inString, char deliminator) {
  vector<string> split;
  stringstream line(inString);
  string item;
  while(getline(line, item, deliminator)) {
    split.push_back(item);
  }
  return split;
}

// Usage: root "make_split_entrylists.cxx+(input_glob, output_directory, apply_chi2_to_higgs_cut)"
// Makes root files with entrylists according to GenModel branches in NanoAOD.
// root filenames: split_entrylist_GenModel_xx.root, 
// entrylist name: entrylist_GenModel_xx
void make_split_entrylists(string input_glob, string output_directory, int apply_chi2_to_higgs_cut) {
  time_t begtime, endtime;
  time(&begtime);

  TChain * chain = new TChain("Events");
  chain->Add(input_glob.c_str());
  JTreeReaderHelper helper;
  JTreeReader reader(chain);
  map<string, string> branchNameToType = helper.getBranchNameToType(chain);
  reader.setBranches(branchNameToType, &helper);

  // JTreeReader members
  // When looping over trees, need to have access to
  //   - TChain
  //   - firstEntries of trees
  //   - branch names
  //   - storage => JTreeReaderHelper
  //   - currentTTreeReader
  //   - current entry, current tree
  // storage(JTreeReaderHelper) is needed externally.
  // branchNames is needed externally be able to set things dynamically.

  // Loop over branches to get GenModel branches
  set<string> genModelNames;
  for (auto const & it : branchNameToType) {
    if (it.first.find("GenModel") != string::npos) genModelNames.insert(it.first);
  }
  //for (auto const & genModel : genModelNames) cout<<genModel<<endl;

  // POI_entryList[GenModel] = TEntryList
  map< string, TEntryList *> POI_entryLists;

  int n_events_pass_veto = -1;
  int iEvent = 0;
  Long64_t previousEntry=-1;
  Long64_t previousLocalEntry=-1;
  ULong64_t previousEvent=0;
  cout<<"[Start] Loop over events"<<endl;
  while (reader.Next()) {
    Long64_t local_entry = reader.getTreeReader()->GetCurrentEntry();
    Long64_t entry = reader.getEntryNumber();
    //if (local_entry == 0) {
    //  cout<<"entry: "<<entry<<" local entry: "<<local_entry<<" event: "<<helper.s_long64("event");
    //  cout<<" | prev entry: "<<previousEntry<<" local entry: "<<previousLocalEntry<<" event: "<<previousEvent<<endl;
    //}

    // Apply cut to event
    // Confirm that each chi_2 decays to a Higgs
    if (apply_chi2_to_higgs_cut==1) {
      bool is_chi2_to_higgs = false;
      int n_decays_found = 0;
      for (unsigned iGenPart = 0; iGenPart < helper.v_int("GenPart_pdgId").GetSize(); ++iGenPart) {
        if (helper.v_int("GenPart_pdgId")[iGenPart] == 25) {
          int iMother = helper.v_int("GenPart_genPartIdxMother")[iGenPart];
          if (helper.v_int("GenPart_pdgId")[iMother] == 1000023) n_decays_found++;
        }
        if (n_decays_found == 2) {
          is_chi2_to_higgs = true;
          break;
        }
      }
      if (!is_chi2_to_higgs) {
        previousLocalEntry = local_entry;
        previousEntry = entry;
        previousEvent = helper.s_long64("event");
        continue;
      }
    }

    // Find GenModel for event
    string eventGenModelName;
    for (string const & genModelName : genModelNames) {
      if (helper.s_bool(genModelName) == 1) {
        eventGenModelName = genModelName;
        break;
      }
    }

    // Add event to POI_entryLists
    if (POI_entryLists.count(eventGenModelName) == 0) {
      string listName = "entrylist_"+eventGenModelName;
      //cout<<listName<<endl;
      POI_entryLists[eventGenModelName] = new TEntryList(listName.c_str(), listName.c_str());
    }
    POI_entryLists[eventGenModelName]->Enter(entry, chain);
    

    previousLocalEntry = local_entry;
    previousEntry = entry;
    previousEvent = helper.s_long64("event");
    if (n_events_pass_veto == ++iEvent) break;
  }
  cout<<"[End] Loop over events"<<endl;

  //// Print entryLists
  //for (auto & POI_entryList : POI_entryLists) {
  //  cout<<"entryList for "<<POI_entryList.first<<endl;
  //  POI_entryList.second->Print("all");
  //}

  cout<<"[Start] Save entrylists"<<endl;
  // Save entryLists to separate root files
  for (auto & POI_entryList : POI_entryLists) {
    string entryListFilename = "split_entrylist_"+POI_entryList.first+".root";
    string entryListFilepath = output_directory+"/"+entryListFilename;
    TFile eventListFile(entryListFilepath.c_str(), "recreate");
    eventListFile.SetCompressionSettings(209);
    POI_entryList.second->Write();
    cout<<"Wrote file "<<entryListFilepath<<endl;
    eventListFile.Close();
  }
  cout<<"[End] Save entrylists"<<endl;

  time(&endtime);
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
