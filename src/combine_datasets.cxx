// combine_datasets: Finds all unique events in a list of datasets
// ported from babymaker

#include <ctime>

#include <vector>
#include <fstream>
#include <iostream>
#include <set>
#include <map>
#include <unistd.h>  // getopt
#include <iomanip>   // setw

#include "TString.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"

#include "utilities.hpp"

using namespace std; 

bool same_era(int run_a, int run_b) {
  if (run_a >= 272007 && run_a <= 275376) {
    if (run_b >= 272007 && run_b <= 275376) {
      return true;
    }
  }
  else if (run_a >= 275657 && run_a <= 276283) {
    if (run_b >= 275657 && run_b <= 276283) {
      return true;
    }
  }
  else if (run_a >= 276315 && run_a <= 276811) {
    if (run_b >= 276315 && run_b <= 276811) {
      return true;
    }
  }
  else if (run_a >= 276831 && run_a <= 277420) {
    if (run_b >= 276831 && run_b <= 277420) {
      return true;
    }
  }
  else if (run_a >= 277772 && run_a <= 278808) {
    if (run_b >= 277772 && run_b <= 278808) {
      return true;
    }
  }
  else if (run_a >= 278820 && run_a <= 280385) {
    if (run_b >= 278820 && run_b <= 280385) {
      return true;
    }
  }
  else if (run_a >= 280919 && run_a <= 284044) {
    if (run_b >= 280919 && run_b <= 284044) {
      return true;
    }
  }
  return false;
}

TString get_era(int run) {
  if (run >= 272007 && run <= 275376) {
    return TString("Run2016B");
  }
  else if (run >= 275657 && run <= 276283) {
    return TString("Run2016C");
  }
  else if (run >= 276315 && run <= 276811) {
    return TString("Run2016D");
  }
  else if (run >= 276831 && run <= 277420) {
    return TString("Run2016E");
  }
  else if (run >= 277772 && run <= 278808) {
    return TString("Run2016F");
  }
  else if (run >= 278820 && run <= 280385) {
    return TString("Run2016G");
  }
  else if (run >= 280919 && run <= 284044) {
    return TString("Run2016H");
  }
  return TString("Error");
}

int main(int argc, char *argv[]){
  time_t startTime, curTime;
  time(&startTime);

  TString file_datasets("txt/combine_datasets/singlelep.txt"), infolder(""), outfolder("out/"),basename(""),yearname("");
  int begrun(-1), endrun(-1);
  
  int c(0);
  while((c=getopt(argc, argv, "f:i:o:n:b:e:y:"))!=-1){
    switch(c){
    case 'i':
      infolder=optarg;
      break;
    case 'b':
      begrun=atoi(optarg);
      break;
    case 'e':
      endrun=atoi(optarg);
      break;
    case 'o':
      outfolder=optarg;
      break;
    case 'n':
      basename=optarg;
      break;
    case 'f':
      file_datasets=optarg;
      break;
    case 'y':
      yearname=optarg;
      break;
    default:
      break;
    }
  }
  if(file_datasets=="" || infolder==""){
    cout<<endl<<"Specify input folder and datasets: "
        <<"./run/combine_datasets.exe -i <infolder> -o <outfolder=out> -f <file_datasets=txt/combine_datasets/singlelep.txt> -b  <begrun=-1> -e <endrun=-1> -n <basename>"<<endl<<endl;
    return 1;
  }

  //make run portion of output name and check that run range is sensible
  //begrun=endrun=-1 uses all runs
  //endrun=-1 uses only begrun
  bool is_single_era = same_era(begrun, endrun);
  TString run_s="_runs"; run_s += begrun; 
  if(endrun>begrun){
    run_s += "-"; run_s += endrun;
  }
  if(begrun>0){
    if(endrun<begrun){
      cout<<"You set begrun to "<<begrun<<", and endrun to "<<endrun
          <<", but endrun has to be >= to begrun. Exiting"<<endl<<endl;
      return 1;
    }
    cout<<"Combining "<<run_s<<" of ntuples in "<<infolder<<endl;
  }

  //add names of datasets to variable datasets and construct basename for output files
  vector<TString> datasets;
  TString buffer;
  ifstream indata(file_datasets);
  while(indata){
    indata >> buffer;
    if(buffer!=""){
      datasets.push_back(buffer);
    }
  }
  for (int dataset_idx = datasets.size() - 1; dataset_idx >= 0; dataset_idx--) {
      basename = (datasets[dataset_idx]+"_"+basename);
  }

  map<int, map<int, set<Long64_t> > > runs;
  Long64_t event;
  int run, lumiblock;

  for(unsigned idata(0); idata < datasets.size(); idata++){
    //make chain of input files
    TChain chain("tree");
    TString filename(infolder+"/*"+datasets[idata]+"*.root");
    if (is_single_era) {
            filename = (infolder+"/*"+datasets[idata]+"*"+get_era(begrun)+"*.root");
    }
    int files = chain.Add(filename);
    if(files<1) {
      cout<<"No files found for "<<filename<<endl;
      continue;
    }

    //generate output files
    //gSystem->mkdir(outfolder, kTRUE);
    TString outname(outfolder+"/pico_");
    if (is_single_era) 
      outname += get_era(begrun) + "_";
    else
    outname += yearname + "_";
    outname += idata;
    outname += "_"+basename;
    if (begrun > 0) outname += run_s;
    outname += ".root";
    TFile outfile(outname, "RECREATE");
    outfile.cd();

    TTree *outtree(chain.CloneTree(0));

    // TBranch *b_event = chain.Branch("event", &event);
    // TBranch *b_run = chain.Branch("run", &run);
    TBranch *b_event(nullptr), *b_lumiblock(nullptr), *b_run(nullptr);
    chain.SetBranchAddress("event", &event, &b_event);
    chain.SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
    chain.SetBranchAddress("run", &run, &b_run);

    long entries(chain.GetEntries()), tree_entry;

    cout<<endl<<"Doing "<<files<<" files in "<<filename<<" with "<<entries<<" entries"<<endl;
    time(&startTime);
    for(int entry(0); entry<entries; entry++){
      if(entry!=0 && entry%250000==0) {
        time(&curTime);
        int seconds(difftime(curTime,startTime));
        
        cout<<"Doing entry "<<setw(10)<<addCommas(entry)<<" of "<<addCommas(entries)
            <<"    Took "<<setw(6)<<seconds<<" seconds at "
            <<setw(4)<<roundNumber(entry,1,seconds*1000.)<<" kHz"<<endl;
      }
      
      // Load "run" first, and check if it's in the range we care about
      tree_entry = chain.LoadTree(entry);
      b_run->GetEntry(tree_entry);
      if(begrun>0 && (run<begrun || run>endrun)) continue;
      b_lumiblock->GetEntry(tree_entry);
      b_event->GetEntry(tree_entry);

      if(runs.find(run) == runs.end()) runs.emplace(run, map<int, set<Long64_t> >{}); // New run
      auto &lumiblocks = runs.at(run);
      if(lumiblocks.find(lumiblock) == lumiblocks.end()) lumiblocks.emplace(lumiblock, set<Long64_t>{}); // New lumiblock
      auto &events = lumiblocks.at(lumiblock);
      if(events.find(event) == events.end()){ // New event
        events.emplace(event);
        // You need to load all branches to copy them into outtree
        chain.GetEntry(entry);
        outtree->Fill();
      } 
    } // Loop over entries
    outtree->Write();
    outfile.Write();
    outfile.Close();
    time(&curTime);
    cout<<"Took "<<difftime(curTime,startTime) <<" seconds to write "<<outname<<endl;

  } // Loop over datasets

  // for(auto it = events.cbegin(); it != events.cend(); ++it) {
  //   cout << it->first  <<", ";
  // } // Needs c++11

  //if(false){
  //  TString txtname(outfolder+"/runs_"+basename+".txt");
  //  ofstream txtfile(txtname);
  //  int prevrun(0);
  //  for(map<int, map<int, set<Long64_t> > >::const_iterator it = runs.begin(); it != runs.end(); ++it) {
  //    run = it->first;
  //    if(run/1000 != prevrun){
  //      prevrun = run/1000;
  //      txtfile<<endl;
  //    }
  //    txtfile << run << "  ";
  //  }
  //  txtfile<<endl;
  //  txtfile.close();
  //  cout<<endl<<"Written run numbers in "<<txtname<<endl;
  //}
  cout<<endl<<endl;

  return 0;
}
