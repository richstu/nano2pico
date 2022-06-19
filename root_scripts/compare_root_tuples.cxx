#include <iostream>
#include <vector>
#include <string>
#include <ctime>

#include "TChain.h"

#include "JTreeReaderHelper.h"

// Usage: root "compare_root_tuples.cxx+(root_filename_a, root_filename_b, tree_name)"
void compare_root_tuples(string root_filename_a="", string root_filename_b="", string tree_name="tree") {
  time_t begtime, endtime;
  time(&begtime);

  if (root_filename_a == "" || root_filename_b == "") {
    cout<<"Usage: root \"compare_root_tuples.cxx+(\\\"root_filename_a\\\", \\\"root_filename_b\\\", \\\"tree_name\\\")\""<<endl;
    return;
  }

  //string root_filename_a = "/homes/jbkim/analysis/nano2pico/unit_test_higgsino/mc/unskimmed/pico_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv7__PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1__120000__1E421EBC-226E-DC47-A7CE-BFBCC3760D67.root";
  //string root_filename_b = "/homes/jbkim/analysis/nano2pico.variables/unit_test_higgsino/mc/unskimmed/pico_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv7__PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1__120000__1E421EBC-226E-DC47-A7CE-BFBCC3760D67.root";
  //string tree_name = "tree";

  TChain * chain_a = new TChain(tree_name.c_str());
  chain_a->Add(root_filename_a.c_str());
  JTreeReaderHelper helper_a;
  JTreeReader reader_a(chain_a);
  map<string, string> branchNameToType_a = helper_a.getBranchNameToType(chain_a);
  reader_a.setBranches(branchNameToType_a, &helper_a);

  TChain * chain_b = new TChain(tree_name.c_str());
  chain_b->Add(root_filename_b.c_str());
  JTreeReaderHelper helper_b;
  JTreeReader reader_b(chain_b);
  map<string, string> branchNameToType_b = helper_b.getBranchNameToType(chain_b);
  reader_b.setBranches(branchNameToType_b, &helper_b);

  // Get all branch names
  set<string> all_branch_names;
  for (auto itBranch : branchNameToType_a) all_branch_names.insert(itBranch.first);
  for (auto itBranch : branchNameToType_b) all_branch_names.insert(itBranch.first);
  // Compare branch names and types
  set<string> missing_branches_a;
  set<string> missing_branches_b;
  set<string> common_branch_names;
  for (auto branch_name : all_branch_names) {
    bool is_branch_a = branchNameToType_a.find(branch_name) != branchNameToType_a.end();
    bool is_branch_b = branchNameToType_b.find(branch_name) != branchNameToType_b.end();
    if (!is_branch_a) missing_branches_a.insert(branch_name);
    if (!is_branch_b) missing_branches_b.insert(branch_name);
    if (is_branch_a && is_branch_b) {
      if (branchNameToType_a[branch_name] != branchNameToType_a[branch_name]) {
        cout<<"[Info] type of branch "<<branch_name<<" are different. ("<<branchNameToType_a[branch_name]<<","<<branchNameToType_a[branch_name]<<"). Will not be compared."<<endl;
      } else {
        common_branch_names.insert(branch_name);
      }
    }
  }
  // Print missing branches
  if (missing_branches_a.size() !=0) cout<<"[Info] root file ("<<root_filename_a<<") is missing following branches: ";
  for (auto branch_name : missing_branches_a) cout<<branch_name<<", ";
  if (missing_branches_a.size() !=0) cout<<endl;
  if (missing_branches_b.size() !=0) cout<<"[Info] root file ("<<root_filename_b<<") is missing following branches: ";
  for (auto branch_name : missing_branches_b) cout<<branch_name<<", ";
  if (missing_branches_b.size() !=0) cout<<endl;

  // different_branches[branch_name] = [(iEntry, iValue)], where iValue is -1 for scalar case and -2 is for vector case where size is different
  map<string, vector< tuple<Long64_t, int> > > different_branches;

  // Compare values for each event
  while (1) {
    bool is_next_a = reader_a.Next();
    bool is_next_b = reader_b.Next();
    if ((is_next_a && is_next_b)==0) break;

    //Long64_t local_entry_a = reader_a.getTreeReader()->GetCurrentEntry();
    Long64_t entry_a = reader_a.getEntryNumber();
    //cout<<"local entry a: "<<local_entry_a<<" global entry a: "<<entry_a<<endl;
    //Long64_t local_entry_b = reader_b.getTreeReader()->GetCurrentEntry();
    //Long64_t entry_b = reader_b.getEntryNumber();
    //cout<<"local entry b: "<<local_entry_b<<" global entry b: "<<entry_b<<endl;

    for (auto branch_name : common_branch_names) { // Loop over branches
      //if (branch_name != "ll_charge") continue;
      //cout<<"Comapre "<<branch_name<<" "<<branchNameToType_a[branch_name]<<endl;
      if (branchNameToType_a[branch_name] == "Bool_t") { // Compare Bool_t 
        if (!std::isnan(helper_a.s_bool(branch_name)) || !std::isnan(helper_b.s_bool(branch_name))) {
          if (helper_a.s_bool(branch_name) != helper_b.s_bool(branch_name)) {
            if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<" has different values: "<<helper_a.s_bool(branch_name)<<" "<<helper_b.s_bool(branch_name)<<endl;
            different_branches[branch_name].push_back({entry_a, -1});
          }
        }
      } else if (branchNameToType_a[branch_name] == "UInt_t") { // Compare UInt_t
        if (!std::isnan(helper_a.s_uint(branch_name)) || !std::isnan(helper_b.s_uint(branch_name))) {
          if (helper_a.s_uint(branch_name) != helper_b.s_uint(branch_name)) {
            if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<" has different values: "<<helper_a.s_uint(branch_name)<<" "<<helper_b.s_uint(branch_name)<<endl;
            different_branches[branch_name].push_back({entry_a, -1});
          }
        }
      } else if (branchNameToType_a[branch_name] == "ULong64_t") { // Compare ULong64_t
        if (!std::isnan(helper_a.s_ulong64(branch_name)) || !std::isnan(helper_b.s_ulong64(branch_name))) {
          if (helper_a.s_ulong64(branch_name) != helper_b.s_ulong64(branch_name)) {
            if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<" has different values: "<<helper_a.s_ulong64(branch_name)<<" "<<helper_b.s_ulong64(branch_name)<<endl;
            different_branches[branch_name].push_back({entry_a, -1});
          }
        }
      } else if (branchNameToType_a[branch_name] == "Long64_t") { // Compare Long64_t
        if (!std::isnan(helper_a.s_long64(branch_name)) || !std::isnan(helper_b.s_long64(branch_name))) {
          if (helper_a.s_long64(branch_name) != helper_b.s_long64(branch_name)) {
            if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<" has different values: "<<helper_a.s_long64(branch_name)<<" "<<helper_b.s_long64(branch_name)<<endl;
            different_branches[branch_name].push_back({entry_a, -1});
          }
        }
      } else if (branchNameToType_a[branch_name] == "Int_t") { // Compare Int_t
        if (!std::isnan(helper_a.s_int(branch_name)) || !std::isnan(helper_b.s_int(branch_name))) {
         if (helper_a.s_int(branch_name) != helper_b.s_int(branch_name)) {
           if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<" has different values: "<<helper_a.s_int(branch_name)<<" "<<helper_b.s_int(branch_name)<<endl;
           different_branches[branch_name].push_back({entry_a, -1});
         }
        }
      } else if (branchNameToType_a[branch_name] == "Float_t") { // Compare Float_t
        if (!std::isnan(helper_a.s_float(branch_name)) || !std::isnan(helper_b.s_float(branch_name))) {
          if (helper_a.s_float(branch_name) != helper_b.s_float(branch_name)) {
            if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<" has different values: "<<helper_a.s_float(branch_name)<<" "<<helper_b.s_float(branch_name)<<endl;
            different_branches[branch_name].push_back({entry_a, -1});
          }
        }
      } else if (branchNameToType_a[branch_name] == "Double_t") { // Compare Double_t
        if (!std::isnan(helper_a.s_double(branch_name)) || !std::isnan(helper_b.s_double(branch_name))) {
          if (helper_a.s_double(branch_name) != helper_b.s_double(branch_name)) {
            if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<" has different values: "<<helper_a.s_double(branch_name)<<" "<<helper_b.s_double(branch_name)<<endl;
            different_branches[branch_name].push_back({entry_a, -1});
          }
        }
      } else if (branchNameToType_a[branch_name] == "vector<bool>") { // Compare vector<bool>
        if (helper_a.v_bool(branch_name).GetSize() != helper_b.v_bool(branch_name).GetSize()) { // Compare size
          if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<" has different sizes: "<<helper_a.v_bool(branch_name).GetSize()<<" "<<helper_b.v_bool(branch_name).GetSize()<<endl;
          different_branches[branch_name].push_back({entry_a, -2});
        } else { // Compare values
          for (unsigned iValue=0; iValue < helper_a.v_bool(branch_name).GetSize(); ++iValue) {
            if (!std::isnan(helper_a.v_bool(branch_name)[iValue]) || !std::isnan(helper_b.v_bool(branch_name)[iValue])) {
              if (helper_a.v_bool(branch_name)[iValue] != helper_b.v_bool(branch_name)[iValue]) {
                if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<"["<<iValue<<"] has different values: "<<helper_a.v_bool(branch_name)[iValue]<<" "<<helper_b.v_bool(branch_name)[iValue]<<endl;
                different_branches[branch_name].push_back({entry_a, iValue});
              }
            }
          }
        }
      } else if (branchNameToType_a[branch_name] == "vector<char>") { // Compare vector<char>
        if (helper_a.v_char(branch_name).GetSize() != helper_b.v_char(branch_name).GetSize()) { // Compare size
          if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<" has different sizes: "<<helper_a.v_char(branch_name).GetSize()<<" "<<helper_b.v_char(branch_name).GetSize()<<endl;
          different_branches[branch_name].push_back({entry_a, -2});
        } else { // Compare values
          for (unsigned iValue=0; iValue < helper_a.v_char(branch_name).GetSize(); ++iValue) {
            if (!std::isnan(helper_a.v_char(branch_name)[iValue]) || !std::isnan(helper_b.v_char(branch_name)[iValue])) {
             if (helper_a.v_char(branch_name)[iValue] != helper_b.v_char(branch_name)[iValue]) {
               if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<"["<<iValue<<"] has different values: "<<helper_a.v_char(branch_name)[iValue]<<" "<<helper_b.v_char(branch_name)[iValue]<<endl;
               different_branches[branch_name].push_back({entry_a, iValue});
             }
            }
          }
        }
      } else if (branchNameToType_a[branch_name] == "vector<int>") { // Compare vector<int>
        if (helper_a.v_int(branch_name).GetSize() != helper_b.v_int(branch_name).GetSize()) { // Compare size
          if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<" has different sizes: "<<helper_a.v_int(branch_name).GetSize()<<" "<<helper_b.v_int(branch_name).GetSize()<<endl;
          different_branches[branch_name].push_back({entry_a, -2});
        } else { // Compare values
          for (unsigned iValue=0; iValue < helper_a.v_int(branch_name).GetSize(); ++iValue) {
            if (!std::isnan(helper_a.v_int(branch_name)[iValue]) || !std::isnan(helper_b.v_int(branch_name)[iValue])) {
             if (helper_a.v_int(branch_name)[iValue] != helper_b.v_int(branch_name)[iValue]) {
               if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<"["<<iValue<<"] has different values: "<<helper_a.v_int(branch_name)[iValue]<<" "<<helper_b.v_int(branch_name)[iValue]<<endl;
               different_branches[branch_name].push_back({entry_a, iValue});
             }
            }
          }
        }
      } else if (branchNameToType_a[branch_name] == "vector<float>") { // Compare vector<float>
        if (helper_a.v_float(branch_name).GetSize() != helper_b.v_float(branch_name).GetSize()) { // Compare size
          if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<" has different sizes: "<<helper_a.v_float(branch_name).GetSize()<<" "<<helper_b.v_float(branch_name).GetSize()<<endl;
          different_branches[branch_name].push_back({entry_a, -2});
        } else { // Compare values
          for (unsigned iValue=0; iValue < helper_a.v_float(branch_name).GetSize(); ++iValue) {
            if (!std::isnan(helper_a.v_float(branch_name)[iValue]) || !std::isnan(helper_b.v_float(branch_name)[iValue])) {
              if (helper_a.v_float(branch_name)[iValue] != helper_b.v_float(branch_name)[iValue]) {
                if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<"["<<iValue<<"] has different values: "<<helper_a.v_float(branch_name)[iValue]<<" "<<helper_b.v_float(branch_name)[iValue]<<endl;
                different_branches[branch_name].push_back({entry_a, iValue});
              }
            }
          }
        }
      } else if (branchNameToType_a[branch_name] == "vector<double>") { // Compare vector<double>
        if (helper_a.v_double(branch_name).GetSize() != helper_b.v_double(branch_name).GetSize()) { // Compare size
          if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<" has different sizes: "<<helper_a.v_double(branch_name).GetSize()<<" "<<helper_b.v_double(branch_name).GetSize()<<endl;
          different_branches[branch_name].push_back({entry_a, -2});
        } else { // Compare values
          for (unsigned iValue=0; iValue < helper_a.v_double(branch_name).GetSize(); ++iValue) {
            if (!std::isnan(helper_a.v_double(branch_name)[iValue]) || !std::isnan(helper_b.v_double(branch_name)[iValue])) {
              if (helper_a.v_double(branch_name)[iValue] != helper_b.v_double(branch_name)[iValue]) {
                if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<"["<<iValue<<"] has different values: "<<helper_a.v_double(branch_name)[iValue]<<" "<<helper_b.v_double(branch_name)[iValue]<<endl;
                different_branches[branch_name].push_back({entry_a, iValue});
              }
            }
          }
        }
      } else if (branchNameToType_a[branch_name] == "vector<TLorentzVector>") { // Compare vector<TLorentzVector>
        if (helper_a.v_TLorentzVector(branch_name).GetSize() != helper_b.v_TLorentzVector(branch_name).GetSize()) { // Compare size
          if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<" has different sizes: "<<helper_a.v_TLorentzVector(branch_name).GetSize()<<" "<<helper_b.v_TLorentzVector(branch_name).GetSize()<<endl;
          different_branches[branch_name].push_back({entry_a, -2});
        } else { // Compare values
          for (unsigned iValue=0; iValue < helper_a.v_TLorentzVector(branch_name).GetSize(); ++iValue) {
            if (helper_a.v_TLorentzVector(branch_name)[iValue] != helper_b.v_TLorentzVector(branch_name)[iValue]) {
              if (different_branches.find(branch_name) == different_branches.end()) cout<<"[Info] "<<branch_name<<"["<<iValue<<"] has different values. "<<helper_a.v_TLorentzVector(branch_name)[iValue].Pt()<<" "<<helper_b.v_TLorentzVector(branch_name)[iValue].Pt()<<endl;
              different_branches[branch_name].push_back({entry_a, iValue});
            }
          }
        }
      } else {
        cout<<"[Info] Unknown branch type ("<<branchNameToType_a[branch_name]<<") for "<<branch_name<<endl;
      }
    } // Loop over branches
    //break;
  } // Compare for each event

  // Print results of different events
  if (different_branches.size() != 0) {
    cout<<"[Info] Differences in following files."<<endl;
    cout<<"  - "<<root_filename_a<<endl;
    cout<<"  - "<<root_filename_b<<endl;
    cout<<"(entry, iValue), where iValue=-1 means the branch is a scalar and -2 means the vector size is different."<<endl;
    for (auto mBranch : different_branches) {
      cout<<"Branch :"<<mBranch.first<<" is different ("<<mBranch.second.size()<<" times): ";
      int iEntry = 0;
      for (auto entry_value: mBranch.second) {
        cout<<"("<<get<0>(entry_value)<<" "<<get<1>(entry_value)<<"), ";
        iEntry++;
        if (iEntry==50) break;
      }
      cout<<endl;
    }
  } else {
    cout<<"[Info] There were no differences in common branches."<<endl;
  }

  time(&endtime);
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
