#ifndef JTREEREADERHELPER_H
#define JTREEREADERHELPER_H
#include <iostream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

using std::cout;
using std::endl;

class JTreeReaderHelper {
  public:
    ~JTreeReaderHelper() {
      clearStorage();

      for (auto it = m_variable.begin(); it != m_variable.end(); ++it) {
        string const & branchName = it->first;
        string const & branchType = m_variableNameToType[branchName];
        if (branchType == "ULong64_t") {
          if (it->second) delete static_cast<ULong64_t *> (it->second);
        } else if (branchType == "Long64_t") {
          if (it->second) delete static_cast<Long64_t *> (it->second);
        } else if (branchType == "Int_t") {
          if (it->second) delete static_cast<int *> (it->second);
        } else if (branchType == "Double_t") {
          if (it->second) delete static_cast<double *> (it->second);

        } else if (branchType == "vector<char>") {
          if (it->second) delete static_cast<vector<char> *> (it->second);
        } else if (branchType == "vector<int>") {
          if (it->second) delete static_cast<vector<int> *> (it->second);
        } else if (branchType == "vector<double>") {
          if (it->second) delete static_cast<vector<double> *> (it->second);
        } else if (branchType == "vector<TLorentzVector>") {
          if (it->second) delete static_cast<vector<TLorentzVector> *> (it->second);
        }
      }

    }

    void clearStorage() {
      for (auto it = m_empty.begin(); it != m_empty.end(); ++it) {
        if (it->second) {
          delete static_cast<bool*>(it->second);
          it->second = 0;
        } 
      }
      m_empty.clear();

      for (auto it = m_storage.begin(); it != m_storage.end(); ++it) {
        string const & branchName = it->first;
        string const & branchType = m_branchNameToType[branchName];
        if (branchType == "Bool_t") {
          if (it->second) {
            delete static_cast<TTreeReaderValue<bool> *> (it->second);
            it->second = 0;
          }
        } 
        else if (branchType == "UInt_t") {
          if (it->second) {
            delete static_cast<TTreeReaderValue<UInt_t> *> (it->second);
            it->second = 0;
          }
        } else if (branchType == "ULong64_t") {
          if (it->second) {
            delete static_cast<TTreeReaderValue<ULong64_t> *> (it->second);
            it->second = 0;
          }
        } else if (branchType == "Long64_t") {
          if (it->second) {
            delete static_cast<TTreeReaderValue<Long64_t> *> (it->second);
            it->second = 0;
          }
        } else if (branchType == "Int_t") {
          if (it->second) {
            delete static_cast<TTreeReaderValue<int> *> (it->second);
            it->second = 0;
          }
        } else if (branchType == "Float_t") {
          if (it->second) {
            delete static_cast<TTreeReaderValue<float> *> (it->second);
            it->second = 0;
          }
        } else if (branchType == "Double_t") {
          if (it->second) {
            delete static_cast<TTreeReaderValue<double> *> (it->second);
            it->second = 0;
          }

        } else if (branchType == "vector<bool>") {
          if (it->second) {
            delete static_cast<TTreeReaderValue<vector<bool>> *> (it->second);
            it->second = 0;
          }
        } else if (branchType == "vector<char>") {
          if (it->second) {
            delete static_cast<TTreeReaderValue<vector<char>> *> (it->second);
            it->second = 0;
          }
        } else if (branchType == "vector<int>") {
          if (it->second) {
            delete static_cast<TTreeReaderValue<vector<int>> *> (it->second);
            it->second = 0;
          }
        } else if (branchType == "vector<float>") {
          if (it->second) {
            delete static_cast<TTreeReaderValue<vector<float>> *> (it->second);
            it->second = 0;
          }
        } else if (branchType == "vector<double>") {
          if (it->second) {
            delete static_cast<TTreeReaderValue<vector<double>> *> (it->second);
            it->second = 0;
          }
        } else if (branchType == "vector<TLorentzVector>") {
          if (it->second) {
            delete static_cast<TTreeReaderValue<vector<TLorentzVector>> *> (it->second);
            it->second = 0;
          }
        }
      }
      m_storage.clear();
    }

    // Loop over all trees to get all branchNameToType
    static map<string, string> getBranchNameToType(TChain * chain) {
      map<string, string> branchNameToType;
      // Loop over all trees
      int nTree = chain->GetNtrees();
      Long64_t currentEntry = 0;
      for (int iTree = 0; iTree < nTree; ++iTree) {
        // Load tree
        chain->LoadTree(currentEntry);
        TTree * currentTree = chain->GetTree(); // tree address might not change.
        //cout<<"Current File: "<<currentTree->GetCurrentFile()->GetName()<<endl;

        map<string, string> treeBranchNameToType = getBranchNameToType(currentTree);
        branchNameToType.insert(treeBranchNameToType.begin(), treeBranchNameToType.end());

        // Go to next tree
        currentEntry += currentTree->GetEntries();
      }

      return branchNameToType;
    }

    static map<string, string> getBranchNameToType(TTree * tree) {
      map<string, string> branchNameToType;

      TObjArray * leaves = tree->GetListOfLeaves();
      TLeaf * leaf=0; string leafType; string leafName; int leafLen;
      for (int iLeaf = 0; iLeaf < leaves->GetSize(); ++iLeaf) {
        leaf = static_cast<TLeaf*>(leaves->At(iLeaf));
        leafType = leaf->GetTypeName();
        leafName = leaf->GetName();
        leafLen = leaf->GetLen();
        if (leafLen == 1) branchNameToType[leafName] = leafType;
        else {
          //cout<<"leaf: "<<leafName<<" "<<leafType<<" "<<leafLen<<endl;
          branchNameToType[leafName] = "vector<"+leafType+">";
        }
      }

      // Loop over all files to find all bool leaf types. 
      if (string(tree->ClassName()) == "TChain") {
        TChain * chain = static_cast<TChain*> (tree);
        TObjArray * files = chain->GetListOfFiles();
        for (int iFile = 0; iFile < files->GetEntries(); ++iFile) {
          TChainElement *chainElement = static_cast<TChainElement*>(files->At(iFile));
          TFile subFile(chainElement->GetTitle());
          TTree * subTree = static_cast<TTree*>(subFile.Get("Events"));
          TObjArray * leaves = subTree->GetListOfLeaves();
          TLeaf * leaf=0; string leafType; string leafName; int leafLen;
          for (int iLeaf = 0; iLeaf < leaves->GetSize(); ++iLeaf) {
            leaf = static_cast<TLeaf*>(leaves->At(iLeaf));
            leafName = leaf->GetName();
            leafType = leaf->GetTypeName();
            leafLen = leaf->GetLen();
            if (leafLen == 1 && leafType == "Bool_t") branchNameToType[leafName] = leafType;
          }
        }
      }

      return branchNameToType;
    }

    void setBranches(TTreeReader & treeReader, map<string, string> branchNameToType) {
      for (auto it = branchNameToType.begin(); it != branchNameToType.end(); ++it) {
        string const & branchName = it->first;
        string const & branchType = it->second;
        m_branchNameToType[branchName] = branchType; // Used to delete memory
        //cout<<"iTree: "<<iTree<<" branchName: "<<branchName<<" "<<readerTree->FindLeaf(branchName.c_str())<<endl;
        if (!(treeReader.GetTree()->FindLeaf(branchName.c_str()))) { // Check if tree has branch
          m_empty[branchName] = new bool(0);
        } else if (branchType == "Bool_t") {
          m_storage[branchName] = new TTreeReaderValue<bool> (treeReader, branchName.c_str());
        } else if (branchType == "UInt_t") {
          m_storage[branchName] = new TTreeReaderValue<UInt_t> (treeReader, branchName.c_str());
        } else if (branchType == "ULong64_t") {
          m_storage[branchName] = new TTreeReaderValue<ULong64_t> (treeReader, branchName.c_str());
        } else if (branchType == "Long64_t") {
          m_storage[branchName] = new TTreeReaderValue<Long64_t> (treeReader, branchName.c_str());
        } else if (branchType == "Int_t") {
          m_storage[branchName] = new TTreeReaderValue<int> (treeReader, branchName.c_str());
        } else if (branchType == "Float_t") {
          m_storage[branchName] = new TTreeReaderValue<float> (treeReader, branchName.c_str());
        } else if (branchType == "Double_t") {
          m_storage[branchName] = new TTreeReaderValue<double> (treeReader, branchName.c_str());

        } else if (branchType == "vector<bool>" || branchType == "vector<Bool_t>") {
          m_storage[branchName] = new TTreeReaderArray<bool> (treeReader, branchName.c_str());
        } else if (branchType == "vector<char>" || branchType == "vector<Char_t>") {
          m_storage[branchName] = new TTreeReaderArray<char> (treeReader, branchName.c_str());
        } else if (branchType == "vector<int>" || branchType == "vector<Int_t>") {
          m_storage[branchName] = new TTreeReaderArray<int> (treeReader, branchName.c_str());
        } else if (branchType == "vector<float>" || branchType == "vector<Float_t>") {
          m_storage[branchName] = new TTreeReaderArray<float> (treeReader, branchName.c_str());
        } else if (branchType == "vector<double>" || branchType == "vector<Double_t>") {
          m_storage[branchName] = new TTreeReaderArray<double> (treeReader, branchName.c_str());
        } else if (branchType == "vector<TLorentzVector>" || branchType == "vector<TLorentzVector>") {
          m_storage[branchName] = new TTreeReaderArray<TLorentzVector> (treeReader, branchName.c_str());
        } else {
          //cout<<"[Info] JTreeReaderHelper::setBranches(): Ignore unknown type: "<<branchType<<" "<<branchName<<endl;
        }
      }
    }

    bool const & s_bool(string const & member) const {
      #ifdef JDEBUG
      if (m_branchNameToType.find(member) == m_branchNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::s_bool(): No branch called "<<member<<endl;
        throw;
      }
      if (m_branchNameToType.at(member) != "Bool_t") {
        cout<<"[Error] JTreeReaderHelper::s_bool(): "<<member<<" type is "<<m_branchNameToType.at(member)<<endl;
        throw;
      }
      #endif
      // To handle case branch does not exist. Specific to handle FastSim GenModel branches
      if (m_empty.count(member)) {
        return *(static_cast<bool*>(m_empty.at(member)));
      }
      return **(static_cast<TTreeReaderValue<bool> *> (m_storage.at(member)));
    }
    UInt_t const & s_uint(string const & member) const {
      #ifdef JDEBUG
      if (m_branchNameToType.find(member) == m_branchNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::s_uint(): No branch called "<<member<<endl;
        throw;
      }
      if (m_branchNameToType.at(member) != "UInt_t") {
        cout<<"[Error] JTreeReaderHelper::s_uint(): "<<member<<" type is "<<m_branchNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return **(static_cast<TTreeReaderValue<UInt_t> *> (m_storage.at(member)));
    }
    ULong64_t const & s_ulong64(string const & member) const {
      #ifdef JDEBUG
      if (m_branchNameToType.find(member) == m_branchNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::s_ulong64(): No branch called "<<member<<endl;
        throw;
      }
      if (m_branchNameToType.at(member) != "ULong64_t") {
        cout<<"[Error] JTreeReaderHelper::s_ulong64(): "<<member<<" type is "<<m_branchNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return **(static_cast<TTreeReaderValue<ULong64_t> *> (m_storage.at(member)));
    }
    Long64_t const & s_long64(string const & member) const {
      #ifdef JDEBUG
      if (m_branchNameToType.find(member) == m_branchNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::s_long64(): No branch called "<<member<<endl;
        throw;
      }
      if (m_branchNameToType.at(member) != "Long64_t") {
        cout<<"[Error] JTreeReaderHelper::s_long64(): "<<member<<" type is "<<m_branchNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return **(static_cast<TTreeReaderValue<Long64_t> *> (m_storage.at(member)));
    }
    int const & s_int(string const & member) const {
      #ifdef JDEBUG
      if (m_branchNameToType.find(member) == m_branchNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::s_int(): No branch called "<<member<<endl;
        throw;
      }
      if (m_branchNameToType.at(member) != "Int_t") {
        cout<<"[Error] JTreeReaderHelper::s_int(): "<<member<<" type is "<<m_branchNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return **(static_cast<TTreeReaderValue<int> *> (m_storage.at(member)));
    }
    float const & s_float(string const & member) const {
      #ifdef JDEBUG
      if (m_branchNameToType.find(member) == m_branchNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::s_float(): No branch called "<<member<<endl;
        throw;
      }
      if (m_branchNameToType.at(member) != "Float_t") {
        cout<<"[Error] JTreeReaderHelper::s_float(): "<<member<<" type is "<<m_branchNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return **(static_cast<TTreeReaderValue<float> *> (m_storage.at(member)));
    }
    double const & s_double(string const & member) const {
      #ifdef JDEBUG
      if (m_branchNameToType.find(member) == m_branchNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::s_double(): No branch called "<<member<<endl;
        throw;
      }
      if (m_branchNameToType.at(member) != "Double_t") {
        cout<<"[Error] JTreeReaderHelper::s_double() "<<member<<" type is "<<m_branchNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return **(static_cast<TTreeReaderValue<double> *> (m_storage.at(member)));
    }

    TTreeReaderArray<bool> const & v_bool(string const & member) const {
      #ifdef JDEBUG
      if (m_branchNameToType.find(member) == m_branchNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::v_bool(): No branch called "<<member<<endl;
        throw;
      }
      if (m_branchNameToType.at(member) != "vector<Bool_t>" && m_branchNameToType.at(member) != "vector<bool>") {
        cout<<"[Error] JTreeReaderHelper::v_bool(): "<<member<<" type is "<<m_branchNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return *(static_cast<TTreeReaderArray<bool> *> (m_storage.at(member)));
    }
    TTreeReaderArray<char> const & v_char(string const & member) const {
      #ifdef JDEBUG
      if (m_branchNameToType.find(member) == m_branchNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::v_char(): No branch called "<<member<<endl;
        throw;
      }
      if (m_branchNameToType.at(member) != "vector<Char_t>" && m_branchNameToType.at(member) != "vector<char>") {
        cout<<"[Error] JTreeReaderHelper::v_char(): "<<member<<" type is "<<m_branchNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return *(static_cast<TTreeReaderArray<char> *> (m_storage.at(member)));
    }
    TTreeReaderArray<int> const & v_int(string const & member) const {
      #ifdef JDEBUG
      if (m_branchNameToType.find(member) == m_branchNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::v_int(): No branch called "<<member<<endl;
        throw;
      }
      if (m_branchNameToType.at(member) != "vector<Int_t>" && m_branchNameToType.at(member) != "vector<int>") {
        cout<<"[Error] JTreeReaderHelper::v_int(): "<<member<<" type is "<<m_branchNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return *(static_cast<TTreeReaderArray<int> *> (m_storage.at(member)));
    }
    TTreeReaderArray<float> const & v_float(string const & member) const {
      #ifdef JDEBUG
      if (m_branchNameToType.find(member) == m_branchNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::v_float(): No branch called "<<member<<endl;
        throw;
      }
      if (m_branchNameToType.at(member) != "vector<Float_t>" && m_branchNameToType.at(member) != "vector<float>") {
        cout<<"[Error] JTreeReaderHelper::v_float(): "<<member<<" type is "<<m_branchNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return *(static_cast<TTreeReaderArray<float> *> (m_storage.at(member)));
    }
    TTreeReaderArray<double> const & v_double(string const & member) const {
      #ifdef JDEBUG
      if (m_branchNameToType.find(member) == m_branchNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::v_double(): No branch called "<<member<<endl;
        throw;
      }
      if (m_branchNameToType.at(member) != "vector<Double_t>" && m_branchNameToType.at(member) != "vector<double>") {
        cout<<"[Error] JTreeReaderHelper::v_double(): "<<member<<" type is "<<m_branchNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return *(static_cast<TTreeReaderArray<double> *> (m_storage.at(member)));
    }
    TTreeReaderArray<TLorentzVector> const & v_TLorentzVector(string const & member) const {
      #ifdef JDEBUG
      if (m_branchNameToType.find(member) == m_branchNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::v_TLorentzVector(): No branch called "<<member<<endl;
        throw;
      }
      if (m_branchNameToType.at(member) != "vector<TLorentzVector>") {
        cout<<"[Error] JTreeReaderHelper::v_TLorentzVector(): "<<member<<" type is "<<m_branchNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return *(static_cast<TTreeReaderArray<TLorentzVector> *> (m_storage.at(member)));
    }


    // Set variable
    void usr_s_int(string const & member, int const & value) {
      // Check if already defined
      if (m_variableNameToType.find(member) == m_variableNameToType.end()) {
        m_variable[member] = new int;
        m_variableNameToType[member] = "Int_t";
      }
      // Set variable
      (*static_cast<int *>(m_variable[member])) = value;
    }
    void usr_s_ulong64(string const & member, ULong64_t const & value) {
      // Check if already defined
      if (m_variableNameToType.find(member) == m_variableNameToType.end()) {
        m_variable[member] = new ULong64_t;
        m_variableNameToType[member] = "ULong64_t";
      }
      // Set variable
      (*static_cast<ULong64_t *>(m_variable[member])) = value;
    }
    void usr_s_long64(string const & member, Long64_t const & value) {
      // Check if already defined
      if (m_variableNameToType.find(member) == m_variableNameToType.end()) {
        m_variable[member] = new Long64_t;
        m_variableNameToType[member] = "Long64_t";
      }
      // Set variable
      (*static_cast<Long64_t *>(m_variable[member])) = value;
    }
    void usr_s_double(string const & member, double const & value) {
      // Check if already defined
      if (m_variableNameToType.find(member) == m_variableNameToType.end()) {
        m_variable[member] = new double;
        m_variableNameToType[member] = "Double_t";
      }
      // Set variable
      (*static_cast<double *>(m_variable[member])) = value;
    }
    void usr_v_char(string const & member, vector<char> const & value) {
      // Check if already defined
      if (m_variableNameToType.find(member) == m_variableNameToType.end()) {
        m_variable[member] = new vector<char>;
        m_variableNameToType[member] = "vector<char>";
      }
      // Set variable
      (*static_cast<vector<char> *>(m_variable[member])) = value;
    }
    void usr_v_int(string const & member, vector<int> const & value) {
      // Check if already defined
      if (m_variableNameToType.find(member) == m_variableNameToType.end()) {
        m_variable[member] = new vector<int>;
        m_variableNameToType[member] = "vector<int>";
      }
      // Set variable
      (*static_cast<vector<int> *>(m_variable[member])) = value;
    }
    void usr_v_double(string const & member, vector<double> const & value) {
      // Check if already defined
      if (m_variableNameToType.find(member) == m_variableNameToType.end()) {
        m_variable[member] = new vector<double>;
        m_variableNameToType[member] = "vector<double>";
      }
      // Set variable
      (*static_cast<vector<double> *>(m_variable[member])) = value;
    }
    void usr_v_TLorentzVector(string const & member, vector<TLorentzVector> const & value) {
      // Check if already defined
      if (m_variableNameToType.find(member) == m_variableNameToType.end()) {
        m_variable[member] = new vector<TLorentzVector>;
        m_variableNameToType[member] = "vector<TLorentzVector>";
      }
      // Set variable
      (*static_cast<vector<TLorentzVector> *>(m_variable[member])) = value;
    }

    // Get variable
    int const & usr_s_int(string const & member) const {
      #ifdef JDEBUG
      if (m_variableNameToType.find(member) == m_variableNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::usr_s_int(): No branch called "<<member<<endl;
        throw;
      }
      if (m_variableNameToType.at(member) != "Int_t") {
        cout<<"[Error] JTreeReaderHelper::usr_s_int(): "<<member<<" type is "<<m_variableNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return *(static_cast<int *> (m_variable.at(member)));
    }
    ULong64_t const & usr_s_ulong64(string const & member) const {
      #ifdef JDEBUG
      if (m_variableNameToType.find(member) == m_variableNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::usr_s_ulong64(): No branch called "<<member<<endl;
        throw;
      }
      if (m_variableNameToType.at(member) != "ULong64_t") {
        cout<<"[Error] JTreeReaderHelper::usr_s_ulong64(): "<<member<<" type is "<<m_variableNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return *(static_cast<ULong64_t *> (m_variable.at(member)));
    }
    Long64_t const & usr_s_long64(string const & member) const {
      #ifdef JDEBUG
      if (m_variableNameToType.find(member) == m_variableNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::usr_s_long64(): No branch called "<<member<<endl;
        throw;
      }
      if (m_variableNameToType.at(member) != "Long64_t") {
        cout<<"[Error] JTreeReaderHelper::usr_s_long64(): "<<member<<" type is "<<m_variableNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return *(static_cast<Long64_t *> (m_variable.at(member)));
    }
    double const & usr_s_double(string const & member) const {
      #ifdef JDEBUG
      if (m_variableNameToType.find(member) == m_variableNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::usr_s_double(): No branch called "<<member<<endl;
        throw;
      }
      if (m_variableNameToType.at(member) != "Double_t") {
        cout<<"[Error] JTreeReaderHelper::usr_s_double(): "<<member<<" type is "<<m_variableNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return *(static_cast<double *> (m_variable.at(member)));
    }
    vector<char> const & usr_v_char(string const & member) const {
      #ifdef JDEBUG
      if (m_variableNameToType.find(member) == m_variableNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::usr_v_char(): No branch called "<<member<<endl;
        throw;
      }
      if (m_variableNameToType.at(member) != "vector<int>") {
        cout<<"[Error] JTreeReaderHelper::usr_v_char(): "<<member<<" type is "<<m_variableNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return *(static_cast<vector<char> *> (m_variable.at(member)));
    }
    vector<int> const & usr_v_int(string const & member) const {
      #ifdef JDEBUG
      if (m_variableNameToType.find(member) == m_variableNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::usr_v_int(): No branch called "<<member<<endl;
        throw;
      }
      if (m_variableNameToType.at(member) != "vector<int>") {
        cout<<"[Error] JTreeReaderHelper::usr_v_int(): "<<member<<" type is "<<m_variableNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return *(static_cast<vector<int> *> (m_variable.at(member)));
    }
    vector<double> const & usr_v_double(string const & member) const {
      #ifdef JDEBUG
      if (m_variableNameToType.find(member) == m_variableNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::usr_v_double(): No branch called "<<member<<endl;
        throw;
      }
      if (m_variableNameToType.at(member) != "vector<double>") {
        cout<<"[Error] JTreeReaderHelper::usr_v_double(): "<<member<<" type is "<<m_variableNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return *(static_cast<vector<double> *> (m_variable.at(member)));
    }
    vector<TLorentzVector> const & usr_v_TLorentzVector(string const & member) const {
      #ifdef JDEBUG
      if (m_variableNameToType.find(member) == m_variableNameToType.end()) {
        cout<<"[Error] JTreeReaderHelper::usr_v_TLorentzVector(): No branch called "<<member<<endl;
        throw;
      }
      if (m_variableNameToType.at(member) != "vector<TLorentzVector>") {
        cout<<"[Error] JTreeReaderHelper::usr_v_TLorentzVector(): "<<member<<" type is "<<m_variableNameToType.at(member)<<endl;
        throw;
      }
      #endif
      return *(static_cast<vector<TLorentzVector> *> (m_variable.at(member)));
    }

  private:
    // To handle when a branch does not exist. Specific to handle FastSim GenModel branches
    map<string, void * > m_empty;
    map<string, void * > m_storage;
    map<string, string> m_branchNameToType; // Used for deleting memory

    map<string, void *> m_variable; // Values defined in process
    map<string, string> m_variableNameToType;
};

class JTreeReader {
  public:
    JTreeReader(TChain * chain) {
      m_chain = chain;
      m_currentEntry = -1;
      m_currentTreeIndex = -1;
      m_currentReader = 0;
      m_helper = 0;
      m_currentReader = 0;

      // Fill m_firstEntryOfTrees
      int nTree = m_chain->GetNtrees();
      m_firstEntryOfTrees.reserve(nTree);
      Long64_t currentEntry = 0;
      // Loop over trees
      for (int iTree = 0; iTree < nTree; ++iTree) {
        // Save iTree entry
        m_firstEntryOfTrees.push_back(currentEntry);

        // Load iTree+1
        m_chain->LoadTree(currentEntry);
        TTree * currentTree = m_chain->GetTree(); // tree address might not change.
        //cout<<"Current File: "<<currentTree->GetCurrentFile()->GetName()<<endl;

        // Go to iTree+1 first entry
        currentEntry += currentTree->GetEntries();
      }
    };
    ~JTreeReader() {
      if (!m_currentReader) delete m_currentReader;
    }
    TTreeReader * getTreeReader() {
      return m_currentReader;
    }
    Long64_t getEntryNumber() {
      return m_currentEntry;
    }

    bool Next() {
      // m_currentEntry = 0, m_currentTreeIndex = 0, TTReader = -1
      m_currentEntry++;
      // Switch tree and make reader according to current entry
      if (m_currentEntry == 0) {
        m_currentTreeIndex = 0;

        m_chain->LoadTree(m_currentEntry);
        m_currentReader = new TTreeReader(m_chain->GetTree());
        m_helper->setBranches(*m_currentReader, m_branchNameToType);
        //cout<<"Filename: "<<m_chain->GetTree()->GetCurrentFile()->GetName()<<" entries: "<<m_chain->GetTree()->GetEntries()<<endl;
        cout<<"Processing ["<<m_currentTreeIndex+1<<"/"<<m_firstEntryOfTrees.size()<<"]: "<<m_chain->GetTree()->GetCurrentFile()->GetName()<<endl;
      } else if (m_currentEntry == m_firstEntryOfTrees[m_currentTreeIndex+1]) {
        m_currentTreeIndex++;

        delete m_currentReader;
        m_helper->clearStorage();

        m_chain->LoadTree(m_currentEntry);
        m_currentReader = new TTreeReader(m_chain->GetTree());
        m_helper->setBranches(*m_currentReader, m_branchNameToType);
        //cout<<"Filename: "<<m_chain->GetTree()->GetCurrentFile()->GetName()<<" entries: "<<m_chain->GetTree()->GetEntries()<<endl;
        cout<<"Processing ["<<m_currentTreeIndex+1<<"/"<<m_firstEntryOfTrees.size()<<"]: "<<m_chain->GetTree()->GetCurrentFile()->GetName()<<endl;
      }

      return m_currentReader->Next();
    }

    void setBranches(map<string, string> branchNameToType, JTreeReaderHelper * helper) {
      m_branchNameToType = branchNameToType;
      m_helper = helper;
    }

  private:
    TChain * m_chain; // not owner
    // m_firstEntryOfTrees[tree] = firstEntry
    vector<Long64_t> m_firstEntryOfTrees;
    map<string, string> m_branchNameToType;
    JTreeReaderHelper * m_helper; // not owner
    TTreeReader * m_currentReader; // owner
    Long64_t m_currentEntry;
    int m_currentTreeIndex;
};

#endif
