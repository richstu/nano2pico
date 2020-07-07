#include "generate_tree_classes.hpp"

#include <cstring>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <set>

#include <unistd.h>

using namespace std;

string ToCaps(string str){
  for(string::iterator it = str.begin();
      it != str.end();
      ++it){
    *it = toupper(*it);
  }
  return str;
}

string execute(const string &cmd){
  FILE *pipe = popen(cmd.c_str(), "r");
  if(!pipe) throw runtime_error("Could not open pipe.");
  const size_t buffer_size = 128;
  char buffer[buffer_size];
  string result = "";
  while(!feof(pipe)){
    if(fgets(buffer, buffer_size, pipe) != NULL) result += buffer;
  }

  pclose(pipe);
  return result;
}

vector<string> Tokenize(const string& input, const string& tokens = " "){
  char* ipt(new char[input.size()+1]);
  memcpy(ipt, input.data(), input.size());
  ipt[input.size()]=static_cast<char>(0);
  char* ptr(strtok(ipt, tokens.c_str()));
  vector<string> output(0);
  while(ptr!=NULL){
    output.push_back(ptr);
    ptr=strtok(NULL, tokens.c_str());
  }
  return output;
}

string FixName(string name){
  //Variable names can have alphanumeric characters and underscores only
  string allowed = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890_";

  //Remove illegal characters
  size_t pos = name.find_first_not_of(allowed);
  while(pos < name.size()){
    name.at(pos) = '_';
    pos = name.find_first_not_of(allowed);
  }

  //Replace double underscore with single underscore
  pos = name.rfind("__");
  while(pos < name.size()){
    name.replace(pos, 2, "_");
    pos = name.rfind("__");
  }

  //Remove leading and trailing spaces
  pos = 0;
  for(pos = 0; pos < name.size(); ++pos){
    if(name.at(pos) != ' ') break;
  }
  size_t endpos = name.size();
  for(endpos = name.size(); endpos != 0; --endpos){
    if(name.at(endpos-1) != ' ') break;
  }

  return name.substr(pos, endpos-pos);
}

int GetArrayLength(const std::string var_name){
  const string var_type = var_name.substr(0,var_name.find_first_of('_'));
  if      (var_type == "Electron") return 40;
  else if (var_type == "Muon") return 40;
  else if (var_type == "Tau") return 30;
  else if (var_type == "IsoTrack") return 30;
  else if (var_type == "Photon") return 30;
  else if (var_type == "Jet") return 200;
  else if (var_type == "SoftActivityJet") return 30;
  else if (var_type == "FatJet") return 30;
  else if (var_type == "SubJet") return 60;
  else if (var_type == "CorrT1METJet") return 40;
  else if (var_type == "GenPart") return 500;
  else if (var_type == "GenDressedLepton") return 20;
  else if (var_type == "GenVisTau") return 20;
  else if (var_type == "GenJet") return 50;
  else if (var_type == "GenJetAK8") return 50;
  else if (var_type == "SubGenJetAK8") return 100;
  else if (var_type == "LHEPart") return 40;
  else if (var_type == "LHEPdfWeight") return 100;
  else if (var_type == "LHEReweightingWeight") return 100;
  else if (var_type == "LHEScaleWeight") return 10;
  else if (var_type == "PSWeight") return 20;
  else if (var_type == "OtherPV") return 20;
  else if (var_type == "SV") return 40;
  else if (var_type == "TrigObj") return 40;
  else return -1;
}

vector<Variable> GetVariables(const string &file_name){
  string allowed = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890_";
  vector<Variable> vars;

  ifstream infile(("txt/variables/"+file_name).c_str());
  string line;
  while(getline(infile, line)){
    size_t start = line.find_first_not_of(" ");
    if(start >= line.size() || line.at(start) == '#' || line.at(start) == '/') continue;

    //Replace double space with single space
    size_t pos = line.rfind("  ");
    while(pos < line.size()){
      line.replace(pos, 2, " ");
      pos = line.rfind("  ");
    }
    size_t end = line.find_last_of(allowed)+1;
    size_t split = line.rfind(' ', end)+1;
    string name_ = line.substr(split, end-split);
    string type_ = line.substr(start, split-start-1);
    if (line.find("array")!=string::npos) {
      type_ = type_.substr(0,type_.length()-1) + "," + to_string(GetArrayLength(name_))+">";
    }

    string base_type_ = "";
    if (line.find_first_of('<')!=string::npos) {
      size_t base_type_start = line.find_first_of('<')+1;
      size_t base_type_end = line.find_first_of('>');
      base_type_ = line.substr(base_type_start, base_type_end-base_type_start);
    } else {
      base_type_ = line.substr(start, split-start);
    }

    vars.push_back(Variable(type_, base_type_,name_));
  }
  infile.close();

  return vars;
}

int main(){

  // creating a read-only tree for Nano specifically because it has arrays...
  vector<Variable> nano_vars = GetVariables("nano");  
  WriteNanoHeader(nano_vars);
  WriteNanoSource(nano_vars);

  vector<Variable> baby_vars = GetVariables("baby");  
  WriteHeader(baby_vars, "baby");
  WriteSource(baby_vars, "baby");

  vector<Variable> pico_vars = GetVariables("pico");  
  WriteHeader(pico_vars, "pico");
  WriteSource(pico_vars, "pico");

  vector<Variable> atto_vars = GetVariables("atto");  
  WriteHeader(atto_vars, "atto");
  WriteSource(atto_vars, "atto");

  vector<Variable> higfeats_vars = GetVariables("higfeats");  
  WriteHeader(higfeats_vars, "higfeats");
  WriteSource(higfeats_vars, "higfeats");

  vector<Variable> zgfeats_vars = GetVariables("zgfeats");  
  WriteHeader(zgfeats_vars, "zgfeats");
  WriteSource(zgfeats_vars, "zgfeats");

  vector<Variable> corr_vars = GetVariables("corrections");
  WriteHeader(corr_vars, "corrections");
  WriteSource(corr_vars, "corrections");
}

bool Contains(const string &text, const string &pattern){
  return text.find(pattern) != string::npos;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// ---------------------------------------------------------------------------
//              NanoAOD reader, READ-ONLY
// ---------------------------------------------------------------------------
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

void WriteNanoHeader(const vector<Variable> &vars){

  ofstream file("inc/nano_tree.hpp");
  file << "// File generated with generate_tree_classes.exe\n\n";

  file << "#ifndef H_NANO_TREE\n";
  file << "#define H_NANO_TREE\n\n";

  file << "#include <vector>\n";
  file << "#include <string>\n";
  file << "#include <cmath>\n\n";
  file << "#include \"TTree.h\"\n";
  file << "#include \"TFile.h\"\n\n";
  file << "#include \"TChain.h\"\n\n";
  file << "#include \"TString.h\"\n\n";

  file << "class nano_tree{\n";
  file << "public:\n";
  file << "  nano_tree(TString infile = \"\"); \n\n";

  file << "  long GetEntries() const;\n";
  file << "  void GetEntry(const long entry);\n";
  file << "  double bad_val_;\n\n";
  file << "  ~nano_tree();\n\n";

  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector") || Contains(var->type_, "array")){
      file << "  std::vector<" << var->base_type_ << "> &" << var->name_ << "();\n";
    } else {
      file << "  " << var->type_ << " &" << var->name_ << "();\n";
    }
  }
  file << '\n';

  file << "  TChain* intree_;\n\n";

  file << "protected:\n";
  file << "private:\n";
  file << "  class VectorLoader{\n";
  file << "  public:\n";
  file << "    VectorLoader();\n";
  file << "  private:\n";
  file << "    static bool loaded_;\n";
  file << "  };\n\n";

  file << "  static VectorLoader vl_;\n";
  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "  std::" << var->type_ << ' ' << var->name_ << "_;\n";
      file << "  std::" << var->type_ << " *p_" << var->name_ << "_;\n";
    } else if(Contains(var->type_, "array")){
      file << "  std::" << var->type_ << " arr_" << var->name_ << "_;\n";
      file << "  std::vector<" << var->base_type_ << "> " << var->name_ << "_;\n";
    } else { 
      file << "  " << var->type_ << ' ' << var->name_ << "_;\n";
    }
    file << "  TBranch *b_" << var->name_ << "_;\n";
    file << "  mutable bool c_" << var->name_ << "_;\n";
  }

  file << "  long entry_;\n";

  file << "};\n\n";

  file << "#endif" << endl;

  file.close();
}

void WriteNanoSource(const vector<Variable> &vars){

  ofstream file("src/nano_tree.cpp");
  file << "//File generated with generate_tree_classes.exe\n\n";

  file << "#include \"nano_tree.hpp\"\n\n";

  file << "#include <stdexcept>\n";
  file << "#include <string>\n";
  file << "#include <iostream>\n";
  file << "#include <vector>\n\n";

  file << "#include \"TROOT.h\"\n";
  file << "#include \"TTree.h\"\n";
  file << "#include \"TBranch.h\"\n";
  file << "#include \"TChain.h\"\n";
  file << "#include \"TString.h\"\n";

  file << "using namespace std;\n\n";

  file << "bool nano_tree::VectorLoader::loaded_ = false;\n\n";

  file << "nano_tree::VectorLoader nano_tree::vl_ = nano_tree::VectorLoader();\n\n";

  file << "nano_tree::VectorLoader::VectorLoader(){\n";
  file << "  if(!loaded_){\n";
  file << "    gROOT->ProcessLine(\"#include <vector>\");\n";
  file << "    loaded_ = true;\n";
  file << "  }\n";
  file << "}\n\n";

  file << "nano_tree::nano_tree(TString infile):\n";
  file << "  bad_val_(-999.),\n";

  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "  " << var->name_ << "_(0),\n";
    } else if (Contains(var->type_, "array")) {
      file << "  arr_" << var->name_ << "_("<< var->type_ <<"{}),\n";
      file << "  " << var->name_ << "_(0),\n";
    } else if(Contains(var->type_, "tring")){
      file << "  " << var->name_ << "_(\"\"),\n";
    } else if(Contains(var->type_, "bool")){
      file << "  " << var->name_ << "_(false),\n";
    } else{
      file << "  " << var->name_ << "_(static_cast<" << var->type_ << ">(bad_val_)),\n";
    }
    if(Contains(var->type_, "vector")){
      file << "  p_" << var->name_ << "_(&" << var->name_ << "_),\n";
    }
    file << "  b_" << var->name_ << "_(NULL),\n";
    file << "  c_" << var->name_ << "_(false),\n";
  }
  file << "  entry_(0){\n";

  file << "  intree_ = new TChain(\"Events\");\n";
  file << "  intree_->Add(infile);\n\n";

  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "  intree_->SetBranchAddress(\"" << var->name_ << "\", &p_" << var->name_ << "_, &b_" << var->name_ << "_);\n";
    } else if (Contains(var->type_, "array")){
      file << "  intree_->SetBranchAddress(\"" << var->name_ << "\", &arr_" << var->name_ << "_[0], &b_" << var->name_ << "_);\n";
    } else {
      file << "  intree_->SetBranchAddress(\"" << var->name_ << "\", &" << var->name_ << "_, &b_" << var->name_ << "_);\n";
    }
  }
  file << "}\n\n"; // closing constructor bracket

  file << "nano_tree::~nano_tree(){\n";
  file << "}\n\n";
  
  file << "long nano_tree::GetEntries() const{\n";
  file << "  return intree_->GetEntries();\n";
  file << "}\n\n";

  file << "void nano_tree::GetEntry(const long entry){\n";
  file << "  // Reset read-trackers for input branches \n";
  for(vector<Variable>::const_iterator var = vars.begin(); var!= vars.end(); ++var){
    file << "  c_" << var->name_ << "_ = false;\n";
  }
  file << "  entry_ = intree_->LoadTree(entry);\n";
  file << "}\n\n";

  // Writing methods called when accessing input tree variable
  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if (Contains(var->type_, "array")){
      file << "vector<"<< var->base_type_ << "> & nano_tree::" << var->name_ << "(){\n";
    } else {
      file << var->type_ << " & nano_tree::" << var->name_ << "(){\n";
    }
    file << "  if(!c_" << var->name_ << "_ && b_" << var->name_ <<"_){\n";
    if (Contains(var->type_, "array")){
      file << "    int bytes = b_" << var->name_ << "_->GetEntry(entry_);\n";
      file << "    "<< var->name_ <<"_ = vector<"<< var->base_type_ <<">(arr_"<< var->name_ <<"_.begin(), arr_"<< var->name_ <<"_.begin()+bytes/sizeof(arr_"<< var->name_ <<"_[0]));\n";
    } else {
      file << "    b_" << var->name_ << "_->GetEntry(entry_);\n";
    }
    file << "    c_" << var->name_ << "_ = true;\n";
    file << "  }\n";
    file << "  return " << var->name_ << "_;\n";
    file << "}\n\n";
  }

  file.close();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// ---------------------------------------------------------------------------
//              STANDARD USCB TREE READER/WRITER, i.e. with vectors
// ---------------------------------------------------------------------------
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

void WriteHeader(const vector<Variable> &vars, const string name){

  ofstream file("inc/"+ name +"_tree.hpp");
  file << "// File generated with generate_tree_classes.exe\n\n";

  file << "#ifndef H_"<< name <<"_tree\n";
  file << "#define H_"<< name <<"_tree\n\n";

  file << "#include <vector>\n";
  file << "#include <string>\n";
  file << "#include <cmath>\n\n";
  file << "#include \"TTree.h\"\n";
  file << "#include \"TFile.h\"\n\n";
  file << "#include \"TChain.h\"\n\n";
  file << "#include \"TString.h\"\n\n";

  file << "class "<< name <<"_tree{\n";
  file << "public:\n";
  file << "  "<< name <<"_tree(TString infile = \"\", TString outfile = \"\"); \n\n";

  file << "  long GetEntries() const;\n";
  file << "  void GetEntry(const long entry);\n";

  file << "  void Fill();\n";
  file << "  void Write();\n\n";

  file << "  bool readOnly_;\n\n";
  file << "  bool writeOnly_;\n\n";
  file << "  bool copyUntouched_;\n\n";
  file << "  double bad_val_;\n\n";

  file << "  ~"<< name <<"_tree();\n\n";

  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "  std::" << var->type_ << " " << var->name_ << "();\n";
    } else {
      file << "  " << var->type_ << " " << var->name_ << "();\n";
    }
  }
  file << '\n';

  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "  std::" << var->type_ << "& out_" << var->name_ << "();\n";
    } else {
      file << "  " << var->type_ << "& out_" << var->name_ << "();\n";
    }
  }
  file << '\n';

  file << "  TFile* outfile_;\n";
  file << "  TChain* intree_;\n";
  file << "  TTree* outtree_;\n";
  file << '\n';

  file << "protected:\n";

  file << "private:\n";
  file << "  class VectorLoader{\n";
  file << "  public:\n";
  file << "    VectorLoader();\n";
  file << "  private:\n";
  file << "    static bool loaded_;\n";
  file << "  };\n\n";

  file << "  static VectorLoader vl_;\n";
  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "  std::" << var->type_ << ' ' << var->name_ << "_;\n";
      file << "  std::" << var->type_ << " *p_" << var->name_ << "_;\n";
    } else { 
      file << "  " << var->type_ << ' ' << var->name_ << "_;\n";
    }
    file << "  TBranch *b_" << var->name_ << "_;\n";
    file << "  mutable bool c_" << var->name_ << "_;\n";
  }

  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "  std::" << var->type_ << " out_" << var->name_ << "_;\n";
      file << "  std::" << var->type_ << " *p_out_" << var->name_ << "_;\n";
    } else { 
      file << "  " << var->type_ << " out_" << var->name_ << "_;\n";
    }
    file << "  mutable bool c_out_" << var->name_ << "_;\n";
  }

  file << "  long entry_;\n";

  file << "};\n\n";

  file << "#endif" << endl;

  file.close();
}

void WriteSource(const vector<Variable> &vars, const string name){

  ofstream file("src/"+ name +"_tree.cpp");
  file << "//File generated with generate_tree_classes.exe\n\n";

  file << "#include \""<< name <<"_tree.hpp\"\n\n";

  file << "#include <stdexcept>\n";
  file << "#include <string>\n";
  file << "#include <iostream>\n";
  file << "#include <vector>\n\n";

  file << "#include \"TROOT.h\"\n";
  file << "#include \"TTree.h\"\n";
  file << "#include \"TBranch.h\"\n";
  file << "#include \"TChain.h\"\n";
  file << "#include \"TString.h\"\n";
  file << "#include \"TObject.h\"\n";

  file << "using namespace std;\n\n";

  file << "bool "<< name <<"_tree::VectorLoader::loaded_ = false;\n\n";

  file << name <<"_tree::VectorLoader "<< name <<"_tree::vl_ = "<< name <<"_tree::VectorLoader();\n\n";

  file << name <<"_tree::VectorLoader::VectorLoader(){\n";
  file << "  if(!loaded_){\n";
  file << "    gROOT->ProcessLine(\"#include <vector>\");\n";
  file << "    loaded_ = true;\n";
  file << "  }\n";
  file << "}\n\n";

  file << name <<"_tree::"<< name <<"_tree(TString infile, TString outfile):\n";
  file << "  readOnly_(infile!=\"\" && outfile==\"\"),\n";
  file << "  writeOnly_(infile==\"\" && outfile!=\"\"),\n";
  file << "  copyUntouched_(infile==\"\" && outfile==\"\"),\n";
  file << "  bad_val_(-999.),\n";

  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "  " << var->name_ << "_(0),\n";
    }else if(Contains(var->type_, "tring")){
      file << "  " << var->name_ << "_(\"\"),\n";
    }else if(Contains(var->type_, "bool")){
      file << "  " << var->name_ << "_(false),\n";
    }else{
      file << "  " << var->name_ << "_(static_cast<" << var->type_ << ">(bad_val_)),\n";
    }
    if(Contains(var->type_, "vector")){
      file << "  p_" << var->name_ << "_(&" << var->name_ << "_),\n";
    }
    file << "  b_" << var->name_ << "_(NULL),\n";
    file << "  c_" << var->name_ << "_(false),\n";
  }

  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "  out_" << var->name_ << "_(0),\n";
    }else if(Contains(var->type_, "tring")){
      file << "  out_" << var->name_ << "_(\"\"),\n";
    }else if(Contains(var->type_, "bool")){
      file << "  out_" << var->name_ << "_(false),\n";
    }else{
      file << "  out_" << var->name_ << "_(static_cast<" << var->type_ << ">(bad_val_)),\n";
    }
    if(Contains(var->type_, "vector")){
      file << "  p_out_" << var->name_ << "_(&out_" << var->name_ << "_),\n";
    }
    file << "  c_out_" << var->name_ << "_(false),\n";
  }
  file << "  entry_(0){\n";

  file << "  if (!writeOnly_) {\n";
  file << "    intree_ = new TChain(\"tree\");\n";
  file << "    intree_->Add(infile);\n";
  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "    intree_->SetBranchAddress(\"" << var->name_ << "\", &p_" << var->name_ << "_, &b_" << var->name_ << "_);\n";
    }else{
      file << "    intree_->SetBranchAddress(\"" << var->name_ << "\", &" << var->name_ << "_, &b_" << var->name_ << "_);\n";
    }
  }
  file << "  }\n\n";
  
  file << "  if (!readOnly_) {\n";
  file << "    outfile_ = new TFile(outfile, \"recreate\", outfile, 209); // compression LZMA level 9 \n";
  file << "    if(!outfile_->IsOpen()) cout << \"Could not open output file \"<<outfile.Data()<<endl;\n";
  file << "    outfile_->cd();\n";
  file << "    if (!writeOnly_) {\n";
  file << "      outtree_ = intree_->CloneTree(0);\n";
  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "      outtree_->SetBranchAddress(\"" << var->name_ << "\", &p_out_" << var->name_ << "_);\n";
    }else{
      file << "      outtree_->SetBranchAddress(\"" << var->name_ << "\", &out_" << var->name_ << "_);\n";
    }
  }
  file << "    } else {\n";
  file << "      outtree_ = new TTree(\"tree\",\"tree\");\n";
  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "      outtree_->Branch(\"" << var->name_ << "\", &p_out_" << var->name_ << "_);\n";
    }else{
      file << "      outtree_->Branch(\"" << var->name_ << "\", &out_" << var->name_ << "_);\n";
    }
  }
  file << "    }\n"; //else
  file << "  }\n\n"; 

  file << "}\n\n"; // closing constructor bracket

  file << "void "<< name <<"_tree::Fill(){\n";
  file << "  //Loading unfilled branches so their values are copied to the output tree\n";
  file << "  if (!readOnly_ && !writeOnly_) {\n";
  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    file << "    if (!c_out_"+var->name_+"_) out_" << var->name_ << "_ = " << var->name_ << "()" <<";\n";
  }
  file << "  }\n\n";
  file << "  outtree_->Fill();\n";

  file << "  //Resetting input tree variables\n";
  file << "  if (!writeOnly_) {\n";
  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "    " << var->name_ << "_.clear();\n";
    }else if(Contains(var->type_, "tring")){
      file << "    " << var->name_ << "_ = \"\";\n";
    }else if(Contains(var->type_, "bool")){
      file << "    " << var->name_ << "_ = false;\n";
    }else{
      file << "    " << var->name_ << "_ = static_cast<" << var->type_ << ">(bad_val_);\n";
    }
  }
  file << "  }\n";

  file << "  //Resetting output tree variables\n";
  file << "  if (!readOnly_) {\n";
  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "    out_" << var->name_ << "_.clear();\n";
    }else if(Contains(var->type_, "tring")){
      file << "    out_" << var->name_ << "_ = \"\";\n";
    }else if(Contains(var->type_, "bool")){
      file << "    out_" << var->name_ << "_ = false;\n";
    }else{ 
      file << "    out_" << var->name_ << "_ = static_cast<" << var->type_ << ">(bad_val_);\n";
    }
  }

  file << "  // Reset modification trackers for output branches \n";
  for(vector<Variable>::const_iterator var = vars.begin(); var!= vars.end(); ++var){
    file << "    c_out_" << var->name_ << "_ = false;\n";
  }
  file << "  }\n";
  file << "}\n\n"; // close Fill method bracket

  file << "void "<< name <<"_tree::Write(){\n";
  file << "  outfile_->cd();\n";
  file << "  outtree_->Write(\"\",TObject::kWriteDelete);\n"; //kWriteDelete to avoid writing multiple trees
  file << "}\n\n";

  file << name <<"_tree::~"<< name <<"_tree(){\n";
  file << "  if (!readOnly_) outfile_->Close();\n";
  file << "}\n\n";
  
  file << "long "<< name <<"_tree::GetEntries() const{\n";
  file << "  return intree_->GetEntries();\n";
  file << "}\n\n";

  file << "void "<< name <<"_tree::GetEntry(const long entry){\n";
  file << "  // Reset read-trackers for input branches \n";
  for(vector<Variable>::const_iterator var = vars.begin(); var!= vars.end(); ++var){
    file << "  c_" << var->name_ << "_ = false;\n";
  }
  file << "  entry_ = intree_->LoadTree(entry);\n";
  file << "}\n\n";

  // Writing methods called when accessing input tree variable
  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    file << var->type_ << " "<< name <<"_tree::" << var->name_ << "(){\n";
    file << "  if(!c_" << var->name_ << "_ && b_" << var->name_ <<"_){\n";
    file << "    b_" << var->name_ << "_->GetEntry(entry_);\n";
    file << "    c_" << var->name_ << "_ = true;\n";
    file << "  }\n";
    file << "  return " << var->name_ << "_;\n";
    file << "}\n\n";
  }

  // Writing methods called when accessing output tree variable
  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    file << var->type_ << " & "<< name <<"_tree::out_" << var->name_ << "(){\n";
    file << "  c_out_" << var->name_ << "_ = true;\n";
    file << "  return out_" << var->name_ << "_;\n";
    file << "}\n\n";
  }

  file.close();
}



