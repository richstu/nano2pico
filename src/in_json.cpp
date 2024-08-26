//Check if your run is inJSON by calling bool inJSON(VRunLumi,run,lumiblock) in the event loop..
//ported from babymaker

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "utilities.hpp"
#include "in_json.hpp"

using namespace std;

std::vector< std::vector<int> > MakeVRunLumi(std::string input){
  std::ifstream orgJSON;
  std::string fullpath = "";
  if(input == "golden2016"){
    fullpath = "txt/json/golden_Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16.json";
  } else if(input == "golden2017"){
    fullpath = "txt/json/golden_Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17.json";
  } else if(input == "golden2018"){
    fullpath = "txt/json/golden_Cert_314472-325175_13TeV_PromptReco_Collisions18.json";
  } else if(input == "goldenUL2016") {
    fullpath = "txt/json/golden_Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt";
  } else if(input == "goldenUL2017") {
    fullpath = "txt/json/golden_Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt";
  } else if(input == "goldenUL2018") {
    fullpath = "txt/json/golden_Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt";
  } else if(input == "golden2022") {
    fullpath = "txt/json/Cert_Collisions2022_355100_362760_Golden.json";
  } else if(input == "golden2023") {
    fullpath = "txt/json/Cert_Collisions2023_366442_370790_Golden.json";
  } else{
    fullpath = input;
  }
  orgJSON.open(fullpath.c_str());
  std::vector<int> VRunLumi;
  if(orgJSON.is_open()){
    char inChar;
    int inInt;
    std::string str;
    while(!orgJSON.eof()){
      char next = orgJSON.peek();
      if( next == '1' || next == '2' || next == '3' ||
          next == '4' || next == '5' || next == '6' ||
          next == '7' || next == '8' || next == '9' || 
          next == '0'){     
        orgJSON >>inInt;
        VRunLumi.push_back(inInt);        
      }
      else if(next == ' '){
        getline(orgJSON,str,' ');
      }
      else{
        orgJSON>>inChar;
      }
    }
  }//check if the file opened.
  else{
    std::cout<<"Invalid JSON File:"<<fullpath<<"!\n";
  }
  orgJSON.close();
  if(VRunLumi.size() == 0){
    std::cout<<"No Lumiblock found in JSON file\n";
  }
  std::vector< std::vector<int> > VVRunLumi;
  for(unsigned int i = 0; i+2 < VRunLumi.size();){
    if(VRunLumi[i] > 130000){
      std::vector<int> RunLumi;
      RunLumi.push_back(VRunLumi[i]);
      while(VRunLumi[i+1] < 130000 && i+1 < VRunLumi.size()){
        RunLumi.push_back(VRunLumi[i+1]);
        ++i;
      }
      VVRunLumi.push_back(RunLumi);
      ++i;
    }
  }
  return VVRunLumi;
}

bool inJSON(std::vector< std::vector<int> > VVRunLumi, int Run, int LS){
  bool answer = false;
  if(Run < 120000){
    answer = true;
  }
  else{
    for(unsigned int i = 0; i < VVRunLumi.size();++i){
      if(Run == VVRunLumi[i][0]){
        for(unsigned int j = 1; j+1 < VVRunLumi[i].size();j=j+2){
          if(LS >= VVRunLumi[i][j] && LS <= VVRunLumi[i][j+1]){
            answer = true;
          }
        }
      }
    }
  }
  return answer;
}

//void CheckVRunLumi(std::vector< std::vector<int> > VVRunLumi){
//  for(unsigned int i = 0; i < VVRunLumi.size();++i){
//    std::cout<<"Run:"<<VVRunLumi[i][0]<<" LS: ";
//    for(unsigned int j = 1; j+1 < VVRunLumi[i].size();j=j+2){
//      std::cout<<VVRunLumi[i][j]<<"-"<<VVRunLumi[i][j+1]<<" ";
//    }
//    std::cout<<std::endl;
//  }
//}
//
//void CheckVRunLumi2(std::vector< std::vector<int> > VVRunLumi){
//  for(unsigned int i = 0; i < VVRunLumi.size();++i){
//    for(unsigned int j = 1; j+1 < VVRunLumi[i].size();j=j+2){
//      if(VVRunLumi[i][j] == VVRunLumi[i][j+1]){
//        std::cout<<VVRunLumi[i][0]<<" "<<VVRunLumi[i][j]<<std::endl;
//      }
//      else{
//        for(int k=VVRunLumi[i][j];k<=VVRunLumi[i][j+1];++k){
//          std::cout<<VVRunLumi[i][0]<<" "<<k<<std::endl;
//        }
//      }
//    }
//    std::cout<<std::endl;
//  }
//}
