#include <iostream>
#include <getopt.h>
#include "cross_sections.hpp"
#include "utilities.hpp"

using namespace std;

int GetHiggsinoMass(const string &path);
int GetGluinoMass(const string &path);
int GetSUSYMass(const string & patth, const string & tag);
int GetHiggsinoMass(const string &path){
  string key = "_mChi-";
  // if (Contains(path, "T2tt")) key = "_mStop-"; 
  auto pos1 = path.rfind(key)+key.size();
  auto pos2 = path.find("_", pos1);
  string mass_string = path.substr(pos1, pos2-pos1);
  int unrounded_mass = stoi(mass_string);
  int rounded_mass = unrounded_mass;
  if (unrounded_mass != 127)
    rounded_mass = ((unrounded_mass+12)/25)*25;
  return rounded_mass;
}
// Tag examples: "_mChi-","-mGluino-" 
int GetSUSYMass(const string & path, const string & tag) {
  auto pos1 = path.rfind(tag)+tag.size();
  auto pos2 = path.find("_", pos1);
  string mass_string = path.substr(pos1, pos2-pos1);
  int unrounded_mass = stoi(mass_string);
  return unrounded_mass;
}
int GetGluinoMass(const string & path) {
  return GetSUSYMass(path, "-mGluino-");
}

namespace {
  std::string filename = "";
  int year = -1;
}

void GetOptions(int argc, char *argv[]);
int main(int argc, char *argv[]){
  GetOptions(argc, argv);
  if (filename=="" || year==-1) {
    std::cout<<"ERROR: filename is not given or year is not specified"<<std::endl;
    exit(1);
  }

  double xsec = -9999;
  if (Contains(filename, "SMS-TChiHH_HToAll")){
    double exsec(0.);
    int mglu = GetHiggsinoMass(filename);
    xsec::higgsinoCrossSection(mglu, xsec, exsec, year);
    xsec = xsec / .5824/.5824; // Remove H to bb branch ratio
    exsec = exsec / .5824/.5824; // Remove H to bb branch ratio
  } else if (Contains(filename, "SMS-TChi")){
    double exsec(0.);
    int mglu = GetHiggsinoMass(filename);
    xsec::higgsinoCrossSection(mglu, xsec, exsec, year);
  } else if (Contains(filename, "SMS-T5qqqqZH_HToBB")) {
    double exsec(0.);
    int mglu = GetGluinoMass(filename);
    xsec::gluinoCrossSection(mglu, xsec, exsec, year);
    xsec = xsec * .5824*.5824; // Add in H to bb branch ratio
    exsec = exsec * .5824*.5824; // Add in H to bb branch ratio
  } else if (Contains(filename, "SMS-T5qqqqZH-")) {
    double exsec(0.);
    int mglu = GetGluinoMass(filename);
    xsec::gluinoCrossSection(mglu, xsec, exsec, year);
  }else{
    xsec = xsec::crossSection(filename, year);  
  }

  std::cout<<"Cross-section for "<<filename<<" is "<<xsec<<" pb"<<std::endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"filename", required_argument, 0,'f'}, 
      {"year",  required_argument, 0,'y'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "f:y:", long_options, &option_index);
    if(opt == -1) break;

    std::string optname;
    switch(opt){
    case 'f':
      filename = optarg;
      break;
    case 'y':
      year = std::atoi(optarg);
      break;
    case 0:
      optname = long_options[option_index].name;
      printf("Bad option! Found option name %s\n", optname.c_str());
      exit(1);
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      exit(1);
      break;
    }
  }
}
