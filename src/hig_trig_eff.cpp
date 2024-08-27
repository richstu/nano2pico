#include "hig_trig_eff.hpp"

// #include <cmath>

// #include <deque>
#include <iostream>
// #include <fstream>
// #include <string>
// #include <stdexcept>
// #include <iomanip>   // setw

// #include <libgen.h>

// #include "TCollection.h"
// #include "TFile.h"
// #include "TGraph.h"
// #include "TH1D.h"
// #include "TList.h"
// #include "TString.h"
// #include "TSystemDirectory.h"
// #include "TSystemFile.h"
// #include "TSystem.h"
// #include "TTree.h"
// #include "TChain.h"
// #include "TRegexp.h"

using namespace std;

namespace hig_trig_eff{

  float eff(pico_tree &pico){
      float errup(0), errdown(0); // Not used, but for reference
      errup+=errdown;
      float eff = 1., met = pico.out_met(), ht = pico.out_ht(); //Note that these are floats being compared to literals (default double. . .)
      if(ht>   0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.532; errup = 0.013; errdown = 0.013;}
      else if(ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff = 0.612; errup = 0.005; errdown = 0.005;}
      else if(ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.589; errup = 0.023; errdown = 0.024;}
      else if(ht> 800 && ht<=1000 && met> 150 && met<= 155) {eff = 0.515; errup = 0.042; errdown = 0.042;}
      else if(ht>1000 && ht<=9999 && met> 150 && met<= 155) {eff = 0.588; errup = 0.052; errdown = 0.054;}
      else if(ht>   0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.591; errup = 0.014; errdown = 0.014;}
      else if(ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff = 0.684; errup = 0.005; errdown = 0.005;}
      else if(ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.678; errup = 0.022; errdown = 0.022;}
      else if(ht> 800 && ht<=1000 && met> 155 && met<= 160) {eff = 0.537; errup = 0.042; errdown = 0.042;}
      else if(ht>1000 && ht<=9999 && met> 155 && met<= 160) {eff = 0.511; errup = 0.057; errdown = 0.057;}
      else if(ht>   0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.619; errup = 0.016; errdown = 0.016;}
      else if(ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff = 0.727; errup = 0.006; errdown = 0.006;}
      else if(ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.699; errup = 0.022; errdown = 0.023;}
      else if(ht> 800 && ht<=1000 && met> 160 && met<= 165) {eff = 0.690; errup = 0.039; errdown = 0.041;}
      else if(ht>1000 && ht<=9999 && met> 160 && met<= 165) {eff = 0.568; errup = 0.048; errdown = 0.049;}
      else if(ht>   0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.678; errup = 0.017; errdown = 0.018;}
      else if(ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff = 0.769; errup = 0.006; errdown = 0.006;}
      else if(ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.732; errup = 0.022; errdown = 0.023;}
      else if(ht> 800 && ht<=1000 && met> 165 && met<= 170) {eff = 0.609; errup = 0.042; errdown = 0.044;}
      else if(ht>1000 && ht<=9999 && met> 165 && met<= 170) {eff = 0.685; errup = 0.058; errdown = 0.064;}
      else if(ht>   0 && ht<= 200 && met> 170 && met<= 175) {eff = 0.670; errup = 0.019; errdown = 0.020;}
      else if(ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff = 0.811; errup = 0.005; errdown = 0.006;}
      else if(ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.779; errup = 0.022; errdown = 0.024;}
      else if(ht> 800 && ht<=1000 && met> 170 && met<= 175) {eff = 0.736; errup = 0.041; errdown = 0.045;}
      else if(ht>1000 && ht<=9999 && met> 170 && met<= 175) {eff = 0.663; errup = 0.056; errdown = 0.061;}
      else if(ht>   0 && ht<= 200 && met> 175 && met<= 180) {eff = 0.730; errup = 0.020; errdown = 0.021;}
      else if(ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff = 0.838; errup = 0.005; errdown = 0.006;}
      else if(ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.820; errup = 0.022; errdown = 0.024;}
      else if(ht> 800 && ht<=1000 && met> 175 && met<= 180) {eff = 0.819; errup = 0.037; errdown = 0.043;}
      else if(ht>1000 && ht<=9999 && met> 175 && met<= 180) {eff = 0.736; errup = 0.055; errdown = 0.062;}
      else if(ht>   0 && ht<= 200 && met> 180 && met<= 185) {eff = 0.745; errup = 0.023; errdown = 0.024;}
      else if(ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff = 0.874; errup = 0.005; errdown = 0.005;}
      else if(ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.848; errup = 0.021; errdown = 0.023;}
      else if(ht> 800 && ht<=1000 && met> 180 && met<= 185) {eff = 0.869; errup = 0.038; errdown = 0.048;}
      else if(ht>1000 && ht<=9999 && met> 180 && met<= 185) {eff = 0.759; errup = 0.048; errdown = 0.055;}
      else if(ht>   0 && ht<= 200 && met> 185 && met<= 190) {eff = 0.777; errup = 0.024; errdown = 0.026;}
      else if(ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff = 0.903; errup = 0.005; errdown = 0.005;}
      else if(ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.850; errup = 0.021; errdown = 0.023;}
      else if(ht> 800 && ht<=1000 && met> 185 && met<= 190) {eff = 0.839; errup = 0.041; errdown = 0.049;}
      else if(ht>1000 && ht<=9999 && met> 185 && met<= 190) {eff = 0.847; errup = 0.044; errdown = 0.055;}
      else if(ht>   0 && ht<= 200 && met> 190 && met<= 195) {eff = 0.792; errup = 0.026; errdown = 0.028;}
      else if(ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff = 0.907; errup = 0.005; errdown = 0.005;}
      else if(ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.884; errup = 0.020; errdown = 0.023;}
      else if(ht> 800 && ht<=1000 && met> 190 && met<= 195) {eff = 0.870; errup = 0.036; errdown = 0.045;}
      else if(ht>1000 && ht<=9999 && met> 190 && met<= 195) {eff = 0.781; errup = 0.051; errdown = 0.059;}
      else if(ht>   0 && ht<= 200 && met> 195 && met<= 200) {eff = 0.757; errup = 0.033; errdown = 0.036;}
      else if(ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff = 0.924; errup = 0.005; errdown = 0.005;}
      else if(ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.921; errup = 0.016; errdown = 0.020;}
      else if(ht> 800 && ht<=1000 && met> 195 && met<= 200) {eff = 0.936; errup = 0.027; errdown = 0.041;}
      else if(ht>1000 && ht<=9999 && met> 195 && met<= 200) {eff = 0.803; errup = 0.049; errdown = 0.059;}
      else if(ht>   0 && ht<= 200 && met> 200 && met<= 210) {eff = 0.841; errup = 0.022; errdown = 0.025;}
      else if(ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff = 0.949; errup = 0.003; errdown = 0.003;}
      else if(ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.927; errup = 0.013; errdown = 0.015;}
      else if(ht> 800 && ht<=1000 && met> 200 && met<= 210) {eff = 0.894; errup = 0.023; errdown = 0.027;}
      else if(ht>1000 && ht<=9999 && met> 200 && met<= 210) {eff = 0.839; errup = 0.036; errdown = 0.042;}
      else if(ht>   0 && ht<= 200 && met> 210 && met<= 220) {eff = 0.850; errup = 0.028; errdown = 0.032;}
      else if(ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff = 0.966; errup = 0.003; errdown = 0.003;}
      else if(ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.952; errup = 0.011; errdown = 0.013;}
      else if(ht> 800 && ht<=1000 && met> 210 && met<= 220) {eff = 0.919; errup = 0.024; errdown = 0.031;}
      else if(ht>1000 && ht<=9999 && met> 210 && met<= 220) {eff = 0.959; errup = 0.018; errdown = 0.027;}
      else if(ht>   0 && ht<= 200 && met> 220 && met<= 230) {eff = 0.896; errup = 0.029; errdown = 0.037;}
      else if(ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff = 0.973; errup = 0.003; errdown = 0.003;}
      else if(ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.979; errup = 0.008; errdown = 0.011;}
      else if(ht> 800 && ht<=1000 && met> 220 && met<= 230) {eff = 0.956; errup = 0.017; errdown = 0.025;}
      else if(ht>1000 && ht<=9999 && met> 220 && met<= 230) {eff = 0.971; errup = 0.016; errdown = 0.027;}
      else if(ht>   0 && ht<= 200 && met> 230 && met<= 240) {eff = 0.844; errup = 0.042; errdown = 0.053;}
      else if(ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff = 0.983; errup = 0.003; errdown = 0.003;}
      else if(ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.976; errup = 0.009; errdown = 0.013;}
      else if(ht> 800 && ht<=1000 && met> 230 && met<= 240) {eff = 0.983; errup = 0.011; errdown = 0.022;}
      else if(ht>1000 && ht<=9999 && met> 230 && met<= 240) {eff = 0.942; errup = 0.028; errdown = 0.043;}
      else if(ht>   0 && ht<= 200 && met> 240 && met<= 250) {eff = 0.880; errup = 0.047; errdown = 0.065;}
      else if(ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff = 0.985; errup = 0.003; errdown = 0.003;}
      else if(ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.992; errup = 0.005; errdown = 0.010;}
      else if(ht> 800 && ht<=1000 && met> 240 && met<= 250) {eff = 0.989; errup = 0.010; errdown = 0.026;}
      else if(ht>1000 && ht<=9999 && met> 240 && met<= 250) {eff = 0.931; errup = 0.030; errdown = 0.044;}
      else if(ht>   0 && ht<= 200 && met> 250 && met<= 275) {eff = 0.915; errup = 0.033; errdown = 0.047;}
      else if(ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff = 0.989; errup = 0.002; errdown = 0.002;}
      else if(ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.992; errup = 0.004; errdown = 0.006;}
      else if(ht> 800 && ht<=1000 && met> 250 && met<= 275) {eff = 0.984; errup = 0.009; errdown = 0.016;}
      else if(ht>1000 && ht<=9999 && met> 250 && met<= 275) {eff = 0.965; errup = 0.015; errdown = 0.023;}
      else if(ht>   0 && ht<= 200 && met> 275 && met<= 300) {eff = 0.862; errup = 0.065; errdown = 0.096;}
      else if(ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff = 0.992; errup = 0.002; errdown = 0.003;}
      else if(ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.989; errup = 0.005; errdown = 0.008;}
      else if(ht> 800 && ht<=1000 && met> 275 && met<= 300) {eff = 0.963; errup = 0.016; errdown = 0.024;}
      else if(ht>1000 && ht<=9999 && met> 275 && met<= 300) {eff = 0.991; errup = 0.007; errdown = 0.020;}
      else if(ht>   0 && ht<= 200 && met> 300 && met<=9999) {eff = 0.744; errup = 0.074; errdown = 0.089;}
      else if(ht> 200 && ht<= 600 && met> 300 && met<=9999) {eff = 0.994; errup = 0.001; errdown = 0.002;}
      else if(ht> 600 && ht<= 800 && met> 300 && met<=9999) {eff = 0.996; errup = 0.002; errdown = 0.003;}
      else if(ht> 800 && ht<=1000 && met> 300 && met<=9999) {eff = 1.000; errup = 0.000; errdown = 0.003;}
      else if(ht>1000 && ht<=9999 && met> 300 && met<=9999) {eff = 0.987; errup = 0.005; errdown = 0.008;}
      std::cout<<ht<<" "<<met<<" "<<eff<<std::endl;

      return eff;
  }

  float eff_unc(pico_tree &pico){
      float errup(0), errdown(0); // Not used, but for reference
      errup+=errdown;
      float uncert = 0., met = pico.out_met(), ht = pico.out_ht();
      if(ht>   0 && ht<= 200 && met> 150 && met<= 155) {uncert = 0.072; errup = 0.013; errdown = 0.013;}
      else if(ht> 200 && ht<= 600 && met> 150 && met<= 155) {uncert = 0.071; errup = 0.005; errdown = 0.005;}
      else if(ht> 600 && ht<= 800 && met> 150 && met<= 155) {uncert = 0.075; errup = 0.023; errdown = 0.024;}
      else if(ht> 800 && ht<=1000 && met> 150 && met<= 155) {uncert = 0.083; errup = 0.042; errdown = 0.042;}
      else if(ht>1000 && ht<=9999 && met> 150 && met<= 155) {uncert = 0.089; errup = 0.052; errdown = 0.054;}
      else if(ht>   0 && ht<= 200 && met> 155 && met<= 160) {uncert = 0.072; errup = 0.014; errdown = 0.014;}
      else if(ht> 200 && ht<= 600 && met> 155 && met<= 160) {uncert = 0.071; errup = 0.005; errdown = 0.005;}
      else if(ht> 600 && ht<= 800 && met> 155 && met<= 160) {uncert = 0.074; errup = 0.022; errdown = 0.022;}
      else if(ht> 800 && ht<=1000 && met> 155 && met<= 160) {uncert = 0.083; errup = 0.042; errdown = 0.042;}
      else if(ht>1000 && ht<=9999 && met> 155 && met<= 160) {uncert = 0.091; errup = 0.057; errdown = 0.057;}
      else if(ht>   0 && ht<= 200 && met> 160 && met<= 165) {uncert = 0.047; errup = 0.016; errdown = 0.016;}
      else if(ht> 200 && ht<= 600 && met> 160 && met<= 165) {uncert = 0.044; errup = 0.006; errdown = 0.006;}
      else if(ht> 600 && ht<= 800 && met> 160 && met<= 165) {uncert = 0.050; errup = 0.022; errdown = 0.023;}
      else if(ht> 800 && ht<=1000 && met> 160 && met<= 165) {uncert = 0.060; errup = 0.039; errdown = 0.041;}
      else if(ht>1000 && ht<=9999 && met> 160 && met<= 165) {uncert = 0.066; errup = 0.048; errdown = 0.049;}
      else if(ht>   0 && ht<= 200 && met> 165 && met<= 170) {uncert = 0.047; errup = 0.017; errdown = 0.018;}
      else if(ht> 200 && ht<= 600 && met> 165 && met<= 170) {uncert = 0.044; errup = 0.006; errdown = 0.006;}
      else if(ht> 600 && ht<= 800 && met> 165 && met<= 170) {uncert = 0.050; errup = 0.022; errdown = 0.023;}
      else if(ht> 800 && ht<=1000 && met> 165 && met<= 170) {uncert = 0.062; errup = 0.042; errdown = 0.044;}
      else if(ht>1000 && ht<=9999 && met> 165 && met<= 170) {uncert = 0.077; errup = 0.058; errdown = 0.064;}
      else if(ht>   0 && ht<= 200 && met> 170 && met<= 175) {uncert = 0.048; errup = 0.019; errdown = 0.020;}
      else if(ht> 200 && ht<= 600 && met> 170 && met<= 175) {uncert = 0.044; errup = 0.005; errdown = 0.006;}
      else if(ht> 600 && ht<= 800 && met> 170 && met<= 175) {uncert = 0.050; errup = 0.022; errdown = 0.024;}
      else if(ht> 800 && ht<=1000 && met> 170 && met<= 175) {uncert = 0.063; errup = 0.041; errdown = 0.045;}
      else if(ht>1000 && ht<=9999 && met> 170 && met<= 175) {uncert = 0.075; errup = 0.056; errdown = 0.061;}
      else if(ht>   0 && ht<= 200 && met> 175 && met<= 180) {uncert = 0.049; errup = 0.020; errdown = 0.021;}
      else if(ht> 200 && ht<= 600 && met> 175 && met<= 180) {uncert = 0.044; errup = 0.005; errdown = 0.006;}
      else if(ht> 600 && ht<= 800 && met> 175 && met<= 180) {uncert = 0.050; errup = 0.022; errdown = 0.024;}
      else if(ht> 800 && ht<=1000 && met> 175 && met<= 180) {uncert = 0.062; errup = 0.037; errdown = 0.043;}
      else if(ht>1000 && ht<=9999 && met> 175 && met<= 180) {uncert = 0.076; errup = 0.055; errdown = 0.062;}
      else if(ht>   0 && ht<= 200 && met> 180 && met<= 185) {uncert = 0.049; errup = 0.023; errdown = 0.024;}
      else if(ht> 200 && ht<= 600 && met> 180 && met<= 185) {uncert = 0.042; errup = 0.005; errdown = 0.005;}
      else if(ht> 600 && ht<= 800 && met> 180 && met<= 185) {uncert = 0.048; errup = 0.021; errdown = 0.023;}
      else if(ht> 800 && ht<=1000 && met> 180 && met<= 185) {uncert = 0.064; errup = 0.038; errdown = 0.048;}
      else if(ht>1000 && ht<=9999 && met> 180 && met<= 185) {uncert = 0.069; errup = 0.048; errdown = 0.055;}
      else if(ht>   0 && ht<= 200 && met> 185 && met<= 190) {uncert = 0.049; errup = 0.024; errdown = 0.026;}
      else if(ht> 200 && ht<= 600 && met> 185 && met<= 190) {uncert = 0.042; errup = 0.005; errdown = 0.005;}
      else if(ht> 600 && ht<= 800 && met> 185 && met<= 190) {uncert = 0.048; errup = 0.021; errdown = 0.023;}
      else if(ht> 800 && ht<=1000 && met> 185 && met<= 190) {uncert = 0.065; errup = 0.041; errdown = 0.049;}
      else if(ht>1000 && ht<=9999 && met> 185 && met<= 190) {uncert = 0.069; errup = 0.044; errdown = 0.055;}
      else if(ht>   0 && ht<= 200 && met> 190 && met<= 195) {uncert = 0.051; errup = 0.026; errdown = 0.028;}
      else if(ht> 200 && ht<= 600 && met> 190 && met<= 195) {uncert = 0.042; errup = 0.005; errdown = 0.005;}
      else if(ht> 600 && ht<= 800 && met> 190 && met<= 195) {uncert = 0.048; errup = 0.020; errdown = 0.023;}
      else if(ht> 800 && ht<=1000 && met> 190 && met<= 195) {uncert = 0.062; errup = 0.036; errdown = 0.045;}
      else if(ht>1000 && ht<=9999 && met> 190 && met<= 195) {uncert = 0.073; errup = 0.051; errdown = 0.059;}
      else if(ht>   0 && ht<= 200 && met> 195 && met<= 200) {uncert = 0.055; errup = 0.033; errdown = 0.036;}
      else if(ht> 200 && ht<= 600 && met> 195 && met<= 200) {uncert = 0.042; errup = 0.005; errdown = 0.005;}
      else if(ht> 600 && ht<= 800 && met> 195 && met<= 200) {uncert = 0.046; errup = 0.016; errdown = 0.020;}
      else if(ht> 800 && ht<=1000 && met> 195 && met<= 200) {uncert = 0.059; errup = 0.027; errdown = 0.041;}
      else if(ht>1000 && ht<=9999 && met> 195 && met<= 200) {uncert = 0.072; errup = 0.049; errdown = 0.059;}
      else if(ht>   0 && ht<= 200 && met> 200 && met<= 210) {uncert = 0.035; errup = 0.022; errdown = 0.025;}
      else if(ht> 200 && ht<= 600 && met> 200 && met<= 210) {uncert = 0.024; errup = 0.003; errdown = 0.003;}
      else if(ht> 600 && ht<= 800 && met> 200 && met<= 210) {uncert = 0.028; errup = 0.013; errdown = 0.015;}
      else if(ht> 800 && ht<=1000 && met> 200 && met<= 210) {uncert = 0.036; errup = 0.023; errdown = 0.027;}
      else if(ht>1000 && ht<=9999 && met> 200 && met<= 210) {uncert = 0.049; errup = 0.036; errdown = 0.042;}
      else if(ht>   0 && ht<= 200 && met> 210 && met<= 220) {uncert = 0.040; errup = 0.028; errdown = 0.032;}
      else if(ht> 200 && ht<= 600 && met> 210 && met<= 220) {uncert = 0.024; errup = 0.003; errdown = 0.003;}
      else if(ht> 600 && ht<= 800 && met> 210 && met<= 220) {uncert = 0.027; errup = 0.011; errdown = 0.013;}
      else if(ht> 800 && ht<=1000 && met> 210 && met<= 220) {uncert = 0.039; errup = 0.024; errdown = 0.031;}
      else if(ht>1000 && ht<=9999 && met> 210 && met<= 220) {uncert = 0.036; errup = 0.018; errdown = 0.027;}
      else if(ht>   0 && ht<= 200 && met> 220 && met<= 230) {uncert = 0.044; errup = 0.029; errdown = 0.037;}
      else if(ht> 200 && ht<= 600 && met> 220 && met<= 230) {uncert = 0.024; errup = 0.003; errdown = 0.003;}
      else if(ht> 600 && ht<= 800 && met> 220 && met<= 230) {uncert = 0.026; errup = 0.008; errdown = 0.011;}
      else if(ht> 800 && ht<=1000 && met> 220 && met<= 230) {uncert = 0.035; errup = 0.017; errdown = 0.025;}
      else if(ht>1000 && ht<=9999 && met> 220 && met<= 230) {uncert = 0.036; errup = 0.016; errdown = 0.027;}
      else if(ht>   0 && ht<= 200 && met> 230 && met<= 240) {uncert = 0.055; errup = 0.042; errdown = 0.053;}
      else if(ht> 200 && ht<= 600 && met> 230 && met<= 240) {uncert = 0.015; errup = 0.003; errdown = 0.003;}
      else if(ht> 600 && ht<= 800 && met> 230 && met<= 240) {uncert = 0.020; errup = 0.009; errdown = 0.013;}
      else if(ht> 800 && ht<=1000 && met> 230 && met<= 240) {uncert = 0.027; errup = 0.011; errdown = 0.022;}
      else if(ht>1000 && ht<=9999 && met> 230 && met<= 240) {uncert = 0.046; errup = 0.028; errdown = 0.043;}
      else if(ht>   0 && ht<= 200 && met> 240 && met<= 250) {uncert = 0.067; errup = 0.047; errdown = 0.065;}
      else if(ht> 200 && ht<= 600 && met> 240 && met<= 250) {uncert = 0.015; errup = 0.003; errdown = 0.003;}
      else if(ht> 600 && ht<= 800 && met> 240 && met<= 250) {uncert = 0.018; errup = 0.005; errdown = 0.010;}
      else if(ht> 800 && ht<=1000 && met> 240 && met<= 250) {uncert = 0.030; errup = 0.010; errdown = 0.026;}
      else if(ht>1000 && ht<=9999 && met> 240 && met<= 250) {uncert = 0.047; errup = 0.030; errdown = 0.044;}
      else if(ht>   0 && ht<= 200 && met> 250 && met<= 275) {uncert = 0.049; errup = 0.033; errdown = 0.047;}
      else if(ht> 200 && ht<= 600 && met> 250 && met<= 275) {uncert = 0.012; errup = 0.002; errdown = 0.002;}
      else if(ht> 600 && ht<= 800 && met> 250 && met<= 275) {uncert = 0.013; errup = 0.004; errdown = 0.006;}
      else if(ht> 800 && ht<=1000 && met> 250 && met<= 275) {uncert = 0.020; errup = 0.009; errdown = 0.016;}
      else if(ht>1000 && ht<=9999 && met> 250 && met<= 275) {uncert = 0.026; errup = 0.015; errdown = 0.023;}
      else if(ht>   0 && ht<= 200 && met> 275 && met<= 300) {uncert = 0.096; errup = 0.065; errdown = 0.096;}
      else if(ht> 200 && ht<= 600 && met> 275 && met<= 300) {uncert = 0.012; errup = 0.002; errdown = 0.003;}
      else if(ht> 600 && ht<= 800 && met> 275 && met<= 300) {uncert = 0.015; errup = 0.005; errdown = 0.008;}
      else if(ht> 800 && ht<=1000 && met> 275 && met<= 300) {uncert = 0.027; errup = 0.016; errdown = 0.024;}
      else if(ht>1000 && ht<=9999 && met> 275 && met<= 300) {uncert = 0.023; errup = 0.007; errdown = 0.020;}
      else if(ht>   0 && ht<= 200 && met> 300 && met<=9999) {uncert = 0.089; errup = 0.074; errdown = 0.089;}
      else if(ht> 200 && ht<= 600 && met> 300 && met<=9999) {uncert = 0.006; errup = 0.001; errdown = 0.002;}
      else if(ht> 600 && ht<= 800 && met> 300 && met<=9999) {uncert = 0.007; errup = 0.002; errdown = 0.003;}
      else if(ht> 800 && ht<=1000 && met> 300 && met<=9999) {uncert = 0.007; errup = 0.000; errdown = 0.003;}
      else if(ht>1000 && ht<=9999 && met> 300 && met<=9999) {uncert = 0.010; errup = 0.005; errdown = 0.008;}

      return uncert;
  }
}
