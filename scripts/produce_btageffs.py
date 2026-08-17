#!/usr/bin/env python3
"""@package docstring
Generates MC photon efficiencies
"""

from argparse import ArgumentParser
from array import array
from correctionlib import schemav2 
from math import sqrt, hypot
import json
import ROOT
from glob import glob
import os
#constants
ETA_BINS = [-2.4,-1.92,-1.44,-0.96,-0.48,0.0,0.48,0.96,1.44,1.92,2.4]
#PT_BINS = [30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,500.0,600.0,800.0,1000.0]
PT_BINS = [20.0,30.0,50.0,70.0,100.0,140.0,200.0,300.0,600.0,1000.0] # changed to match BTV recommendations
DF_WPS = {'2016APV' : [0.0508, 0.2598, 0.6502],
          '2016' : [0.0480, 0.2489, 0.6377],
          '2017' : [0.0532, 0.3040, 0.7476],
          '2018' : [0.0490, 0.2783, 0.7100],
          '2022' : [0.0583, 0.3086, 0.7183],
          '2022EE' : [0.0614, 0.3196, 0.73],
          '2023' : [0.0479, 0.2431, 0.6553],
          '2023BPix' : [0.048, 0.2435, 0.6563]}
PNET_WPS = {'2022' : [0.047, 0.245, 0.6734],
            '2022EE' : [0.0499, 0.2605, 0.6915],
            '2023' : [0.0358, 0.1917, 0.6172],
            '2023BPix' : [0.0359, 0.1919, 0.6133]}
UPART_WPS = {'2024' : [0.0246, 0.1272, 0.4648],
             '2025' : [0.0246, 0.1272, 0.4648]}
WPS = {'deepflav' : DF_WPS,
       'btagpnetb' : PNET_WPS,
       'btaguptb' : UPART_WPS}

#ROOT JIT'ed C++ definitions
ROOT.gInterpreter.Declare('''
template <class C>
using RVec = ROOT::VecOps::RVec<C>;
using std::vector;

float dr(float eta1, float phi1, float eta2, float phi2) {
  float dphi = fmod(fabs(phi2-phi1), 2.0*M_PI);
  if (dphi > M_PI) {
    dphi = 2*M_PI-dphi;
  }
  return sqrt((eta2-eta1)*(eta2-eta1)+dphi*dphi);
}

RVec<bool> found_ph(RVec<bool> ph_sig, RVec<float> ph_eta, 
                    RVec<float> ph_phi, RVec<float> mc_eta, 
                    RVec<float> mc_phi) {
  vector<bool> found_photon;
  for (unsigned imc = 0; imc < mc_eta.size(); imc++) {
    bool match = false;
    for (unsigned iph = 0; iph < ph_sig.size(); iph++) {
      if (ph_sig[iph]) {
        if (dr(ph_eta[iph], ph_phi[iph], mc_eta[imc], mc_phi[imc])<0.2) {
          match = true;
        }
      }
    }
    found_photon.push_back(match);
  }
  return found_photon;
}

//checks there is a separation of at least 0.3 from the nearest (prompt) lepton
RVec<bool> mc_leptonsep(RVec<float> mc_eta, RVec<float> mc_phi, RVec<int> mc_id,
                        RVec<int> mc_statusflag) {
  vector<bool> no_nearby_lepton;
  for (unsigned imc = 0; imc < mc_eta.size(); imc++) {
    bool nolepton = true;
    for (unsigned imc2 = 0; imc2 < mc_eta.size(); imc2++) {
      if ((mc_statusflag[imc2]&0x1)!=0&&(mc_statusflag[imc2]&0x1000)!=0) {
        if (abs(mc_id[imc2])==11 || abs(mc_id[imc2])==13) {
          if (dr(mc_eta[imc], mc_phi[imc], mc_eta[imc2], mc_phi[imc2]) < 0.3) {
            nolepton = false;
          }
        }
      }
    }
    no_nearby_lepton.push_back(nolepton);
  }
  return no_nearby_lepton;
}

''')

def fix_correctionlib_json(json_texts):
  '''Fixes the format of correctionlib json created using corr.json, since 
  it is not properly formatted by default
  '''
  corr = []
  for json_text in json_texts:
    corr.append(json.loads(json_text))
  json_dict = {
    'schema_version' : 2,
    'description' : '',
    'corrections' : corr
  }
  return json.dumps(json_dict,indent=2)

def book_btagwp_hists(df, prefix, flavor_selector, algo, year):
  '''Books histograms for given flavor type

  prefix           string used to prefix this flavor
  flavor_selector  string describing how to select flavor
  algo             algorithm name in pico (deepcsv, deepflav, btaguptb)
  year             string year for b-tag WPs
  '''
  jettype = 'tru{}jet'.format(prefix)
  df = df.Define(jettype+'_eta',
                 'jet_eta[{}]'.format(flavor_selector))
  df = df.Define(jettype+'_pt',
                 'jet_pt[{}]'.format(flavor_selector))
  df = df.Define(jettype+'_{}'.format(algo),
                 'jet_{}[{}]'.format(algo,flavor_selector))
  df = df.Define(jettype+'loose_eta',
                 '{0}_eta[{0}_{2}>{1}]'.format(jettype,
                                               WPS[algo][year][0],algo))
  df = df.Define(jettype+'loose_pt',
                 '{0}_pt[{0}_{2}>{1}]'.format(jettype,
                                              WPS[algo][year][0],algo))
  df = df.Define(jettype+'med_eta',
                 '{0}_eta[{0}_{2}>{1}]'.format(jettype,
                                               WPS[algo][year][1],algo))
  df = df.Define(jettype+'med_pt',
                 '{0}_pt[{0}_{2}>{1}]'.format(jettype,
                                              WPS[algo][year][1],algo))
  df = df.Define(jettype+'tight_eta',
                 '{0}_eta[{0}_{2}>{1}]'.format(jettype,
                                               WPS[algo][year][2],algo))
  df = df.Define(jettype+'tight_pt',
                 '{0}_pt[{0}_{2}>{1}]'.format(jettype,
                                              WPS[algo][year][2],algo))
  den_hist_ptr = df.Histo2D((prefix+'den_hist',prefix+'den_hist',
                             len(ETA_BINS)-1,
                             array('d',ETA_BINS),len(PT_BINS)-1,
                             array('d',PT_BINS)),jettype+'_eta',
                             jettype+'_pt','w_lumi')
  loo_hist_ptr = df.Histo2D((prefix+'loo_hist',prefix+'loo_hist',
                             len(ETA_BINS)-1,
                             array('d',ETA_BINS),len(PT_BINS)-1,
                             array('d',PT_BINS)),jettype+'loose_eta',
                             jettype+'loose_pt','w_lumi')
  med_hist_ptr = df.Histo2D((prefix+'med_hist',prefix+'med_hist',
                             len(ETA_BINS)-1,
                             array('d',ETA_BINS),len(PT_BINS)-1,
                             array('d',PT_BINS)),jettype+'med_eta',
                             jettype+'med_pt','w_lumi')
  tig_hist_ptr = df.Histo2D((prefix+'tig_hist',prefix+'tig_hist',
                             len(ETA_BINS)-1,
                             array('d',ETA_BINS),len(PT_BINS)-1,
                             array('d',PT_BINS)),jettype+'tight_eta',
                             jettype+'tight_pt','w_lumi')
  return den_hist_ptr, loo_hist_ptr, med_hist_ptr, tig_hist_ptr

def generate_mchists(pico_filename, algo, year):
  '''Generate MC efficiency ratio plot

  pico_filename     filename of ZG pico(s) for appropriate era
  algo              btag algorithm
  year              year
  '''
  #get MC efficiencies from picos
  df = ROOT.RDataFrame('tree',pico_filename)
  df = df.Define('jet_absflavor','ROOT::VecOps::abs(jet_hflavor)')

  bden_histptr, bloo_histptr, bmed_histptr, btig_histptr = book_btagwp_hists(
          df, 'b', 'jet_absflavor==5', algo, year)
  cden_histptr, cloo_histptr, cmed_histptr, ctig_histptr = book_btagwp_hists(
          df, 'c', 'jet_absflavor==4', algo, year)
  lden_histptr, lloo_histptr, lmed_histptr, ltig_histptr = book_btagwp_hists(
          df, 'l', 'jet_absflavor<=3', algo, year)
  print('Processing n-tuples, this may take a while.')

  bden_hist = bden_histptr.GetValue()
  cden_hist = cden_histptr.GetValue()
  lden_hist = lden_histptr.GetValue()
  bloo_hist = bloo_histptr.GetValue()
  cloo_hist = cloo_histptr.GetValue()
  lloo_hist = lloo_histptr.GetValue()
  bmed_hist = bmed_histptr.GetValue()
  cmed_hist = cmed_histptr.GetValue()
  lmed_hist = lmed_histptr.GetValue()
  btig_hist = btig_histptr.GetValue()
  ctig_hist = ctig_histptr.GetValue()
  ltig_hist = ltig_histptr.GetValue()
  bloo_hist.Divide(bden_hist)
  bmed_hist.Divide(bden_hist)
  btig_hist.Divide(bden_hist)
  cloo_hist.Divide(cden_hist)
  cmed_hist.Divide(cden_hist)
  ctig_hist.Divide(cden_hist)
  lloo_hist.Divide(lden_hist)
  lmed_hist.Divide(lden_hist)
  ltig_hist.Divide(lden_hist)
  return { 'b' : [bloo_hist, bmed_hist, btig_hist],
           'c' : [cloo_hist, cmed_hist, ctig_hist],
           'l' : [lloo_hist, lmed_hist, ltig_hist] }

def make_correction(name, eff_values, unc_values):
  '''Generates correctionlib corrections with given name and values
  '''
  effs = schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['eta','pt'],
          edges=[ETA_BINS,PT_BINS],
          content=eff_values,
          flow='clamp',
          )
  uncs = schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['eta','pt'],
          edges=[ETA_BINS,PT_BINS],
          content=unc_values,
          flow='clamp',
          )
  return schemav2.Correction(
      name=name,
      version=1,
      inputs=[schemav2.Variable(name='ValType', type='string', description='effmc/systmc'),
              schemav2.Variable(name='eta', type='real', description='jet eta'),
              schemav2.Variable(name='pt', type='real', description='jet pt')],
      output=schemav2.Variable(name='sf', type='real', description=name),
      data=schemav2.Category(
          nodetype='category',
          input='ValType',
          content=[
              schemav2.CategoryItem(
                  key='effmc',
                  value=effs,
              ),
              schemav2.CategoryItem(
                  key='systmc',
                  value=uncs,
              ),
          ],
      ))

def generate_correctionlib(pico_filename, out_filename, algo, year):
  '''Generate correctionlib file with btag efficiencies

  pico_filename  filename of Z->mumu pico(s) for appropriate era
  out_filename   output filename
  '''

  #get inputs
  mc_hists = generate_mchists(pico_filename, algo, year)

  name_conversion = {'b' : 'b',
                     'c' : 'c',
                     'l' : 'uds'}

  #calculate effs and generate output
  eff_l = {'b':[],'c':[],'l':[]}
  unc_l = {'b':[],'c':[],'l':[]}
  eff_m = {'b':[],'c':[],'l':[]}
  unc_m = {'b':[],'c':[],'l':[]}
  eff_t = {'b':[],'c':[],'l':[]}
  unc_t = {'b':[],'c':[],'l':[]}
  json_l = {'b':[],'c':[],'l':[]}
  json_m = {'b':[],'c':[],'l':[]}
  json_t = {'b':[],'c':[],'l':[]}
  for jettype in ('b','c','l'):
    for ieta in range(len(ETA_BINS)-1):
      for ipt in range(len(PT_BINS)-1):
        eff_l[jettype].append(mc_hists[jettype][0].GetBinContent(ieta+1,ipt+1))
        unc_l[jettype].append(mc_hists[jettype][0].GetBinError(ieta+1,ipt+1))
        eff_m[jettype].append(mc_hists[jettype][1].GetBinContent(ieta+1,ipt+1))
        unc_m[jettype].append(mc_hists[jettype][1].GetBinError(ieta+1,ipt+1))
        eff_t[jettype].append(mc_hists[jettype][2].GetBinContent(ieta+1,ipt+1))
        unc_t[jettype].append(mc_hists[jettype][2].GetBinError(ieta+1,ipt+1))
    json_l[jettype] = make_correction('Btag_{}_WPloose_MCeff'.format(
        name_conversion[jettype]),eff_l[jettype],unc_l[jettype])
    json_m[jettype] = make_correction('Btag_{}_WPmedium_MCeff'.format(
        name_conversion[jettype]),eff_m[jettype],unc_m[jettype])
    json_t[jettype] = make_correction('Btag_{}_WPtight_MCeff'.format(
        name_conversion[jettype]),eff_t[jettype],unc_t[jettype])

  with open(out_filename,'w') as out_file:
    out_file.write(fix_correctionlib_json([
        json_l['b'].json(exclude_unset=True),
        json_l['c'].json(exclude_unset=True),
        json_l['l'].json(exclude_unset=True),
        json_m['b'].json(exclude_unset=True),
        json_m['c'].json(exclude_unset=True),
        json_m['l'].json(exclude_unset=True),
        json_t['b'].json(exclude_unset=True),
        json_t['c'].json(exclude_unset=True),
        json_t['l'].json(exclude_unset=True)]))

if __name__ == '__main__':
  parser = ArgumentParser(prog='produce_btageffs',
                          description='Generates btag efficiency files')
  parser.add_argument('-i','--input_filename')
  parser.add_argument('-y','--year')
  parser.add_argument('-a','--algorithm', choices=['deepflav','btagpnetb', 'btaguptb'])
  parser.add_argument('-o','--output_filename')
  args  = parser.parse_args()

  ROOT.EnableImplicitMT()
  
  in_file_paths = glob(os.path.join(args.input_filename+"/*.root"))
  print(in_file_paths)
  generate_correctionlib(in_file_paths, args.output_filename, 
                         args.algorithm, args.year) 
#  generate_correctionlib(args.input_filename, args.output_filename, 
# 