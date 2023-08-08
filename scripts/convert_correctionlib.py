#!/usr/bin/env python
"""Python script to generate correctionlib style corrections for trigger 
scale factors from the various formats output by the Egamma and Muon POG
official tools. 

Additional functionality to create pretty histograms may eventually also be 
added
-MO
"""

import argparse
import json
import math
import ROOT

def make_multibinning_node(abseta_bins, pt_bins, content):
  """
  Creates a JSON multibinning node

  abseta_bins   list of floats  bin boundaries for abseta
  pt_bins       list of floats  bin boundaries for pt
  content       list of floats  values for each pt-abseta bin
  """
  return {
    'nodetype' : 'multibinning',
    'inputs' : ['pt','abseta'],
    'edges' : [pt_bins,abseta_bins],
    'content' : content,
    'flow' : 'clamp',
  }

def make_data_category(category_name, values, content):
  """
  Creates a JSON category node

  category_name string         input variable name
  values        list           possible values of variable
  content       list of nodes  nodes associated to each value
  """
  content_list = []
  for i in range(len(values)):
    content_list.append({
      'key' : values[i],
      'value' : content[i]
    })
  return {
    'nodetype' : 'category',
    'input' : category_name,
    'content' : content_list
  }

#constants
EGAMMA_DESCRIPTION = 'These are electron trigger efficiencies derived using the official egamma POG tools and the nano2pico create_corrections.py script.'
MUON_DESCRIPTION = 'These are muon trigger efficiencies derived using the official muon POG tools and the nano2pico create_corrections.py script.'

if __name__ == '__main__':
  #set up argparse
  parser = argparse.ArgumentParser(
            prog='CorrectionLibGenerator',
            description='Generates correctionlib corrections')
  parser.add_argument('-i','--inputfile',help='input filename')
  parser.add_argument('-o','--outputfile',default='corrections',
            help='output filename (no extension)')
  parser.add_argument('-n','--histname',default='',
            help='histogram name in root file')
  parser.add_argument('-f','--filetype',help='type of input file')
  args = parser.parse_args()
  
  #check that input file type is supported
  supported_inputs = ['egamma_text','muon_root']
  if not args.filetype in supported_inputs:
    raise NotImplementedError('Unsupported filetype')
  
  #handle egammaEffi files
  if args.filetype == 'egamma_text':
    #read input file
    eg_txt_file = open(args.inputfile,'r')
    eg_txt = eg_txt_file.read().split('\n')
    eg_txt_file.close()
    if (eg_txt[0] != '### var1 :  el_sc_abseta ' or
              eg_txt[1] != '### var2 :  el_pt '):
      raise NotImplementedError('Currently only |eta|, pt egammaEffi binning supported')
    #parse text file
    txt_line = 2
    pt_bins = []
    abseta_bins = []
    data_effs = []
    data_uncs = []
    simu_effs = []
    simu_uncs = []
    while len(eg_txt[txt_line])>0:
      eg_txt_line = eg_txt[txt_line].split('\t')
      abseta_lowedge = float(eg_txt_line[0])
      abseta_highedge = float(eg_txt_line[1])
      pt_lowedge = float(eg_txt_line[2])
      pt_highedge = float(eg_txt_line[3])
      data_eff = float(eg_txt_line[4])
      data_statunc = float(eg_txt_line[5])
      simu_eff = float(eg_txt_line[6])
      simu_statunc = float(eg_txt_line[7])
      data_systuncsig = abs(float(eg_txt_line[8])-data_eff)
      data_systuncbkg = abs(float(eg_txt_line[9])-data_eff)
      if len(pt_bins)==0:
        #first entry
        pt_bins.append(pt_lowedge)
        pt_bins.append(pt_highedge)
        abseta_bins.append(abseta_lowedge)
        abseta_bins.append(abseta_highedge)
      else:
        if (pt_highedge > pt_bins[-1]):
          pt_bins.append(pt_highedge)
        if (abseta_highedge > abseta_bins[-1]):
          abseta_bins.append(abseta_highedge)
      data_effs.append(data_eff)
      simu_effs.append(simu_eff)
      data_uncs.append(
                math.sqrt(data_statunc**2+data_systuncsig**2+data_systuncbkg**2))
      simu_uncs.append(simu_statunc)
      txt_line += 1
    #remove directory from name in file
    clean_name = args.outputfile
    if '/' in clean_name:
      clean_name = clean_name[len(clean_name)-(clean_name[::-1]).find('/'):]
    #generate output
    json_dict = {
      'schema_version' : 2,
      'description' : EGAMMA_DESCRIPTION,
      'corrections' : [
        {
          'name' : clean_name,
          'description' : EGAMMA_DESCRIPTION,
          'version' : 1,
          'inputs' : [
            {
              'name' : 'ValType',
              'type' : 'string',
              'description' : 'effdata/systdata/effmc/systmc'
            },
            {
              'name' : 'abseta',
              'type' : 'real',
              'description' : 'absolute value of electron pseudorapidity'
            },
            {
              'name' : 'pt',
              'type' : 'real',
              'description' : 'electron transverse momentum'
            }
          ],
          'output' : {
            'name' : 'efficiency',
            'type' : 'real',
            'description' : 'trigger efficiency (uncertainty)'
          },
          'data' : make_data_category('ValType',['effdata','systdata','effmc','systmc'],
                    [make_multibinning_node(abseta_bins,pt_bins,data_effs),
                    make_multibinning_node(abseta_bins,pt_bins,data_uncs),
                    make_multibinning_node(abseta_bins,pt_bins,simu_effs),
                    make_multibinning_node(abseta_bins,pt_bins,simu_uncs)])
        }
      ]
      }
    corr_file = open(args.outputfile+'.json','w')
    corr_file.write(json.dumps(json_dict,indent=2))
    corr_file.close()
    print('Successfully wrote corrections to '+args.outputfile+'.json')
  #handle spark tnp root files
  elif args.filetype == 'muon_root':
    muon_root_file = ROOT.TFile(args.inputfile,'READ')
    #extract binning and efficiencies from ROOT file
    data_hist_name = args.inputfile[:-5]+'_efficiencyData'
    simu_hist_name = args.inputfile[:-5]+'_efficiencyMC'
    if args.histname != '':
      data_hist_name = args.histname+'_efficiencyData'
      simu_hist_name = args.histname+'_efficiencyMC'
    data_hist = muon_root_file.Get(data_hist_name)
    simu_hist = muon_root_file.Get(simu_hist_name)
    abseta_nbins = data_hist.GetXaxis().GetNbins()
    pt_nbins = data_hist.GetYaxis().GetNbins()
    abseta_bins = []
    pt_bins = []
    for i in range(abseta_nbins+1):
      abseta_bins.append(data_hist.GetXaxis().GetBinUpEdge(i))
    for i in range(pt_nbins+1):
      pt_bins.append(data_hist.GetYaxis().GetBinUpEdge(i))
    #By convention, the uncertainties already include systematic and 
    #statistical effects and furthermore are identical up and down
    data_effs = []
    data_uncs = []
    simu_effs = []
    simu_uncs = []
    for pt_bin in range(1,pt_nbins+1):
      for abseta_bin in range(1,abseta_nbins+1):
        data_effs.append(data_hist.GetBinContent(abseta_bin,pt_bin))
        data_uncs.append(data_hist.GetBinErrorUp(abseta_bin,pt_bin))
        simu_effs.append(simu_hist.GetBinContent(abseta_bin,pt_bin))
        simu_uncs.append(simu_hist.GetBinErrorUp(abseta_bin,pt_bin))
    #remove directory from name in file
    clean_name = args.outputfile
    if '/' in clean_name:
      clean_name = clean_name[len(clean_name)-(clean_name[::-1]).find('/'):]
    #generate output
    json_dict = {
      'schema_version' : 2,
      'description' : MUON_DESCRIPTION,
      'corrections' : [
        {
          'name' : clean_name,
          'description' : MUON_DESCRIPTION,
          'version' : 1,
          'inputs' : [
            {
              'name' : 'ValType',
              'type' : 'string',
              'description' : 'effdata/systdata/effmc/systmc'
            },
            {
              'name' : 'abseta',
              'type' : 'real',
              'description' : 'absolute value of muon pseudorapidity'
            },
            {
              'name' : 'pt',
              'type' : 'real',
              'description' : 'muon transverse momentum'
            }
          ],
          'output' : {
            'name' : 'efficiency',
            'type' : 'real',
            'description' : 'trigger efficiency (uncertainty)'
          },
          'data' : make_data_category('ValType',['effdata','systdata','effmc','systmc'],
                    [make_multibinning_node(abseta_bins,pt_bins,data_effs),
                    make_multibinning_node(abseta_bins,pt_bins,data_uncs),
                    make_multibinning_node(abseta_bins,pt_bins,simu_effs),
                    make_multibinning_node(abseta_bins,pt_bins,simu_uncs)])
        }
      ]
      }
    corr_file = open(args.outputfile+'.json','w')
    corr_file.write(json.dumps(json_dict,indent=2))
    corr_file.close()
    print('Successfully wrote corrections to '+args.outputfile+'.json')
    muon_root_file.Close()


