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

def make_multibinning_node_absetapt(abseta_bins, pt_bins, content):
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

def make_multibinning_node(input_names, bin_edges, content, flow_type='clamp'):
  """
  Creates a JSON multibinning node

  input_names  list of strings         names for input variables
  bin_edges    list of list of floats  bins for input variables
  content      list of floats          values for each bin
  flowtype     string                  overflow/underflow behavior
  """
  return {
    'nodetype' : 'multibinning',
    'inputs' : input_names,
    'edges' : bin_edges,
    'content' : content,
    'flow' : flow_type,
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

def make_valtype_input(valtypes, name='ValType'):
  """
  Creates a JSON ValType input variable node

  valtypes  list of strings  possible inputs for ValType
  name      string           title of input variable
  """
  return {
    'name' : name,
    'type' : 'string',
    'description' : '/'.join(valtypes)
  }

def write_correction_dict(name, input_variables, output_variable, data_contents):
  """
  Writes dictionary that can be converted to correctionlib json with 
  write_correctionlib_json

  name             string                name of correction set
  input_variables  list of dictionaries  input variables
  output_variable  dictionary            output variable
  data             dictionary            content of json from make_data_category or make_multibinning_node
  """
  return {
    'name' : name,
    'description' : name,
    'version' : 1,
    'inputs' : input_variables,
    'output' : output_variable,
    'data' : data
  }

def write_correctionlib_json(correction_list, out_file):
  """
  Writes corretionlib format json

  correction_list  list of dictionaries  list of outputs from write_correction_dict
  out_file         string                name of output file
  """
  json_dict = {
    'schema_version' : 2,
    'description' : '',
    'corrections' : correction_list
    }
  if (out_file[-5:]!='.json'):
    out_file = out_file+'.json'
  corr_file = open(out_file,'w')
  corr_file.write(json.dumps(json_dict,indent=2))
  corr_file.close()
  print('Successfully wrote corrections to '+out_file)

#constants
EGAMMA_DESCRIPTION = 'These are electron trigger efficiencies derived using the official egamma POG tools and the nano2pico create_corrections.py script.'
ELID_DESCRIPTION = 'These are electron ID efficiencies derived using the official egamma POG tools and the nano2pico create_corrections.py script.'
MUON_DESCRIPTION = 'These are muon trigger efficiencies derived using the official muon POG tools and the nano2pico create_corrections.py script.'
ABSETA_INPUT = {'name' : 'abseta', 'type' : 'real', 'description' : 'absolute value of pseudorapidity'}
ETA_INPUT = {'name' : 'eta', 'type' : 'real', 'description' : 'pseudorapidity'}
PT_INPUT = {'name' : 'pt', 'type' : 'real', 'description' : 'transverse momentum'}
SUPPORTED_INPUTS = ['egamma_trigger','egamma_elid','muon_root','muon_pogjson','rootfile_th2']

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
  parser.add_argument('-f','--filetype',help='type of input file: '+'/'.join(SUPPORTED_INPUTS))
  args = parser.parse_args()
  
  #check that input file type is supported
  if not args.filetype in SUPPORTED_INPUTS:
    raise NotImplementedError('Unsupported filetype')
  
  #handle egammaEffi files
  if args.filetype == 'egamma_trigger' or args.filetype == 'egamma_elid':
    #default settings for trigger
    description = EGAMMA_DESCRIPTION 
    etaname = 'abseta'
    etadesc = 'absolute value of electron pseudorapidity'
    effdesc = 'trigger efficiency (uncertainty)'
    #alternate settings for ID
    if args.filetype == 'egamma_elid':
      description = ELID_DESCRIPTION
      etaname = 'eta'
      etadesc = 'electron pseudorapidity'
      effdesc = 'electron WPL ID efficiency (uncertainty)'

    #read input file
    eg_txt_file = open(args.inputfile,'r')
    eg_txt = eg_txt_file.read().split('\n')
    eg_txt_file.close()
    if (args.filetype == 'egamma_trigger' and
        (eg_txt[0] != '### var1 :  el_sc_abseta ' or
         eg_txt[1] != '### var2 :  el_pt ')):
      raise NotImplementedError('Currently only |eta|, pt egammaEffi trigger binning supported')
    if (args.filetype == 'egamma_elid' and
        (eg_txt[0] != '### var1 :  el_sc_eta ' or
         eg_txt[1] != '### var2 :  el_pt ')):
      raise NotImplementedError('Currently only eta, pt egammaEffi ID binning supported')

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
    if args.filetype == 'egamma_elid':
      clean_name = 'ElectronWPL'

    #generate general json inputs
    input_variables = [make_valtype_input(['effdata','systdata','effmc','systmc']), ABSETA_INPUT, PT_INPUT]
    bin_variable_names = ['pt','abseta']
    output_variable = {'name' : 'efficiency', 'type' : 'real', 'description' : effdesc}
    if args.filetype == 'egamma_elid':
      input_variables = [make_valtype_input(['effdata','systdata','effmc','systmc']), ETA_INPUT, PT_INPUT]
      bin_variable_names = ['pt','eta']

    #generate output and write to file
    data = make_data_category('ValType',['effdata','systdata','effmc','systmc'],
        [make_multibinning_node(bin_variable_names,[pt_bins,abseta_bins],data_effs),
        make_multibinning_node(bin_variable_names,[pt_bins,abseta_bins],data_uncs),
        make_multibinning_node(bin_variable_names,[pt_bins,abseta_bins],simu_effs),
        make_multibinning_node(bin_variable_names,[pt_bins,abseta_bins],simu_uncs)])
    correction_list = [write_correction_dict(clean_name, input_variables, output_variable, data)]
    write_correctionlib_json(correction_list, args.outputfile)

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
                    [make_multibinning_node_absetapt(abseta_bins,pt_bins,data_effs),
                    make_multibinning_node_absetapt(abseta_bins,pt_bins,data_uncs),
                    make_multibinning_node_absetapt(abseta_bins,pt_bins,simu_effs),
                    make_multibinning_node_absetapt(abseta_bins,pt_bins,simu_uncs)])
        }
      ]
      }
    corr_file = open(args.outputfile+'.json','w')
    corr_file.write(json.dumps(json_dict,indent=2))
    corr_file.close()
    print('Successfully wrote corrections to '+args.outputfile+'.json')
    muon_root_file.Close()

  #handle general root files
  elif args.filetype == 'rootfile_th2':
    #get histogram from ROOT file
    root_file = ROOT.TFile(args.inputfile,'READ')
    if args.histname == '':
      raise ValueError('Must supply a histogram name with rootfile')
    hist = root_file.Get(args.histname)

    #process histogram to extract binning and efficiencies
    x_nbins = hist.GetXaxis().GetNbins()
    y_nbins = hist.GetYaxis().GetNbins()
    x_bins = []
    y_bins = []
    for i in range(x_nbins+1):
      x_bins.append(hist.GetXaxis().GetBinUpEdge(i))
    for i in range(y_nbins+1):
      y_bins.append(hist.GetYaxis().GetBinUpEdge(i))
    effs = []
    uncs = []
    for y_bin in range(1,y_nbins+1):
      for x_bin in range(1,x_nbins+1):
        effs.append(hist.GetBinContent(x_bin,y_bin))
        uncs.append(hist.GetBinErrorUp(x_bin,y_bin))

    #generate general json inputs
    #NOTE: this is currently hard-coded for a specific histogram input, 
    #haven't got the time to make it general yet
    input_variables = [make_valtype_input(['effmc','systmc']), ETA_INPUT, PT_INPUT]
    output_variable = {'name' : 'eff', 'type' : 'real', 'description' : 'electron reco/ID/Iso efficiency/uncertainty'}
    bin_variable_names = ['pt','eta']

    #generate output
    #NOTE: see above about hardcoded output
    data = make_data_category('ValType',['effmc','systmc'],
            [make_multibinning_node(bin_variable_names,[y_bins,x_bins],effs),
            make_multibinning_node(bin_variable_names,[y_bins,x_bins],uncs)])
    correction_list = [write_correction_dict('Electron_WPL_MCeff', input_variables, output_variable, data)]
    write_correctionlib_json(correction_list, args.outputfile)

  #muon POG json format
  elif args.filetype == 'muon_pogjson':
    #read input file
    muon_json_file = open(args.inputfile,'r')
    muon_json = json.load(muon_json_file)
    muon_json_file.close()

    #turn muon json into structure that can be used to sort data into 
    #correctionlib format
    corrections = []
    for correction_name in muon_json:
      eta_values = []
      for eta_bin in muon_json[correction_name]['abseta_pt']:
        eta_bound_lower = float(eta_bin[eta_bin.find('[')+1:eta_bin.find(',')])
        eta_bound_upper = float(eta_bin[eta_bin.find(',')+1:eta_bin.find(']')])
        pt_values = []
        for pt_bin in muon_json[correction_name]['abseta_pt'][eta_bin]:
          pt_bound_lower = float(pt_bin[pt_bin.find('[')+1:pt_bin.find(',')])
          pt_bound_upper = float(pt_bin[pt_bin.find(',')+1:pt_bin.find(']')])
          sf = muon_json[correction_name]['abseta_pt'][eta_bin][pt_bin]['value']
          unc = muon_json[correction_name]['abseta_pt'][eta_bin][pt_bin]['error']
          pt_values.append([pt_bound_lower, pt_bound_upper, float(sf), float(unc)])
        #sort pt bins
        pt_values.sort(key = lambda x : x[0])
        eta_values.append([eta_bound_lower, eta_bound_upper, pt_values])
      #sort eta bins
      eta_values.sort(key = lambda x : x[0])
      corrections.append([correction_name,eta_values])

    #generate general json inputs
    input_variables = [make_valtype_input(['sf','unc']), ABSETA_INPUT, PT_INPUT]
    output_variable = {'name' : 'sf', 'type' : 'real', 'description' : 'muon scale factor'}
    bin_variable_names = ['pt','abseta']

    #generate output
    correction_list = []
    for correction in corrections:
      abseta_bins = [correction[1][n][0] for n in range(len(correction[1]))] + [correction[1][-1][1]]
      pt_bins = [correction[1][0][2][n][0] for n in range(len(correction[1][0][2]))] + [correction[1][0][2][-1][1]]
      sfs = []
      uncs = []
      for ipt in range(len(correction[1][0][2])):
        for ieta in range(len(correction[1])):
          sfs.append(correction[1][ieta][2][ipt][2])
          uncs.append(correction[1][ieta][2][ipt][3])
      data = make_data_category('ValType',['sf','unc'],
              [make_multibinning_node(bin_variable_names,[pt_bins,abseta_bins],sfs),
              make_multibinning_node(bin_variable_names,[pt_bins,abseta_bins],uncs)])
      correction_list.append(write_correction_dict(correction[0], input_variables, output_variable, data))
    write_correctionlib_json(correction_list, args.outputfile)
