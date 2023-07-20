#!/usr/bin/env python
import os
import ROOT
import argparse
import datetime
import re
import sys

def get_dataset_name(path):
  input_filename = os.path.basename(path)
  dataset_name = re.sub('__[0-9]+?__.*.root', '', input_filename)
  dataset_name = dataset_name.replace('*', '')
  dataset_name = dataset_name.replace('.root', '')
  return dataset_name

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-i','--input_file_glob', required=True, default = "")
  parser.add_argument('-fp', '--input_file_re_pattern', required=True, default= "")
  parser.add_argument('-fr', '--input_file_re_replace', required=True, default= "")
  parser.add_argument('-e', '--entrylist_name', required=True, default="")
  parser.add_argument('-o','--output_directory', required=True, default = "")
  parser.add_argument('-n','--nlsp_mass', required=True, default = "")
  parser.add_argument('-l','--lsp_mass', required=True, default = "")
  parser.add_argument('-ed','--entrylist_directory', default=None)
  args = parser.parse_args()

  startTime = datetime.datetime.now()

  # Find dataset name
  dataset_name = get_dataset_name(args.input_file_glob)
  output_filename = dataset_name+'.root'

  #if model == 'SMS-T5qqqqZH_HToBB-mGluino':
  #  output_filename = re.sub('mGluino.*mLSP[0-9]+to[0-9]+', 'mGluino-'+args.nlsp_mass+'_mChi-'+str(int(args.nlsp_mass)-50)+'_mLSP-'+args.lsp_mass, output_filename)
  #elif model == 'SMS-T5qqqqZH-mGluino':
  #  output_filename = re.sub('mGluino-[0-9]+to[0-9]+', 'mGluino-'+args.nlsp_mass+'_mChi-'+str(int(args.nlsp_mass)-50)+'_mLSP-'+args.lsp_mass, output_filename)
  #elif model == 'SMS-T5qqqqZH_HToBB-mN2':
  #  output_filename = re.sub('mN2.*[0-9]+to[0-9]+', 'mGluino-'+args.nlsp_mass+'_mChi-'+args.lsp_mass+'_mLSP-1', output_filename)

  output_filename = re.sub(args.input_file_re_pattern, args.input_file_re_replace.replace("-NLSP","-"+args.nlsp_mass).replace("-LSP","-"+args.lsp_mass), output_filename)
  output_filepath = args.output_directory+'/'+output_filename

  treename = 'Events'
  chain = ROOT.TChain(treename)
  chain.Add(args.input_file_glob)

  
  #entrylist_GenModel_T5qqqqZH_NLSP_LSP
  #entrylistName = "entrylist_GenModel_T5qqqqZH_"+args.nlsp_mass+"_"+args.lsp_mass
  #entrylistName = "entrylist_GenModel_T5qqqqZH_"+args.nlsp_mass+"_"+args.lsp_mass+"_1"
  entrylistName = args.entrylist_name.replace("NLSP", args.nlsp_mass).replace("LSP", args.lsp_mass)

  entrylistFile = ROOT.TFile(args.entrylist_directory+"/split_"+entrylistName+".root");
  elist = ROOT.gDirectory.Get(entrylistName)
  print('Found {} events satisfying the skim requirements.'.format(elist.GetN()))

  out_file = ROOT.TFile(output_filepath, "recreate") # in case we want compression 
  out_file.SetCompressionSettings(209); # kLZMA algo with level 9 compression
  if args.entrylist_directory:
    chain.SetEntryList(elist)
  else:
    chain.SetEventList(elist)
  print('Time to find events: '+str(datetime.datetime.now() - startTime))
  chain.CopyTree("");
  out_file.Write()
  print('Wrote file '+output_filepath)

  if args.entrylist_directory:
    entrylistFile.Close()
  chain.Reset() # Might prevent crash when closing
  out_file.Close()

  print('Total time: '+str(datetime.datetime.now() - startTime))
