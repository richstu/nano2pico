#!/bin/env python
import os
import ROOT
import argparse

# Input:   root_files/SMS-TChiHH_HToBB_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16NanoAODv5_PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1_0.root
# Output:  SMS-TChiHH_HToBB_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv5__PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1_0.root
def convert_name(file_paths_string, split_strings):
  converted_name = os.path.basename(file_paths_string)
  for split_string in split_strings:
    converted_name = converted_name.replace(split_string, '_'+split_string)
  converted_name = converted_name.replace('_*','')

  return converted_name

if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser.add_argument('mass')
  parser.add_argument('input_paths')
  parser.add_argument('output_directory')
  args = vars(parser.parse_args())
  
  mass = args['mass']
  #file_paths_string = 'root_files/SMS-TChiHH_HToBB_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16NanoAODv5_PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1_*.root'
  file_paths_string = args['input_paths']

  print(file_paths_string)

  # Input:   SMS-TChiHH_HToBB_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16NanoAODv5_PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1_0.root
  # Output: SMS-TChiHH_mChi-1000_mLSP-1_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv5__PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1.root
  out_filename = convert_name(file_paths_string, ['RunII','PUSummer'])
  out_filename = out_filename.replace('SMS-TChiHH_HToBB_HToBB', 'SMS-TChiHH_mChi-'+str(mass)+'_mLSP-1')
  out_file_path = os.path.join(args['output_directory'], out_filename)
  
  chain = ROOT.TChain('Events')
  chain.Add(file_paths_string)
  cut_string = "GenPart_pdgId==1000023&&(GenPart_mass>"+str(int(mass)-10)+"&&GenPart_mass<"+str(int(mass)+10)+")"
  #cut_string = "GenPart_pdgId==1000023&&GenPart_mass=="+mass
  #print(cut_string)
  chain.Draw(">>elist",cut_string)
  elist = ROOT.gDirectory.Get("elist")
  chain.SetEventList(elist)
  out_file = ROOT.TFile(out_file_path, "recreate", "", 209)
  chain.CopyTree("");
  out_file.Write();

  #chain.Scan("GenPart_pdgId:GenPart_mass:GenPart_genPartIdxMother")

  #root_number_entries = chain.GetEntries()
  #print(root_number_entries)
