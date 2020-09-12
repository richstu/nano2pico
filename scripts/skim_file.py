#!/usr/bin/env python
import os
import ROOT
import argparse

# Input:   root_files/SMS-TChiHH_HToBB_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16NanoAODv5_PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1_0.root
# Output:  SMS-TChiHH_HToBB_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv5__PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1_0.root
def convert_name(file_paths_string, split_strings):
  converted_name = os.path.basename(file_paths_string)
  # Below is not needed anymore.
  #for split_string in split_strings:
  #  converted_name = converted_name.replace(split_string, '_'+split_string)
  converted_name = converted_name.replace('_*','')

  return converted_name

def get_cuts(skim_name):
  cuts = ''

  # General use
  pass_1l_trig40 = '(Max$(el_pt*el_sig)>40 || Max$(mu_pt*mu_sig)>40)' # use for 1L CR
  pass_1l_trig30 = '(Max$(el_pt*el_sig)>30 || Max$(mu_pt*mu_sig)>30)' # use for 2L CR, can lower the cut since two leps!
  mllcut = '@ll_m.size()>=1 && Sum$(ll_m>80 && ll_m<100)>=1'
  if(skim_name=='met150'): cuts = 'met>150'
  if(skim_name=='nj45'): cuts = 'njet>=4 && njet<=5'
  if(skim_name=='zcand'): cuts = '&&'.join(['nlep==2', 'nbm==0', pass_1l_trig30, mllcut])
  if(skim_name=='ttisr'): cuts = '&&'.join(['nlep==2', 'nbm==2', pass_1l_trig30])
  if(skim_name=='wisr'):  cuts = '&&'.join(['met>100', 'nbl==0', pass_1l_trig40])
  if(skim_name=='1l2j'): cuts = '&&'.join(['nlep>=1', 'njet>=2'])
  if(skim_name=='1e4j'): cuts = '&&'.join(['nel>=1', 'njet>=4'])
  if(skim_name=='2l'): cuts = 'nlep>=2'
  if(skim_name=='3l'): cuts = 'nlep==3'
  if(skim_name=='llg'): cuts = '(nmu > 1 || nel > 1) && nphoton > 0'

  # Higgsino loose
  nb_or_fjet_cut = '(nbt>=2 || nbdft>=2 || Sum$(fjet_pt>300 && fjet_msoftdrop>50)>0)'
  if(skim_name=='higloose'): cuts = '&&'.join([nb_or_fjet_cut, 'met>150', 'nvlep==0'])
  
  # Higgsino tight
  higtrim = '((Alt$(hig_cand_drmax[0],0)<2.2 && Alt$(hig_cand_dm[0],0)<=40 && Alt$(hig_cand_am[0],0)<=200) ||'
  higtrim += '(Alt$(hig_df_cand_drmax[0],0)<2.2 && Alt$(hig_df_cand_dm[0],0)<=40 && Alt$(hig_df_cand_am[0],0)<=200))'
  resolved = 'nbt>=2 && njet>=4 && njet<=5 &&' + higtrim
  boosted = 'Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1'
  if(skim_name=='higtight'): 
    cuts = '&&'.join(['nvlep==0', 'ntk==0', 'met>150', '!low_dphi_met', '(('+resolved+')||('+boosted+'))'])
    print('Using cut string:  '+cuts.replace('&&',' && '))

  # Higgsino preselection for training DNN, can be done from higloose
  resolved = 'nbt>=2 && njet>=4 && njet<=5'
  boosted = 'Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1'
  if(skim_name=='preselect'): 
    cuts = '&&'.join(['nvlep==0', 'ntk==0', '!low_dphi_met', 'met>150', '(('+resolved+')||('+boosted+'))'])
    print('Using cut string:  '+cuts.replace('&&',' && '))

  # Control regions skims - to be updated when needed
  nbnjet_or_fjet_cut = '(Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || ((nbt>=2 || nbdft>=2) && njet>=4 && njet<=5))' 
  if(skim_name=='higlep1T'):  cuts = '&&'.join([nbnjet_or_fjet_cut, 'nlep==1', pass_1l_trig40])
  # do not cut on Nb since CR goes down to 0b!
  njet_or_fjet_cut = '(Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || (njet>=4 && njet<=5))' 
  if(skim_name=='higlep2T'):  cuts = '&&'.join([njet_or_fjet_cut, 'nlep==2', mllcut,  pass_1l_trig30])

  if(skim_name=='higqcd'):  cuts = '&&'.join([njet_or_fjet_cut, 'nvlep==0', 'ntk==0', 'low_dphi_met'])

  # Loosen up just enough to do systematics - to be updated when needed
  # sys_nbcut = 'max(nbdft,Max$(sys_nbdft))>=2'
  # sys_njcut = '(njet==4||sys_njet[1]==4||sys_njet[2]==4||njet==5||sys_njet[1]==5||sys_njet[2]==5)'
  # sys_higtrim = '&&'.join(['min(hig_cand_drmax,Min$(sys_hig_cand_drmax))<2.2',
  #                          'min(hig_cand_dm,Min$(sys_hig_cand_dm))<=40',
  #                          'min(hig_cand_am,Min$(sys_hig_cand_am))<=200'])
  # if(skim_name=='higsys'):   cuts = '&&'.join([sys_njcut, sys_nbcut, 'max(met,Max$(sys_met))>150', 'nvlep==0', 'ntk==0', sys_higtrim])

  return cuts

if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser.add_argument('-m','--mass', default = -1)
  parser.add_argument('-l','--mass_lsp', default = 0)
  parser.add_argument('-k','--skim_name', default = "")
  parser.add_argument('-i','--input_paths', required=True, default = "")
  parser.add_argument('-o','--output_dir', required=True, default = "")
  args = vars(parser.parse_args())
  
  file_paths_string = args['input_paths']
  cut_string = ''
  out_file_path = ''
  treename = ''

  if args['mass'] !=-1:
    mass = args['mass']
    mass_lsp = args['mass_lsp']
    #cut_string = "GenPart_pdgId==1000023&&(GenPart_mass>"+str(int(mass)-10)+"&&GenPart_mass<"+str(int(mass)+10)+")"  
    cut_string = "(MaxIf$(GenPart_mass,GenPart_pdgId==1000023)>"+str(int(mass)-10)+"&&MaxIf$(GenPart_mass,GenPart_pdgId==1000023)<"+str(int(mass)+10)+")"
    cut_string += "&&(MaxIf$(GenPart_mass,GenPart_pdgId==1000022)>"+str(int(mass_lsp)-10)+"&&MaxIf$(GenPart_mass,GenPart_pdgId==1000022)<"+str(int(mass_lsp)+10)+")"
    # Input:   SMS-TChiHH_HToBB_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16NanoAODv5_PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1_0.root
    # Output: SMS-TChiHH_mChi-1000_mLSP-1_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv5__PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1.root
    out_filename = convert_name(file_paths_string, ['RunII','PUSummer'])
    #out_filename = out_filename.replace('SMS-TChiHH_HToBB_HToBB', 'SMS-TChiHH_mChi-'+str(mass)+'_mLSP-1')
    out_filename = out_filename.replace('SMS-TChiHH', 'SMS-TChiHH_mChi-'+str(mass)+'_mLSP-'+str(mass_lsp))
    out_file_path = os.path.join(args['output_dir'], out_filename)
    treename = 'Events'
  elif args['skim_name'] !="":
    cut_string = get_cuts(args['skim_name'])
    out_filename = os.path.basename(os.path.realpath(file_paths_string))
    out_file_path = os.path.join(args['output_dir'], out_filename.replace('pico_','pico_'+args['skim_name']+'_'))
    treename = 'tree'
  else:
    sys.exit("You have to specify either mass or cut string.")
  

  chain = ROOT.TChain(treename)
  chain.Add(file_paths_string)
  nent_total = chain.GetEntries()
  nent_skim = chain.GetEntries(cut_string)
  print('Found {} events satisfying the skim requirements.'.format(nent_skim))

  out_file = ROOT.TFile(out_file_path, "recreate") # in case we want compression 
  if (args['mass']!=-1): 
    out_file.SetCompressionSettings(209); # kLZMA algo with level 9 compression

  if nent_total>0:
    chain.Draw(">>elist",cut_string)
    elist = ROOT.gDirectory.Get("elist")
    chain.SetEventList(elist)
  chain.CopyTree("");
  out_file.Write();

  print('Wrote file '+out_file_path)
