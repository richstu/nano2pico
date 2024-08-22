#!/usr/bin/env python
# The ./jobscript_check.py should return 'success' or 'fail' or 'to_submit' or 'submitted' for a job_log_string
# The input is given as sys.argv[1] = queue_system.compress_string(job_log_string) sys.argv[2] = queue_system.compress_string(job_argument_string)
import sys
import os
import queue_system
from ROOT import TChain

# job_argument_string = "/net/top/homes/aovcharova/code/nano2pico/run/process_nano.exe -f DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv5__PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1__250000__233E4358-D599-FE4A-B585-A6B18F4DDEF1.root -i /mnt/hadoop/pico/NanoAODv5/nano/2016/mc/ -o /net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/"
job_log_string = queue_system.decompress_string(sys.argv[1])
job_argument_string = queue_system.decompress_string(sys.argv[2])
#job_argument_string = "--command=\"/net/cms37/data1/jbkim/analysis/nano2pico.inyo/run/process_nano.exe -f JetHT__Run2016H__02Apr2020-v1__40000__E316083A-DB8C-484C-8013-C9C3A301ED61.root -i /net/cms25/cms25r5/pico/NanoAODv7/nano/2016/data/ -o /net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/2016/data/\""

print(job_argument_string)

class GoldenJson:

  #dictionary where key is the run # and value is a list of tuples of start and end points of good lumi blocks
  good_lumi_blocks = {}

  def __init__(self):
    #load golden jsons
    for json_filename in os.listdir('txt/json/'):
      json_file = open('txt/json/'+json_filename,'r')
      json_file_text = json_file.read().split('\n')
      json_file.close()
      for run_line in json_file_text:
        if (len(run_line.split('"')) < 2):
          continue
        run_number = int(run_line.split('"')[1])
        lumi_blocks_string = run_line[run_line.find('['):run_line.rfind(']')]
        lumi_blocks_list = lumi_blocks_string.replace('[','').replace(']','').replace(' ','').split(',')
        int_lumi_blocks_list = []
        for index in range(len(lumi_blocks_list)/2):
          int_lumi_blocks_list.append((int(lumi_blocks_list[2*index]),int(lumi_blocks_list[2*index+1])))
        self.good_lumi_blocks[run_number] = int_lumi_blocks_list

  def check(self, run, lumi_block):
    if run in self.good_lumi_blocks:
      for lumi_range in self.good_lumi_blocks[run]:
        if (lumi_block >= lumi_range[0] and lumi_block <= lumi_range[1]):
          return True
    return False

def CheckPassTriggers(chain, trigger_list):
  r_check_pass_triggers = False
  for trigger_name in trigger_list:
    if hasattr(chain, trigger_name):
      r_check_pass_triggers = r_check_pass_triggers or getattr(chain, trigger_name)
  return r_check_pass_triggers

golden_json = GoldenJson()
met_triggers = ['HLT_PFMET90_PFMHT90_IDTight','HLT_PFMETNoMu90_PFMHTNoMu90_IDTight','HLT_PFMET100_PFMHT100_IDTight','HLT_PFMETNoMu100_PFMHTNoMu100_IDTight','HLT_PFMET110_PFMHT110_IDTight','HLT_PFMETNoMu110_PFMHTNoMu110_IDTight','HLT_PFMET120_PFMHT120_IDTight','HLT_PFMETNoMu120_PFMHTNoMu120_IDTight','HLT_PFMET130_PFMHT130_IDTight','HLT_PFMETNoMu130_PFMHTNoMu130_IDTight','HLT_PFMET140_PFMHT140_IDTight','HLT_PFMETNoMu140_PFMHTNoMu140_IDTight','HLT_PFMET100_PFMHT100_IDTight_PFHT60','HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60','HLT_PFMET110_PFMHT110_IDTight_PFHT60','HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60','HLT_PFMET120_PFMHT120_IDTight_PFHT60','HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60','HLT_PFMET130_PFMHT130_IDTight_PFHT60','HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60','HLT_PFMET140_PFMHT140_IDTight_PFHT60','HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_PFHT60','HLT_PFMET120_PFMHT120_IDTight','HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned','HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned','HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1','HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1','HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1']
egamma_triggers = ['HLT_Ele25_WPTight_Gsf','HLT_Ele27_WPTight_Gsf','HLT_Ele28_WPTight_Gsf','HLT_Ele32_WPTight_Gsf','HLT_Ele32_WPTight_Gsf_L1DoubleEG','HLT_Ele35_WPTight_Gsf','HLT_Ele20_WPLoose_Gsf','HLT_Ele45_WPLoose_Gsf','HLT_Ele105_CaloIdVT_GsfTrkIdT','HLT_Ele115_CaloIdVT_GsfTrkIdT','HLT_Ele135_CaloIdVT_GsfTrkIdT','HLT_Ele145_CaloIdVT_GsfTrkIdT','HLT_Ele25_eta2p1_WPTight_Gsf','HLT_Ele27_eta2p1_WPTight_Gsf','HLT_Ele20_eta2p1_WPTight_Gsf','HLT_Ele25_eta2p1_WPTight_Gsf','HLT_Ele27_eta2p1_WPTight_Gsf','HLT_Ele15_IsoVVVL_PFHT350','HLT_Ele15_IsoVVVL_PFHT400','HLT_Ele15_IsoVVVL_PFHT450','HLT_Ele15_IsoVVVL_PFHT600','HLT_Ele50_IsoVVVL_PFHT450']
muon_triggers = ['HLT_IsoMu20','HLT_IsoMu22','HLT_IsoMu24','HLT_IsoMu27','HLT_IsoTkMu20','HLT_IsoTkMu22','HLT_IsoTkMu24','HLT_Mu50','HLT_Mu55','HLT_TkMu50','HLT_IsoMu22_eta2p1','HLT_IsoMu24_eta2p1','HLT_Mu45_eta2p1','HLT_Mu15_IsoVVVL_PFHT350','HLT_Mu15_IsoVVVL_PFHT400','HLT_Mu15_IsoVVVL_PFHT450','HLT_Mu15_IsoVVVL_PFHT600','HLT_Mu50_IsoVVVL_PFHT400','HLT_Mu50_IsoVVVL_PFHT450']
jetht_triggers = ['HLT_PFJet500','HLT_PFHT125','HLT_PFHT200','HLT_PFHT300','HLT_PFHT400','HLT_PFHT475','HLT_PFHT600','HLT_PFHT650','HLT_PFHT800','HLT_PFHT900','HLT_PFHT180','HLT_PFHT370','HLT_PFHT430','HLT_PFHT510','HLT_PFHT590','HLT_PFHT680','HLT_PFHT780','HLT_PFHT890','HLT_PFHT1050','HLT_PFHT250','HLT_PFHT350']

args = job_argument_string.split('--command="')[1].split('"')[0]
tmp = args.split(' ')
infile_path = tmp[4]+'/'+tmp[2]
outfile_path = tmp[6]+'/raw_pico/raw_pico_'+tmp[2]

infile = TChain("Events");
infile.Add(infile_path);
in_nent = 0
if 'data/' in infile_path: # changed 'data' to 'data/' so that nano2pico can be run locally without errors 
  for i in range(0, infile.GetEntries()):
    infile.GetEntry(i)
    #check triggers, matching overlap removal scheme in event_tools.cpp
    #currently fails, need to have a way to check for branches in TChain
    pass_met_trigger = CheckPassTriggers(infile, met_triggers)
    pass_egamma_trigger = CheckPassTriggers(infile, egamma_triggers)
    pass_muon_trigger = CheckPassTriggers(infile, muon_triggers)
    pass_jetht_trigger = CheckPassTriggers(infile, jetht_triggers)

    if ('MET' in infile_path) and (not pass_met_trigger):
      continue
    if (('SingleElectron' in infile_path) or ('EGamma' in infile_path)) and (not pass_egamma_triggers or pass_met_trigger):
      continue
    if (('SingleMuon' in infile_path)) and (not pass_muon_trigger or pass_egamma_trigger or pass_met_trigger):
      continue
    if (('JetHT' in infile_path)) and (not pass_jetht_trigger or pass_muon_trigger or pass_egamma_trigger or pass_met_trigger):
      continue
    #check golden json
    if golden_json.check(infile.run, infile.luminosityBlock):
      in_nent = in_nent + 1
else:
  in_nent = infile.GetEntries()

outfile = TChain("tree");
outfile.Add(outfile_path);
out_nent = outfile.GetEntries()

print('DEBUG: input entries: {} output entries: {}'.format(in_nent, out_nent))

if infile.GetNbranches() == 0:
  print('[For queue_system] fail: input ({}) has no branches.'.format(infile_path))
if outfile.GetNbranches() == 0:
  print('[For queue_system] fail: output ({}) has no branches.'.format(outfile_path))
elif in_nent == out_nent:
  print('[For queue_system] success')
else:
  print('[For queue_system] fail: Input ({}) has {} entries, while output ({}) has {} entries.'.format(infile_path, in_nent, outfile_path, out_nent))
