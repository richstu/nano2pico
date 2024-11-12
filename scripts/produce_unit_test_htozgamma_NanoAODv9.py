#!/bin/sh
''''exec python3 -u -- "$0" ${1+"$@"} # '''
import os
import sys
import subprocess
import threading
import time
import argparse

def output_reader(process, commandOutput):
  for line in iter(process.stdout.readline, b''):
    output = line.decode('utf-8')
    print(output, end='')
    commandOutput[0] += output

def runCommand(command):
  print("\n[Info] Running command: "+command)
  process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)

  commandOutput = ['']
  thread = threading.Thread(target=output_reader, args=(process,commandOutput))
  thread.start()
  thread.join()
  # Try to get poll
  for iTime in range(60):
    if (process.poll() == None): time.sleep(1)
    else: break
  return process.poll(), commandOutput[0]

def makeTestCommands(pico_directory, n_entries, mc_nanoaod_directory="", mc_nanoaod_filename="", data_nanoaod_directory="", data_nanoaod_filename="", year = ""):
  higgsino_test_commands = []
  if (mc_nanoaod_directory !=""):
    test_commands = [
      # mc commands
      "mkdir -p "+pico_directory+"/"+year+"/mc/raw_pico",
      "mkdir "+pico_directory+"/"+year+"/mc/wgt_sums",
      "./run/process_nano.exe -f "+mc_nanoaod_filename+" -i "+mc_nanoaod_directory+" -o "+pico_directory+"/"+year+"/mc --nent "+str(n_entries),
      "mkdir "+pico_directory+"/"+year+"/mc/corrections",
      "./run/merge_corrections.exe "+pico_directory+"/"+year+"/mc/corrections/"+mc_nanoaod_filename+" "+pico_directory+"/"+year+"/mc/wgt_sums/wgt_sums_"+mc_nanoaod_filename,
      "mkdir "+pico_directory+"/"+year+"/mc/unskimmed",
      "./run/apply_corrections.exe -f raw_pico_"+mc_nanoaod_filename+" -i "+pico_directory+"/"+year+"/mc/raw_pico/ -c "+mc_nanoaod_filename,
      "mkdir "+pico_directory+"/"+year+"/mc/llg",
      "./scripts/skim_file.py -k llg -i "+pico_directory+"/"+year+"/mc/unskimmed/pico_"+mc_nanoaod_filename+" -o "+pico_directory+"/"+year+"/mc/llg/",
      "mkdir "+pico_directory+"/"+year+"/mc/merged_llg",
      "./scripts/slim_and_merge.py -s txt/slim_rules/zgmc.txt -i "+pico_directory+"/"+year+"/mc/llg/pico_llg_"+mc_nanoaod_filename+" -o "+pico_directory+"/"+year+"/mc/merged_llg/merged_"+mc_nanoaod_filename,
    ]
    higgsino_test_commands.extend(test_commands)
  if (data_nanoaod_directory!=""):
    test_commands = [
      # data commands
      "mkdir -p "+pico_directory+"/"+year+"/data/raw_pico",
      "mkdir "+pico_directory+"/"+year+"/data/wgt_sums",
      "./run/process_nano.exe -f "+data_nanoaod_filename+" -i "+data_nanoaod_directory+" -o "+pico_directory+"/"+year+"/data --nent "+str(n_entries),
      "mkdir "+pico_directory+"/"+year+"/data/llg",
      "./scripts/skim_file.py -k llg -i "+pico_directory+"/"+year+"/data/raw_pico/raw_pico_"+data_nanoaod_filename+" -o "+pico_directory+"/"+year+"/data/llg/",
      "mkdir "+pico_directory+"/"+year+"/data/merged_llg",
      "./scripts/slim_and_merge.py -s txt/slim_rules/zgdata.txt -i "+pico_directory+"/"+year+"/data/llg/raw_pico_llg_"+data_nanoaod_filename+" -o "+pico_directory+"/"+year+"/data/merged_llg/merged_"+data_nanoaod_filename,
    ]
    higgsino_test_commands.extend(test_commands)
  return higgsino_test_commands

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='''Runs production on a few NanoAOD files.''', formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('-i','--input_pico_folder', default="/net/cms11/cms11r0/pico", help='Input pico folder.')
  parser.add_argument('-f','--output_folder', required=True, help='Output folder containing picos.')
  parser.add_argument('-l','--output_log', required=True, help='Output log filename.')
  args = parser.parse_args()

  if 'zgamma' not in args.output_folder:
    print("[Error] --output_folder needs to include zgamma in filename for zgamma production.")
    sys.exit(1)

  #pico_directory = "unit_test_htozgamma_nanoaodv9"
  #log_filename = "unit_test_htozgamma_nanoaodv9.log"
  input_pico_folder = args.input_pico_folder
  pico_directory = args.output_folder
  log_filename = args.output_log

  n_entries = 3000 # about 1 min for process_nano.exe (~550 Hz). 
  higgsino_test_commands = []
  # For DYJets and EG data
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016APV/mc", 
                                                 mc_nanoaod_filename="DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v1__280000__52E5237A-EB6F-4F43-A91F-20E2EAAD9E7D.root",
                                                 data_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016APV/data", data_nanoaod_filename="DoubleEG__Run2016B__ver1_HIPM_UL2016_MiniAODv2_NanoAODv9-v2__2530000__BD01422E-8128-834F-93FE-07E8F922A3D9.root", 
                                                 year = "2016APV",
                                                 pico_directory=pico_directory, n_entries=n_entries))
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016/mc", mc_nanoaod_filename="DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16NanoAODv9__106X_mcRun2_asymptotic_v17-v1__30000__0082C29D-E74C-024A-BE9B-97B29EE7A4A2.root", 
                                                 data_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016/data", data_nanoaod_filename="DoubleEG__Run2016H__UL2016_MiniAODv2_NanoAODv9-v1__260000__FC1ECF3C-83C1-394B-A69B-9BFB47A1BA61.root", 
                                                 year = "2016",
                                                 pico_directory=pico_directory, n_entries=n_entries))
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2017/mc", mc_nanoaod_filename="DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v2__100000__04620FA2-DB18-7F4A-B30C-BAFD1C5B673D.root", 
                                                 data_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2017/data", data_nanoaod_filename="DoubleEG__Run2017B__UL2017_MiniAODv2_NanoAODv9-v1__70000__04646CDD-F24C-DA40-BCC8-6EB722486EAF.root", 
                                                 year = "2017",
                                                 pico_directory=pico_directory, n_entries=n_entries))
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2018/mc", mc_nanoaod_filename="DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018_realistic_v16_L1v1-v2__230000__00EA9563-5449-D24E-9566-98AE8E2A61AE.root", 
                                                 data_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2018/data", 
                                                 data_nanoaod_filename="EGamma__Run2018A__UL2018_MiniAODv2_NanoAODv9-v1__270000__00B7FFB1-3455-C941-AE3B-CF7085966A41.root", 
                                                 year = "2018",
                                                 pico_directory=pico_directory, n_entries=n_entries))
                                                
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv12/nano/2022/mc", 
                                                 mc_nanoaod_filename="DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8__Run3Summer22NanoAODv12__130X_mcRun3_2022_realistic_v5_ext1-v1__2560000__74ed7194-9623-484f-adf4-d55268f95320.root", 
                                                 data_nanoaod_directory=input_pico_folder+"/NanoAODv12/nano/2022/data", 
                                                 data_nanoaod_filename="Muon__Run2022D__22Sep2023-v1__2520000__2fe0b89f-4fef-4cd0-b81a-eb5c8964194e.root", 
                                                 year = "2022",
                                                 pico_directory=pico_directory, n_entries=n_entries))

  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv12/nano/2022EE/mc", 
                                                 mc_nanoaod_filename="DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8__Run3Summer22EENanoAODv12__forPOG_130X_mcRun3_2022_realistic_postEE_v6-v2__2530000__94a3133e-238e-4b52-8862-1e5bf11cf441.root", 
                                                 data_nanoaod_directory=input_pico_folder+"/NanoAODv12/nano/2022EE/data", 
                                                 data_nanoaod_filename="EGamma__Run2022E__22Sep2023-v1__2530000__0900e874-435c-4ac1-aa69-dd6a5eff46c7.root", 
                                                 year = "2022EE",
                                                 pico_directory=pico_directory, n_entries=n_entries))

  higgsino_test_commands.extend(makeTestCommands(data_nanoaod_directory=input_pico_folder+"/NanoAODv12/nano/2023/data", 
                                                 data_nanoaod_filename="Muon0__Run2023C__22Sep2023_v4-v1__50000__b62545f1-6761-4a06-8d2b-9b02e10685a7.root", 
                                                 year = "2023",
                                                 pico_directory=pico_directory, n_entries=n_entries))

  higgsino_test_commands.extend(makeTestCommands(data_nanoaod_directory=input_pico_folder+"/NanoAODv12/nano/2023BPix/data", 
                                                 data_nanoaod_filename="Muon0__Run2023D__22Sep2023_v1-v1__2530000__45788039-83df-45b9-9c6c-490f2fd73c08.root", 
                                                 year = "2023BPix",
                                                 pico_directory=pico_directory, n_entries=n_entries))
  # For ggH
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016/mc", mc_nanoaod_filename="GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL16NanoAODv9__106X_mcRun2_asymptotic_v17-v1__40000__C3E0FC8E-4157-7E40-956C-F04E32F6B8C7.root", 
                                                 year = "2016",
                                                 pico_directory=pico_directory, n_entries=n_entries))
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016APV/mc", 
                                                 mc_nanoaod_filename="GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v1__40000__FF1A0E0D-369B-FD43-87A0-5E058074F250.root",
                                                 year = "2016APV",
                                                 pico_directory=pico_directory, n_entries=n_entries))
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2017/mc", mc_nanoaod_filename="GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1__70000__93B37B9F-16D5-6F4F-B920-3B9C682CA8A8.root", 
                                                 year = "2017",
                                                 pico_directory=pico_directory, n_entries=n_entries))
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2018/mc", mc_nanoaod_filename="GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018_realistic_v16_L1v1-v1__2550000__32457BD6-6802-DF4C-A62E-75FB77525DB6.root", 
                                                 year = "2018",
                                                 pico_directory=pico_directory, n_entries=n_entries))
  # For Zg
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016/mc", mc_nanoaod_filename="ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16NanoAODv9__106X_mcRun2_asymptotic_v17-v1__280000__6FB79EFC-1B50-F94F-BED8-34BD7935F8CF.root", 
                                                 year = "2016",
                                                 pico_directory=pico_directory, n_entries=n_entries))
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016APV/mc", 
                                                 mc_nanoaod_filename="ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v1__40000__CF762017-6761-1B44-BEE1-FE1B08B6F8A7.root",
                                                 year = "2016APV",
                                                 pico_directory=pico_directory, n_entries=n_entries))
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2017/mc", mc_nanoaod_filename="ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1__70000__EA57A9D4-BB7D-364E-9AE5-C7C6CC0E6E4E.root", 
                                                 year = "2017",
                                                 pico_directory=pico_directory, n_entries=n_entries))
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2018/mc", mc_nanoaod_filename="ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018_realistic_v16_L1v1-v1__70000__783A2FC8-6CAB-CD46-A7D2-F0AA56E72197.root", 
                                                 year = "2018",
                                                 pico_directory=pico_directory, n_entries=n_entries))
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9UCSB1/nano/2018/mc", mc_nanoaod_filename="ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL18NanoAODv9UCSB1__106X_upgrade2018_realistic_v16_L1v1-v1__028.root", 
                                                 year = "2018",
                                                 pico_directory=pico_directory, n_entries=n_entries))

  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv12/nano/2022/mc", 
                                                 mc_nanoaod_filename="DYGto2LG-1Jets_MLL-50_PTG-50to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8__Run3Summer22NanoAODv12__130X_mcRun3_2022_realistic_v5-v3__2560000__c894c764-6dcc-4723-a4d8-456f058e215e.root", 
                                                 year = "2022",
                                                 pico_directory=pico_directory, n_entries=n_entries))

  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv12/nano/2022/mc", 
                                                 mc_nanoaod_filename="GluGluHtoZG_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8__Run3Summer22NanoAODv12__130X_mcRun3_2022_realistic_v5-v2__2530000__6e37d56d-c21d-4d3d-b647-8e941c6fa844.root", 
                                                 year = "2022",
                                                 pico_directory=pico_directory, n_entries=n_entries))


  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv12/nano/2022EE/mc", 
                                                 mc_nanoaod_filename="DYGto2LG-1Jets_MLL-50_PTG-50to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8__Run3Summer22EENanoAODv12__130X_mcRun3_2022_realistic_postEE_v6-v1__40000__9eaff41b-6b12-463b-97a1-b56f798821c5.root", 
                                                 year = "2022EE",
                                                 pico_directory=pico_directory, n_entries=n_entries))

  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv12/nano/2022EE/mc", 
                                                 mc_nanoaod_filename="GluGluHtoZG_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8__Run3Summer22EENanoAODv12__130X_mcRun3_2022_realistic_postEE_v6-v2__2530000__4d0c2375-cd02-4dd3-a4e6-cd2648f3a7e3.root", 
                                                 year = "2022EE",
                                                 pico_directory=pico_directory, n_entries=n_entries))


  os.makedirs(pico_directory)
  log_file = open(pico_directory+"/"+log_filename, 'w')
  for command in higgsino_test_commands:
    t0 = time.time()
    return_code, output = runCommand(command)
    execution_time = time.time()-t0
    log_file.write("[Info] command: "+command+" return code: "+str(return_code)+" execution time: "+"{:.1f}".format(execution_time)+" seconds\n")
    log_file.write(output)
    if "mkdir" not in command:
      if return_code != 0: 
        sys.exit("[Error] Command: ("+command+") has an error")
  log_file.close()
