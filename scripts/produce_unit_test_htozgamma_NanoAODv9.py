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

def makeTestCommands(mc_nanoaod_directory, mc_nanoaod_filename, data_nanoaod_directory, data_nanoaod_filename, pico_directory, n_entries):
  higgsino_test_commands = [
    # mc commands
    "mkdir -p "+pico_directory+"/mc/raw_pico",
    "mkdir "+pico_directory+"/mc/wgt_sums",
    "./run/process_nano.exe -f "+mc_nanoaod_filename+" -i "+mc_nanoaod_directory+" -o "+pico_directory+"/mc --nent "+str(n_entries),
    "mkdir "+pico_directory+"/mc/corrections",
    "./run/merge_corrections.exe "+pico_directory+"/mc/corrections/"+mc_nanoaod_filename+" "+pico_directory+"/mc/wgt_sums/wgt_sums_"+mc_nanoaod_filename,
    "mkdir "+pico_directory+"/mc/unskimmed",
    "./run/apply_corrections.exe -f raw_pico_"+mc_nanoaod_filename+" -i "+pico_directory+"/mc/raw_pico/ -c "+mc_nanoaod_filename,
    "mkdir "+pico_directory+"/mc/llg",
    "./scripts/skim_file.py -k llg -i "+pico_directory+"/mc/unskimmed/pico_"+mc_nanoaod_filename+" -o "+pico_directory+"/mc/llg/",
    "mkdir "+pico_directory+"/mc/merged_llg",
    "./scripts/slim_and_merge.py -s txt/slim_rules/zgmc.txt -i "+pico_directory+"/mc/llg/pico_llg_"+mc_nanoaod_filename+" -o "+pico_directory+"/mc/merged_llg/merged_"+mc_nanoaod_filename,

    # data commands
    "mkdir -p "+pico_directory+"/data/raw_pico",
    "mkdir "+pico_directory+"/data/wgt_sums",
    "./run/process_nano.exe -f "+data_nanoaod_filename+" -i "+data_nanoaod_directory+" -o "+pico_directory+"/data --nent "+str(n_entries),
    "mkdir "+pico_directory+"/data/llg",
    "./scripts/skim_file.py -k llg -i "+pico_directory+"/data/raw_pico/raw_pico_"+data_nanoaod_filename+" -o "+pico_directory+"/data/llg/",
    "mkdir "+pico_directory+"/data/merged_llg",
    "./scripts/slim_and_merge.py -s txt/slim_rules/zgdata.txt -i "+pico_directory+"/data/llg/raw_pico_llg_"+data_nanoaod_filename+" -o "+pico_directory+"/data/merged_llg/merged_"+data_nanoaod_filename,

  ]
  return higgsino_test_commands

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='''Runs production on a few NanoAOD files.''', formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('-f','--output_folder', required=True, help='Output folder containing picos.')
  parser.add_argument('-l','--output_log', required=True, help='Output log filename.')
  args = parser.parse_args()

  if 'zgamma' not in args.output_folder:
    print("[Error] --output_folder needs to include zgamma in filename for zgamma production.")
    sys.exit(1)

  #pico_directory = "unit_test_htozgamma_nanoaodv9"
  #log_filename = "unit_test_htozgamma_nanoaodv9.log"
  pico_directory = args.output_folder
  log_filename = args.output_log

  n_entries = 100000 # about 3 min for process_nano.exe (~550 Hz). 
  higgsino_test_commands = []
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory="/net/cms17/cms17r0/pico/NanoAODv9/nano/2016/mc", mc_nanoaod_filename="DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16NanoAODv9__106X_mcRun2_asymptotic_v17-v1__30000__0082C29D-E74C-024A-BE9B-97B29EE7A4A2.root", 
                                                 data_nanoaod_directory="/net/cms17/cms17r0/pico/NanoAODv9/nano/2016/data", data_nanoaod_filename="DoubleEG__Run2016B__ver1_HIPM_UL2016_MiniAODv2_NanoAODv9-v1__280000__0475EF51-B32B-3842-87A1-4401BB753DD9.root", 
                                                 pico_directory=pico_directory, n_entries=n_entries))
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory="/net/cms17/cms17r0/pico/NanoAODv9/nano/2017/mc", mc_nanoaod_filename="DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v2__100000__04620FA2-DB18-7F4A-B30C-BAFD1C5B673D.root", 
                                                 data_nanoaod_directory="/net/cms17/cms17r0/pico/NanoAODv9/nano/2017/data", data_nanoaod_filename="DoubleEG__Run2017B__UL2017_MiniAODv2_NanoAODv9-v1__70000__04646CDD-F24C-DA40-BCC8-6EB722486EAF.root", 
                                                 pico_directory=pico_directory, n_entries=n_entries))
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory="/net/cms17/cms17r0/pico/NanoAODv9/nano/2018/mc", mc_nanoaod_filename="DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018_realistic_v16_L1v1-v2__230000__00EA9563-5449-D24E-9566-98AE8E2A61AE.root", 
                                                 data_nanoaod_directory="/net/cms17/cms17r0/pico/NanoAODv9/nano/2018/data", data_nanoaod_filename="EGamma__Run2018A__UL2018_MiniAODv2_NanoAODv9-v1__270000__00B7FFB1-3455-C941-AE3B-CF7085966A41.root", 
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
