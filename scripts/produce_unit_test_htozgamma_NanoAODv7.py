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

  #pico_directory = "unit_test_htozgamma_nanoaodv7"
  #log_filename = "unit_test_htozgamma_nanoaodv7.log"
  pico_directory = args.output_folder
  log_filename = args.output_log

  n_entries = 100000 # about 3 min for process_nano.exe (~550 Hz). 
  higgsino_test_commands = []
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory="/net/cms17/cms17r0/pico/NanoAODv7/nano/2016/mc", mc_nanoaod_filename="DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv7__PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1__260000__017F5826-2321-C34C-9525-2EE901E65759.root", 
                                                 data_nanoaod_directory="/net/cms17/cms17r0/pico/NanoAODv7/nano/2016/data", data_nanoaod_filename="DoubleEG__Run2016B__02Apr2020_ver2-v1__100000__03757D92-2C17-6640-A535-D2C6E77A7ECC.root", 
                                                 pico_directory=pico_directory, n_entries=n_entries))
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory="/net/cms17/cms17r0/pico/NanoAODv7/nano/2017/mc", mc_nanoaod_filename="DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIIFall17NanoAODv7__PU2017RECOSIMstep_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8_ext1-v1__100000__0FD3A2B8-9BD9-2442-83C7-B43CC073443F.root", 
                                                 data_nanoaod_directory="/net/cms17/cms17r0/pico/NanoAODv7/nano/2017/data", data_nanoaod_filename="DoubleEG__Run2017B__02Apr2020-v1__30000__08F99D33-012C-1641-BB33-1204571842F9.root", 
                                                 pico_directory=pico_directory, n_entries=n_entries))
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory="/net/cms17/cms17r0/pico/NanoAODv7/nano/2018/mc", mc_nanoaod_filename="DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIIAutumn18NanoAODv7__Nano02Apr2020_102X_upgrade2018_realistic_v21-v1__100000__3391856D-092B-624A-A3BB-F6112936F5D9.root", 
                                                 data_nanoaod_directory="/net/cms17/cms17r0/pico/NanoAODv7/nano/2018/data", data_nanoaod_filename="EGamma__Run2018A__02Apr2020-v1__2410000__1D1F868C-0C62-E447-B2D7-2C37CC7EAAD1.root", 
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
