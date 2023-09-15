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

# SampleTypeFlag is:
# 0: for MC 
# 1: for Signal
# 2: for Data

def makeTestCommands(nanoaod_directory, nanoaod_filename, SampleTypeFlag, pico_directory, n_entries): 
    
    if SampleTypeFlag == 0:

      higgsino_test_commands = [

      # mc commands
      "mkdir -p "+pico_directory+"/mc/raw_pico",
      "mkdir "+pico_directory+"/mc/wgt_sums",
      "./run/process_nano.exe -f "+nanoaod_filename+" -i "+nanoaod_directory+" -o "+pico_directory+"/mc --nent "+str(n_entries),
      "mkdir "+pico_directory+"/mc/corrections",
      "./run/merge_corrections.exe "+pico_directory+"/mc/corrections/"+nanoaod_filename+" "+pico_directory+"/mc/wgt_sums/wgt_sums_"+nanoaod_filename,
      "mkdir "+pico_directory+"/mc/unskimmed",
      "./run/apply_corrections.exe -f raw_pico_"+nanoaod_filename+" -i "+pico_directory+"/mc/raw_pico/ -c "+nanoaod_filename,
      "mkdir "+pico_directory+"/mc/met150",
      "./scripts/skim_file.py -k met150 -i "+pico_directory+"/mc/unskimmed/pico_"+nanoaod_filename+" -o "+pico_directory+"/mc/met150/",
      "mkdir "+pico_directory+"/mc/merged_met150",
      "./scripts/slim_and_merge.py -s txt/slim_rules/higmc.txt -i "+pico_directory+"/mc/met150/pico_met150_"+nanoaod_filename+" -o "+pico_directory+"/mc/merged_met150/merged_"+nanoaod_filename
      
      ]

    elif SampleTypeFlag == 1:

      higgsino_test_commands = [

      # signal commands
      "mkdir -p "+pico_directory+"/signal/raw_pico",
      "mkdir "+pico_directory+"/signal/wgt_sums",
      "./run/process_nano.exe -f "+nanoaod_filename+" -i "+nanoaod_directory+" -o "+pico_directory+"/signal",
      "mkdir "+pico_directory+"/signal/corrections",
      "./run/merge_corrections.exe "+pico_directory+"/signal/corrections/"+nanoaod_filename+" "+pico_directory+"/signal/wgt_sums/wgt_sums_"+nanoaod_filename,
      "mkdir "+pico_directory+"/signal/unskimmed",
      "./run/apply_corrections.exe -f raw_pico_"+nanoaod_filename+" -i "+pico_directory+"/signal/raw_pico/ -c "+nanoaod_filename,
      "mkdir "+pico_directory+"/signal/met150",
      "./scripts/skim_file.py -k met150 -i "+pico_directory+"/signal/unskimmed/pico_"+nanoaod_filename+" -o "+pico_directory+"/signal/met150/",
      "mkdir "+pico_directory+"/signal/merged_met150",
      "./scripts/slim_and_merge.py -s txt/slim_rules/higmc.txt -i "+pico_directory+"/signal/met150/pico_met150_"+nanoaod_filename+" -o "+pico_directory+"/signal/merged_met150/merged_"+nanoaod_filename

      ]

    elif SampleTypeFlag == 2:
      
      higgsino_test_commands = [

      # data commands
      "mkdir -p "+pico_directory+"/data/raw_pico",
      "mkdir "+pico_directory+"/data/wgt_sums",
      "./run/process_nano.exe -f "+nanoaod_filename+" -i "+nanoaod_directory+" -o "+pico_directory+"/data --nent "+str(n_entries),
      "mkdir "+pico_directory+"/data/met150",
      "./scripts/skim_file.py -k met150 -i "+pico_directory+"/data/raw_pico/raw_pico_"+nanoaod_filename+" -o "+pico_directory+"/data/met150/",
      "mkdir "+pico_directory+"/data/merged_met150",
      "./scripts/slim_and_merge.py -s txt/slim_rules/higdata.txt -i "+pico_directory+"/data/met150/raw_pico_met150_"+nanoaod_filename+" -o "+pico_directory+"/data/merged_met150/merged_"+nanoaod_filename

      ]
  
  return higgsino_test_commands

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='''Runs production on a few NanoAOD files.''', formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('-f','--output_folder', required=True, help='Ouptut folder containing picos.')
  parser.add_argument('-l','--output_log', required=True, help='Ouptut log filename.')
  args = parser.parse_args()

  #pico_directory = "unit_test_higgsino"
  #log_filename = "unit_test_higgsino_production.log"
  pico_directory = args.output_folder
  log_filename = args.output_log

  n_entries = 100000 # about 3 min for process_nano.exe (~550 Hz). Note for signal, nent is -1.
  higgsino_test_commands = []
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory="/net/cms37/data1/mhussain/HH-MET/DataSample/Background/DiPhotonJetsBox_MGG-80toInf", mc_nanoaod_filename="*", 
                                                 SampleTypeFlag = 0, pico_directory=pico_directory, n_entries=n_entries))
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory="/net/cms17/cms17r0/pico/NanoAODv7/nano/2017/mc", mc_nanoaod_filename="TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8__RunIIFall17NanoAODv7__PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1__100000__0ED19AF0-B248-8344-91D7-B241CE0729FA.root", 
                                                 signal_nanoaod_directory="/net/cms24/cms24r0/pico/NanoAODv7/nano/2017/SMS-TChiHH_2D_fastSimJmeCorrection", signal_nanoaod_filename="SMS-TChiHH_mChi-500_mLSP-0_HToBB_HToBB_TuneCP2_13TeV-madgraphMLM-pythia8__RunIIFall17NanoAODv7__PUFall17Fast_Nano02Apr2020_102X_mc2017_realistic_v8-v1.root",
                                                 pico_directory=pico_directory, n_entries=n_entries))
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory="/net/cms17/cms17r0/pico/NanoAODv7/nano/2018/mc", mc_nanoaod_filename="TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8__RunIIAutumn18NanoAODv7__Nano02Apr2020_102X_upgrade2018_realistic_v21-v1__100000__0969ED8B-18AE-4F4A-8C75-4637D3C688B7.root", 
                                                 signal_nanoaod_directory="/net/cms24/cms24r0/pico/NanoAODv7/nano/2018/SMS-TChiHH_2D_fastSimJmeCorrection", signal_nanoaod_filename="SMS-TChiHH_mChi-500_mLSP-0_HToBB_HToBB_TuneCP2_13TeV-madgraphMLM-pythia8__RunIIAutumn18NanoAODv7__PUFall18Fast_Nano02Apr2020_102X_upgrade2018_realistic_v21-v1.root",
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
