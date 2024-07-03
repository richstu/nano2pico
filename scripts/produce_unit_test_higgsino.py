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

#def makeTestCommands(mc_nanoaod_directory, mc_nanoaod_filename, data_nanoaod_directory, data_nanoaod_filename, signal_nanoaod_directory, signal_nanoaod_filename, pico_directory, n_entries):
def makeTestCommands(data_nanoaod_directory, data_nanoaod_filename, pico_directory):

  higgsino_test_commands = [
    # mc commands
    #"mkdir -p "+pico_directory+"/mc/raw_pico",
    #"mkdir "+pico_directory+"/mc/wgt_sums",
    #"./run/process_nano.exe -f "+mc_nanoaod_filename+" -i "+mc_nanoaod_directory+" -o "+pico_directory+"/mc --nent "+str(n_entries),
    #"mkdir "+pico_directory+"/mc/corrections",
    #"./run/merge_corrections.exe "+pico_directory+"/mc/corrections/"+mc_nanoaod_filename+" "+pico_directory+"/mc/wgt_sums/wgt_sums_"+mc_nanoaod_filename,
    #"mkdir "+pico_directory+"/mc/unskimmed",
    #"./run/apply_corrections.exe -f raw_pico_"+mc_nanoaod_filename+" -i "+pico_directory+"/mc/raw_pico/ -c "+mc_nanoaod_filename,
    #"mkdir "+pico_directory+"/mc/met150",
    #"./scripts/skim_file.py -k met150 -i "+pico_directory+"/mc/unskimmed/pico_"+mc_nanoaod_filename+" -o "+pico_directory+"/mc/met150/",
    #"mkdir "+pico_directory+"/mc/merged_met150",
    #"./scripts/slim_and_merge.py -s txt/slim_rules/higmc.txt -i "+pico_directory+"/mc/met150/pico_met150_"+mc_nanoaod_filename+" -o "+pico_directory+"/mc/merged_met150/merged_"+mc_nanoaod_filename,

    # data commands
    "mkdir -p "+pico_directory+"/data/raw_pico",
    "mkdir "+pico_directory+"/data/wgt_sums",
    "./run/process_nano.exe -f "+data_nanoaod_filename+" -i "+data_nanoaod_directory+" -o "+pico_directory+"/data",
    "mkdir "+pico_directory+"/data/met150",
    "./scripts/skim_file.py -k met150 -i "+pico_directory+"/data/raw_pico/raw_pico_"+data_nanoaod_filename+" -o "+pico_directory+"/data/met150/",
    "mkdir "+pico_directory+"/data/merged_met150",
    "./scripts/slim_and_merge.py -s txt/slim_rules/higdata.txt -i "+pico_directory+"/data/met150/raw_pico_met150_"+data_nanoaod_filename+" -o "+pico_directory+"/data/merged_met150/merged_"+data_nanoaod_filename,

    # signal commands
    #"mkdir -p "+pico_directory+"/signal/raw_pico",
    #"mkdir "+pico_directory+"/signal/wgt_sums",
    #"./run/process_nano.exe -f "+signal_nanoaod_filename+" -i "+signal_nanoaod_directory+" -o "+pico_directory+"/signal",
    #"mkdir "+pico_directory+"/signal/corrections",
    #"./run/merge_corrections.exe "+pico_directory+"/signal/corrections/"+signal_nanoaod_filename+" "+pico_directory+"/signal/wgt_sums/wgt_sums_"+signal_nanoaod_filename,
    #"mkdir "+pico_directory+"/signal/unskimmed",
    #"./run/apply_corrections.exe -f raw_pico_"+signal_nanoaod_filename+" -i "+pico_directory+"/signal/raw_pico/ -c "+signal_nanoaod_filename,
    #"mkdir "+pico_directory+"/signal/met150",
    #"./scripts/skim_file.py -k met150 -i "+pico_directory+"/signal/unskimmed/pico_"+signal_nanoaod_filename+" -o "+pico_directory+"/signal/met150/",
    #"mkdir "+pico_directory+"/signal/merged_met150",
    #"./scripts/slim_and_merge.py -s txt/slim_rules/higmc.txt -i "+pico_directory+"/signal/met150/pico_met150_"+signal_nanoaod_filename+" -o "+pico_directory+"/signal/merged_met150/merged_"+signal_nanoaod_filename,
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

  n_entries = 1000 # about 3 min for process_nano.exe (~550 Hz). Note for signal, nent is -1.
  higgsino_test_commands = []
  higgsino_test_commands.extend(makeTestCommands(#mc_nanoaod_directory="/net/cms11/cms11r0/pico/NanoAODv9/nano/2016/mc", mc_nanoaod_filename="QCD_Pt-40ToInf_DoubleEMEnriched_MGG-80ToInf_TuneCP5_13TeV-pythia8__RunIISummer20UL16NanoAODv9__106X_mcRun2_asymptotic_v17-v1__70000__F93E4917-0EB7-F344-A545-0DFF3A907AE4.root", 
                                                 data_nanoaod_directory="/net/cms18/cms18r0/pico/NanoAODv9/nano/2016/mc", data_nanoaod_filename="SMS-TChiHH_mChi-1000_mLSP-0_HToGG_2D_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1__privateProduction.root", 
                                                 #signal_nanoaod_directory="/net/cms18/cms18r0/pico/NanoAODv9/nano/2016/mc", signal_nanoaod_filename="SMS-TChiHH_mChi-1000_mLSP-0_HToGG_2D_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1__privateProduction.root",
                                                 pico_directory=pico_directory))
  '''higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory="/net/cms17/cms17r0/pico/NanoAODv7/nano/2017/mc", mc_nanoaod_filename="TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8__RunIIFall17NanoAODv7__PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1__100000__0ED19AF0-B248-8344-91D7-B241CE0729FA.root", 
                                                 data_nanoaod_directory="/net/cms17/cms17r0/pico/NanoAODv7/nano/2017/data", data_nanoaod_filename="MET__Run2017B__02Apr2020-v1__230000__00DCCA4E-F5F1-F84D-A6EC-2956ACAB6E02.root", 
                                                 signal_nanoaod_directory="/net/cms24/cms24r0/pico/NanoAODv7/nano/2017/SMS-TChiHH_2D_fastSimJmeCorrection", signal_nanoaod_filename="SMS-TChiHH_mChi-500_mLSP-0_HToBB_HToBB_TuneCP2_13TeV-madgraphMLM-pythia8__RunIIFall17NanoAODv7__PUFall17Fast_Nano02Apr2020_102X_mc2017_realistic_v8-v1.root",
                                                 pico_directory=pico_directory, n_entries=n_entries))
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory="/net/cms17/cms17r0/pico/NanoAODv7/nano/2018/mc", mc_nanoaod_filename="TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8__RunIIAutumn18NanoAODv7__Nano02Apr2020_102X_upgrade2018_realistic_v21-v1__100000__0969ED8B-18AE-4F4A-8C75-4637D3C688B7.root", 
                                                 data_nanoaod_directory="/net/cms17/cms17r0/pico/NanoAODv7/nano/2018/data", data_nanoaod_filename="MET__Run2018A__02Apr2020-v1__20000__1F59D1E0-0193-E34B-91BE-4DF00AEA4FD6.root", 
                                                 signal_nanoaod_directory="/net/cms24/cms24r0/pico/NanoAODv7/nano/2018/SMS-TChiHH_2D_fastSimJmeCorrection", signal_nanoaod_filename="SMS-TChiHH_mChi-500_mLSP-0_HToBB_HToBB_TuneCP2_13TeV-madgraphMLM-pythia8__RunIIAutumn18NanoAODv7__PUFall18Fast_Nano02Apr2020_102X_upgrade2018_realistic_v21-v1.root",
                                                 pico_directory=pico_directory, n_entries=n_entries))
  '''
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
