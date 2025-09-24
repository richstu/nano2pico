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
      "./run/apply_corrections.exe -f raw_pico_"+mc_nanoaod_filename+" -i "+pico_directory+"/"+year+"/mc/raw_pico/ -c "+mc_nanoaod_filename
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

  n_entries = 6000 # about 1 min for process_nano.exe (~550 Hz). 
  higgsino_test_commands = []
  # For DY Samples
  test_file="DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v1__280000__52E5237A-EB6F-4F43-A91F-20E2EAAD9E7D.root"
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016APV/mc", mc_nanoaod_filename=test_file,
                                                 year = "2016APV", pico_directory=pico_directory, n_entries=n_entries))

  test_file="ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v1__30000__FA4ABBF4-D49C-604A-9ACF-EA5C6E7F9433.root"
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016APV/mc", mc_nanoaod_filename=test_file,
                                                 year = "2016APV", pico_directory=pico_directory, n_entries=n_entries))

  #ttbar samples
  test_file="TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v1__100000__7A786273-E40C-C24F-82D0-21133C646277.root"
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016APV/mc", mc_nanoaod_filename=test_file,
                                                 year = "2016APV", pico_directory=pico_directory, n_entries=n_entries))

  test_file="TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v2__2530000__810C1184-2FB4-FA49-96AE-46ADA95C28D3.root"
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016APV/mc", mc_nanoaod_filename=test_file,
                                                 year = "2016APV", pico_directory=pico_directory, n_entries=n_entries))

  #WW samples
  test_file="WW_TuneCP5_13TeV-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v1__120000__20B5B05D-3629-9446-838A-5A8C9F10F26F.root"
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016APV/mc", mc_nanoaod_filename=test_file,
                                                 year = "2016APV", pico_directory=pico_directory, n_entries=n_entries))

  test_file="WWG_TuneCP5_13TeV-amcatnlo-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v1__70000__388E7C19-DB4F-E54B-9C6C-62E90CB6BC34.root"
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016APV/mc", mc_nanoaod_filename=test_file,
                                                 year = "2016APV", pico_directory=pico_directory, n_entries=n_entries))

  #WZ samples
  test_file="WZ_TuneCP5_13TeV-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v1__280000__591D0E0B-90FC-984B-BCA0-F12266EFD6B4.root"
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016APV/mc", mc_nanoaod_filename=test_file,
                                                 year = "2016APV", pico_directory=pico_directory, n_entries=n_entries))

  test_file="WZG_TuneCP5_13TeV-amcatnlo-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v1__130000__55B79E91-F0E4-7045-98B0-7466223BC645.root"
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016APV/mc", mc_nanoaod_filename=test_file,
                                                 year = "2016APV", pico_directory=pico_directory, n_entries=n_entries))

  #ZZ samples
  test_file="ZZ_TuneCP5_13TeV-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v1__70000__3294B799-AA9D-3247-9FEC-851FA0F312E8.root"
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016APV/mc", mc_nanoaod_filename=test_file,
                                                 year = "2016APV", pico_directory=pico_directory, n_entries=n_entries))

  test_file="ZZGTo4L_TuneCP5_4f_NLO_13TeV-amcatnlo-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v2__2820000__A5D7C6DC-5AEC-3448-8AAF-C7AE5DB34B0C.root"
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016APV/mc", mc_nanoaod_filename=test_file,
                                                 year = "2016APV", pico_directory=pico_directory, n_entries=n_entries))


  #EWKZ samples
  test_file="ZGamma2JToGamma2L2J_EWK_MLL-50_MJJ-120_TuneCP5_13TeV-madgraph-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v2__2530000__B03ED76F-33EB-2543-95F5-55DE04815446.root"
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016APV/mc", mc_nanoaod_filename=test_file,
                                                 year = "2016APV", pico_directory=pico_directory, n_entries=n_entries))


  test_file="EWKZ2Jets_ZToLL_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v1__270000__7693DE5C-FCAE-5D4C-930C-45D4C17A7BCC.root"
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv9/nano/2016APV/mc", mc_nanoaod_filename=test_file,
                                                 year = "2016APV", pico_directory=pico_directory, n_entries=n_entries))

  ###########################################################################################################################################################################################################
  #Run 3 Overlap Removal
  ###########################################################################################################################################################################################################

  test_file="VBFto2L_MLL-50_TuneCP5_13p6TeV_madgraph-pythia8__Run3Summer22NanoAODv12__130X_mcRun3_2022_realistic_v5-v2__40000__288b9e0f-da35-45b3-b122-715e94d3c664.root"
  higgsino_test_commands.extend(makeTestCommands(mc_nanoaod_directory=input_pico_folder+"/NanoAODv12/nano/2022/mc",  mc_nanoaod_filename=test_file,
                                                 year = "2022", pico_directory=pico_directory, n_entries=n_entries))



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
