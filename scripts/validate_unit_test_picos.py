#!/bin/sh
''''exec python3 -u -- "$0" ${1+"$@"} # '''
import sys
import os
import glob
import subprocess
import threading
import time
import multiprocessing
import re
import argparse

def output_reader(process, commandOutput):
  for line in iter(process.stdout.readline, b''):
    output = line.decode('utf-8')
    #print(output, end='')
    commandOutput[0] += output

def runCommand(command):
  print("\n[Info] Running command: "+command)
  process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)

  commandOutput = ['']
  thread = threading.Thread(target=output_reader, args=(process,commandOutput))
  thread.start()
  thread.join()
  # Try to get poll
  for iTime in range(10):
    if (process.poll() == None): time.sleep(1)
    else: break
  #print(commandOutput[0])
  return process.poll(), commandOutput[0]

def getProcessTime(log_path):
  # process_information = [(command, return_code, execution_time (sec))]
  process_information = []
  with open(log_path) as log_file:
    count = 1
    for line in log_file:
      #print(re.search())
      match_object = re.search("\[Info\] command: (.*) return code: (.*) execution time: (.*) seconds", line)
      if match_object:
        command = match_object.group(1)
        return_code = match_object.group(2)
        execution_time = match_object.group(3)
        process_information.append([command, return_code, execution_time])
  return process_information

# ./run/process_nano.exe -f DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16NanoAODv9__106X_mcRun2_asymptotic_v17-v1__30000__0082C29D-E74C-024A-BE9B-97B29EE7A4A2.root -i /net/cms17/cms17r0/pico/NanoAODv9/nano/2016/mc -o unit_test_htozgamma_nanoaodv9/2016/mc --nent 30000
def strip_input(command):
  no_input_list = []
  input_arg = 0
  for arg in command.split():
    if input_arg == 2: input_arg = 0
    if input_arg == 1: input_arg += 1
    if arg == "-i": input_arg = 1
    if input_arg == 0:
      no_input_list.append(arg)
  return ''.join(no_input_list)

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='''Compares unit test pico files and logs. Prints speed of production and root file differences.''', formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('-o','--output_log_filename', required=True, help='Ouptut log filename.')
  parser.add_argument('-u','--unit_test_log_filename', required=True, help='Unit test production log filename.')
  parser.add_argument('-g','--golden_base_folder', required=True, help='Golden base folder for pico files.')
  parser.add_argument('-v','--validate_base_folder', required=True, help='Validate base folder for pico files.')
  args = parser.parse_args()

  validation_log_filename = args.output_log_filename
  log_filename = args.unit_test_log_filename
  golden_base_folder = args.golden_base_folder
  validate_base_folder = args.validate_base_folder
  #validation_log_filename = "validate_unit_test_higgsino_production.log"
  #log_filename = "unit_test_higgsino_production.log"
  #golden_base_folder = "/homes/jbkim/analysis/nano2pico/unit_test_higgsino"
  #validate_base_folder = "/homes/jbkim/analysis/nano2pico.variables/unit_test_higgsino"
  compare_script = os.path.dirname(os.path.realpath(__file__))+"/../root_scripts/compare_root_tuples.cxx"

  t0 = time.time()
  validation_log_file = open(validation_log_filename, 'w')

  # Get process times
  # process_information = [(command, return_code, execution_time (sec))]
  golden_log_path = golden_base_folder+"/"+log_filename
  golden_process_information = getProcessTime(golden_log_path)
  # process_information = [(command, return_code, execution_time (sec))]
  validate_log_path = validate_base_folder+"/"+log_filename
  validate_process_information = getProcessTime(validate_log_path)
  # Print process times
  for process in golden_process_information:
    if 'mkdir' in process[0]: continue
    command = process[0]
    golden_execution_time = process[2]
    validate_execution_time = -1
    for validate_process in validate_process_information:
      # Ignore input directory
      if strip_input(validate_process[0]) == strip_input(command): 
        validate_execution_time = validate_process[2]
        break
    print("Command: "+process[0])
    print("  Execution time: golden "+golden_execution_time+" sec vs validate "+validate_execution_time+" sec. ")
    print("  Difference: "+"{:.1f}".format((float(validate_execution_time)-float(golden_execution_time))*100./float(golden_execution_time))+" %")
    validation_log_file.write("Command: "+process[0]+'\n')
    validation_log_file.write("  Execution time: golden "+golden_execution_time+" sec vs validate "+validate_execution_time+" sec. \n")
    validation_log_file.write("  Difference: "+"{:.1f}".format((float(validate_execution_time)-float(golden_execution_time))*100./float(golden_execution_time))+" %\n")

  # Find folders to compare in golden folder
  compare_folders = []
  for year in os.listdir(golden_base_folder):
    if os.path.isfile(os.path.join(golden_base_folder,year)): continue
    for data_type in os.listdir(os.path.join(golden_base_folder,year)):
      if os.path.isfile(os.path.join(golden_base_folder,year,data_type)): continue
      for production in os.listdir(os.path.join(golden_base_folder,year,data_type)): 
        if production == "wgt_sums": continue
        if production == "corrections": continue
        if production == "raw_pico" and data_type != "data": continue
        path = os.path.join(year,data_type,production)
        compare_folders.append(path)
  print("[Info] Will compare following folders: "+", ".join(compare_folders))

  # Make command list
  command_list = []
  for compare_folder in compare_folders:
    golden_folder = golden_base_folder+"/"+compare_folder
    validate_folder = validate_base_folder+"/"+compare_folder
    for golden_filepath in glob.glob(golden_folder+"/*.root"):
      validate_filepath = validate_folder+"/"+os.path.basename(golden_filepath)
      command = 'root -b -l -q "'+compare_script+'+(\\"'+golden_filepath+'\\",\\"'+validate_filepath+'\\",\\"tree\\")"'
      command_list.append(command)

  # Print commands
  print("====Commands to run====")
  for command in command_list: print(command)
  print("====Commands to run====")
  print(str(len(command_list))+" commands to run")

  # Make sure so file for compare_script exists
  print("Making .so file for compare script")
  os.system("root -b -l -q "+compare_script+"++")

  ## Run command list
  #for command in command_list:
  #  runCommand(command)
  #  break

  # Run commands in multiple processes
  #pool = multiprocessing.Pool(processes=3)
  pool = multiprocessing.Pool()
  command_results = pool.map(runCommand, command_list)

  # Make script organize results
  print("====Command results====")
  for command_result in command_results:
    print(command_result[1])
    validation_log_file.write(command_result[1]+'\n')
  print("====Command results====")

  validation_log_file.close()

  print('\nProgram took %.0fm %.0fs.' % ((time.time()-t0)/60,(time.time()-t0)%60))
