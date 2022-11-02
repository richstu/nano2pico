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


def check_target_production(string):
  if string != "higgsino" and string != "htozgamma" and string != "htogammagamma":
    raise argparse.ArgumentTypeError("Value has to be higgsino or htozgamma or htozgamma")
  return string

def check_file_exists(string):
  if not os.path.isfile(string):
    raise argparse.ArgumentTypeError("Input file does not exist")
  return string

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='''Runs nano2pico for testing''', formatter_class=argparse.RawTextHelpFormatter)
  
  parser.add_argument('-i','--input_file', help='Input NanoAOD file', required=True, type=check_file_exists)
  parser.add_argument('-o','--output_base_folder', help='Output base folder', default='out')
  parser.add_argument('-t','--target_production', help='Target production: higgsino, htozgamma, htogammagamma', required=True, type=check_target_production)
  parser.add_argument('-n','--number_events', help='Number of events', default='10000')
  parser.add_argument('-p', '--run_only_process', action="store_true", help='Run only process_nano')
  parser.add_argument('-m', '--run_only_merge', action="store_true", help='Run only merge_corrections')
  parser.add_argument('-a', '--run_only_apply', action="store_true", help='Run only apply_corrections')
  
  args = parser.parse_args()

  output_folder = os.path.join(args.output_base_folder,args.target_production)
  raw_pico_folder = os.path.join(output_folder,'raw_pico')
  wgt_sums_folder = os.path.join(output_folder,'wgt_sums')
  corrections_folder = os.path.join(output_folder,'corrections')
  unskimmed_folder = os.path.join(output_folder,'unskimmed')
  folder_list = [raw_pico_folder, wgt_sums_folder, corrections_folder, unskimmed_folder]
  for folder in folder_list: 
    if not os.path.isdir(folder): 
      print('Making folders '+folder)
      os.makedirs(folder)

  input_filename = os.path.basename(args.input_file)
  input_folder = os.path.dirname(args.input_file)

  run_all = (args.run_only_process == False) and (args.run_only_merge == False) and (args.run_only_apply == False)
  command_list = []
  if args.run_only_process or run_all:
    command_list.append('scons && run/process_nano.exe --in_file '+input_filename+' --in_dir '+input_folder+' --out_dir '+output_folder+' --nent '+args.number_events)
  if args.run_only_merge or run_all:
    command_list.append('scons && run/merge_corrections.exe '+corrections_folder+'/corr_'+input_filename+' '+wgt_sums_folder+'/wgt_sums_'+input_filename)
  if args.run_only_apply or run_all:
    command_list.append('scons && run/apply_corrections.exe --in_file raw_pico_'+input_filename+' --in_dir '+raw_pico_folder+'/ --corr_file corr_'+input_filename)

  for command in command_list:
    t0 = time.time()
    return_code =''
    output = ''
    return_code, output = runCommand(command)
    execution_time = time.time()-t0
    print("[Info] command: "+command)
    print("[Info] return code: "+str(return_code))
    print("[Info] execution time: "+"{:.1f}".format(execution_time)+" seconds\n")
    #print("[Info] output:")
    #print(output)
