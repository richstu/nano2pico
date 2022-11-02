#!/usr/bin/env python
from __future__ import print_function
import os
import glob
import sys
import argparse
import subprocess
import datetime
import threading
import time

def runCommand(command):
  print("\nRunning command: "+command)
  process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  out, err = process.communicate()
  if err != "": print(err.rstrip())
  return out.rstrip(),err.rstrip()

def output_reader(process, log_file):
  for line in iter(process.stdout.readline, b''):
    output = line.decode('utf-8')
    print(output, end='')
    log_file.write(output)
def runAndLogCommand(command, log_file):
  print("\nRunning command: "+command)
  log_file.write("\nRunning command: "+command+'\n')
  process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
  thread = threading.Thread(target=output_reader, args=(process,log_file))
  thread.start()
  thread.join()
  # Try to get poll
  for iTime in range(60):
    if (process.poll() == None): time.sleep(1)
    else: break
  return int(process.poll())

def processSteps(process_commands, YEAR, PRODUCTION_NAME, STEP_FILEBASENAME, LOG_FILENAME, PICO_DIR, NANOAOD_VERSION, FIRST_COMMAND, NO_RUN, tag):
  # Make command list for book keeping
  commandList = []
  for commands in process_commands:
    commandList.append("&&".join(commands))

  # Read step_file
  # Make file command list that shows done commands
  diskCommandList = []
  if os.path.isfile(STEP_FILEBASENAME):
    with open(STEP_FILEBASENAME) as step_file:
      for line in step_file:
        #if line[0] == "#": continue
        #if len(line.strip()) == "": continue
        diskCommandList.append(line.strip())

  # Make list command status
  didCommandRun = []
  for iStep in xrange(len(process_commands)):
    #print(commandList[iStep], diskCommandList)
    if commandList[iStep] in diskCommandList: didCommandRun.append(True)
    else: didCommandRun.append(False)
  # Comment out commands not in commandList
  for iCommand in xrange(len(diskCommandList)):
    if diskCommandList[iCommand].strip() == "": continue
    if diskCommandList[iCommand][0] == "#": continue
    if diskCommandList[iCommand] not in commandList:
      diskCommandList[iCommand] = '#not_in_command_list '+diskCommandList[iCommand]
  # Make a link between diskCommandList and iStep
  # stepToDiskCommand[iStep] = iLine
  stepToDiskCommand= {}
  for iCommand in xrange(len(diskCommandList)):
    # Search for command in commandList
    for iStep in xrange(len(commandList)):
      if commandList[iStep] == diskCommandList[iCommand]: stepToDiskCommand[iStep] = iCommand

  # Return if all commands have run
  if all(didCommandRun):
    print('All commands have run for '+tag)
    return FIRST_COMMAND

  # Print commands that was run and to run
  if FIRST_COMMAND:
    for iStep in xrange(len(process_commands)):
      if didCommandRun[iStep]: continue
      if iStep != 0:
        print('[Step '+str(iStep-1)+'] Previous commands')
        print('\n'.join(process_commands[iStep-1]))
        print('\n')
      print('[Step '+str(iStep)+'] Next command to run')
      print('\n'.join(process_commands[iStep]))
      raw_input('Press Enter to continue')
      FIRST_COMMAND = False
      break

  # Backup step file
  fileIndex = 1
  if os.path.isfile(STEP_FILEBASENAME):
    while (1):
      filename = STEP_FILEBASENAME+"."+str(fileIndex)
      if not os.path.isfile(filename):
        os.system('mv '+STEP_FILEBASENAME+' '+STEP_FILEBASENAME+"."+str(fileIndex))
        break
      fileIndex += 1

  # If log file exists, move log file to backup
  fileIndex = 1
  if os.path.isfile(LOG_FILENAME):
    while (1):
      filename = LOG_FILENAME+"."+str(fileIndex)
      if not os.path.isfile(filename):
        os.system('mv '+LOG_FILENAME+' '+LOG_FILENAME+"."+str(fileIndex))
        break
      fileIndex += 1
  # Open log file
  if not NO_RUN: log_file = open(LOG_FILENAME, 'w')

  # Do steps
  for iStep in xrange(len(process_commands)):
    if didCommandRun[iStep]: continue
    # Run commands
    error = 0
    for command in process_commands[iStep]:
      print('Will run below command')
      print('  '+command)
      if NO_RUN: continue
      #error += os.system(command)
      error += runAndLogCommand(command, log_file)
    didCommandRun[iStep] = True
    # Log
    if NO_RUN: continue
    # Update diskCommandList
    isAfterStep = False
    insertIndex = -1
    for iLine in xrange(len(diskCommandList)):
      # Adds a step in log file if there is an after step
      for step in stepToDiskCommand:
        if step > iStep:
         isAfterStep = True
         insertIndex = stepToDiskCommand[step]
         break
    if isAfterStep:
      diskCommandList.insert(max(0,insertIndex), commandList[iStep])
    else:
      diskCommandList.append(commandList[iStep])
    # Update stepToDiskCommand
    for iCommand in xrange(len(diskCommandList)):
      # Search for command in commandList
      for t_iStep in xrange(len(commandList)):
        if commandList[t_iStep] == diskCommandList[iCommand]: stepToDiskCommand[t_iStep] = iCommand
    # Update step file
    with open(STEP_FILEBASENAME,'w') as step_file:
      for iLine in xrange(len(diskCommandList)):
        step_file.write(diskCommandList[iLine]+'\n')

    # Check if there was an error
    if error != 0:
      print('There was an error')
      #os.system(notify_script+' "There was an error on step '+str(iStep)+' for '+tag+'"')
      runAndLogCommand(notify_script+' "There was an error on step '+str(iStep)+' for '+tag+'"',log_file)
      sys.exit(1)

  if not NO_RUN: log_file.close()

  return FIRST_COMMAND

def processMc(YEAR, PRODUCTION_NAME, STEP_FILEBASENAME, LOG_FILENAME, PICO_DIR, NANOAOD_VERSION, FIRST_COMMAND, notify_script):
  print("[Info] pico dir: "+PICO_DIR)
  print("[Info] step file: "+STEP_FILEBASENAME)
  print("[Info] log file: "+LOG_FILENAME)
  YEAR = str(YEAR)
  mc_tag=PRODUCTION_NAME+'_'+YEAR+'_mc'
  # Add mc commands
  process_commands = [
    #0
    [notify_script+' "Start process nano '+mc_tag+'"',
    './scripts/write_process_nano_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/nano/'+YEAR+'/mc/ --production '+PRODUCTION_NAME+' --dataset_list txt/datasets/'+NANOAOD_VERSION+'_htozgamma_'+YEAR+'_mc_dataset_paths --tag '+mc_tag,
    'auto_submit_jobs.py process_nano_cmds_'+mc_tag+'.json -c scripts/check_process_nano_job.py -f',
    notify_script+' "Finished process nano '+mc_tag+'"'], 
    
    #1
    [notify_script+' "Start merge corrections '+mc_tag+'"',
    './scripts/merge_corrections.py --wgt_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/wgt_sums/ --corr_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/corrections/',
    notify_script+' "Finished merge corrections '+mc_tag+'"'],
    
    #2
    [notify_script+' "Start applied corrections '+mc_tag+'"',
    './scripts/write_apply_corrections_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/raw_pico/ --tag '+mc_tag,
    'auto_submit_jobs.py '+mc_tag+'_apply_corr_cmds.json -c scripts/check_apply_corrections_job.py -f',
    notify_script+' "Finished applied corrections '+mc_tag+'"'],

    #3
    [notify_script+' "Start skim llg '+mc_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/unskimmed/ --skim_name llg --tag '+mc_tag,
    'auto_submit_jobs.py skim_llg_cmds_'+mc_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim llg '+mc_tag+'"'],
    
    #4
    [notify_script+' "Start merge llg '+mc_tag+'"',
    './scripts/write_slim_and_merge_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/skim_llg/ --slim_name zgmc --tag '+mc_tag,
    'auto_submit_jobs.py '+mc_tag+'_slim_zgmc_llg_cmds.json -f',
    './scripts/confirm_slim.py '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/skim_llg '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/merged_zgmc_llg',
    notify_script+' "Finished merge llg '+mc_tag+'"'],

  ]
  return processSteps(process_commands, YEAR, PRODUCTION_NAME, STEP_FILEBASENAME, LOG_FILENAME, PICO_DIR, NANOAOD_VERSION, FIRST_COMMAND, NO_RUN, mc_tag)

def processData(YEAR, PRODUCTION_NAME, STEP_FILEBASENAME, LOG_FILENAME, PICO_DIR, NANOAOD_VERSION, FIRST_COMMAND, notify_script):
  YEAR = str(YEAR)
  # Add signal commands
  data_tag=PRODUCTION_NAME+'_'+YEAR+'_data'
  print("[Info] pico dir: "+PICO_DIR)
  print("[Info] step file: "+STEP_FILEBASENAME)
  print("[Info] log file: "+LOG_FILENAME)
  process_commands = [
    # signal
    #0
    [notify_script+' "Started process nano '+data_tag+'"',
    './scripts/write_process_nano_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/nano/'+YEAR+'/data/ --production '+PRODUCTION_NAME+' --dataset_list txt/datasets/htozgamma_data_infile_list.txt --list_format filename --tag '+data_tag,
    'auto_submit_jobs.py process_nano_cmds_'+data_tag+'.json -c scripts/check_data_process_nano_job.py -f',
    notify_script+' "Finished process nano '+data_tag+'"',
    'ln -s '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/raw_pico '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/unskimmed'],
    
    #1
    [notify_script+' "Started skim llg '+data_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/raw_pico/ --skim_name llg --tag '+data_tag,
    'auto_submit_jobs.py skim_llg_cmds_'+data_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim llg '+data_tag+'"'],
    
    #2
    [notify_script+' "Started merge llg '+data_tag+'"',
    './scripts/write_slim_and_merge_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/skim_llg/ --slim_name zgdata --tag '+data_tag,
    'auto_submit_jobs.py '+data_tag+'_slim_zgdata_llg_cmds.json -c scripts/check_slim_and_merge.py -f',
    './scripts/confirm_slim.py '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/skim_llg '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/merged_zgdata_llg',
    notify_script+' "Finished merge llg '+data_tag+'"'],

  ]

  return processSteps(process_commands, YEAR, PRODUCTION_NAME, STEP_FILEBASENAME, LOG_FILENAME, PICO_DIR, NANOAOD_VERSION, FIRST_COMMAND, NO_RUN, data_tag)

if __name__ == '__main__':

  parser = argparse.ArgumentParser(description='''\
Script to produce zgamma piocs from NanoAOD. 
Requirements
- A git tag to do production
- dataset_list files that have NanoAOD names to process. 
  - Ex) txt/datasets/NanoAODv7_zgamma_mc_dataset_paths, txt/datasets/zgamma_data_infile_list.txt 
- slim_rule files that has rules for branches
  - Ex) txt/slim_rules/zgmc.txt, txt/slim_rules/zgdata.txt
- Script assumes it will be run in the directory that has .git.

Below output files are saved in BASE_FOLDERNAME/NANOAOD_VERSION/TAG_NAME/
- pico ntuples.
- step files which saves which step production is in. 
- log files which saves log of production.

Below are the steps the script follows.
1. Checkout a git tag to do production. 
2. Run commands for production, where used commands are saved in *.step files.

Note that the script stops when there is an error. 
A notification will be sent by email to $USER@hep.physics.ucsb.edu using scripts/sendEmail.sh.
To send the email elsewhere use the --email option.
After one fixes the error, rerun the script to resume to the next production command.

The folder structure is shown below
NanoAOD files: BASE_FOLDERNAME/NANOAOD_VERSION/nano/(2016,2017,2018)/(data,mc,signal)/*.root
  Ex) /net/cms17/cms17r0/pico/NanoAODv7/nano/2016/data/*.root
Pico files: BASE_FOLDERNAME/NANOAOD_VERSION/TAG_NAME/(2016,2017,2018)/(data,mc,signal)/SKIM_NAME/*.root
''', formatter_class=argparse.RawTextHelpFormatter)

  parser.add_argument('-t','--tag_name', required=True, help='Tag name for production in git')
  parser.add_argument('-n','--nanoaod_version', required=True, help='Nanoaod version for production')
  parser.add_argument('-b','--base_foldername', required=True, help='Base folder for ntuple files. Ex) /net/cms17/cms17r0/pico')
  parser.add_argument('-f', '--fake_run', action="store_true", help='Do not run commands. Only print commands to run.')
  parser.add_argument('--use_telegram', action="store_true", help='Uses telegram script to notify about steps. Requires telegram setup.')
  parser.add_argument('--email', help='Uses email to notify about steps. Type in your email.')
  
  args = parser.parse_args()

  # Check if base_folder exists
  if not os.path.exists(args.base_foldername):
    print('[Error] Could not find base folder: '+args.base_foldername+'. Exiting.')
    sys.exit()

  # Check if dataset_list (nanoaod files to process) exist
  dataset_list_files = [
    'txt/datasets/'+args.nanoaod_version+'_htozgamma_2016_mc_dataset_paths',
    'txt/datasets/'+args.nanoaod_version+'_htozgamma_2017_mc_dataset_paths',
    'txt/datasets/'+args.nanoaod_version+'_htozgamma_2018_mc_dataset_paths',
    'txt/datasets/htozgamma_data_infile_list.txt'
    ]
  for dataset_list_file in dataset_list_files:
    if not os.path.exists(dataset_list_file): 
      print('[Error] '+dataset_list_file+' does not exist. Existing.')
      sys.exit()

  # Check if slim rules exist
  slim_rule_files = [
    'txt/slim_rules/zgmc.txt',
    'txt/slim_rules/zgdata.txt',
    ]
  for slim_rule_file in slim_rule_files:
    if not os.path.exists(slim_rule_file): 
      print('[Error] '+slim_rule_file+' does not exist. Existing.')
      sys.exit()

  # Search if tag exists and is pushed.
  output, error = runCommand('git tag -l '+args.tag_name)
  if args.tag_name not in output:
    print("[Error] tag_name: "+args.tag_name+" does not exist. Exiting.")
    sys.exit()
  #output, error = runCommand('git ls-remote origin refs/tags/'+args.tag_name)
  #if args.tag_name not in output:
  #  print("[Error] tag_name: "+args.tag_name+" does not exist on github. Exiting.")
  #  sys.exit()
  
  # Checkout tag
  output, error = runCommand('git checkout '+args.tag_name)
  if 'Aborting' in error:
    print("[Error] Could not checkout "+args.tag_name+". Exiting.")
    sys.exit()

  # Do compile after checkout
  if os.system('scons') != 0:
    print("[Error] Could not compile. Existing")
    sys.exit()

  # Make output BASE_FOLDERNAME/NANOAOD_VERSION/TAG_NAME folder
  print("[Info] Trying to make "+args.base_foldername+'/'+args.nanoaod_version+'/'+args.tag_name)
  if not os.path.exists(args.base_foldername+'/'+args.nanoaod_version+'/'+args.tag_name):
    os.makedirs(args.base_foldername+'/'+args.nanoaod_version+'/'+args.tag_name)
    print("[Info] Made "+args.base_foldername+'/'+args.nanoaod_version+'/'+args.tag_name+' folder.')
  else:
    print("[Info] "+args.base_foldername+'/'+args.nanoaod_version+'/'+args.tag_name+" exists. No need for making folder.")

  # Set variables
  PRODUCTION_NAME = args.tag_name #Ex) 'htozgamma_klamath'
  PICO_DIR = args.base_foldername #Ex) '/net/cms17/cms17r0/pico'
  NANOAOD_VERSION = args.nanoaod_version # Ex) 'NanoAODv9'
  NO_RUN = args.fake_run
  # To prompt for first command
  FIRST_COMMAND = True

  # Set notification script
  if args.use_telegram: notify_script = 'sendTelegramMessage.py'
  else: 
    if args.email:
      notify_script = os.path.dirname(os.path.realpath(__file__))+'/sendEmail.sh '+args.email
    else:
      notify_script = os.path.dirname(os.path.realpath(__file__))+'/sendEmail.sh '+os.environ.get("USER")+'@hep.physics.ucsb.edu'

  # Main steps
  # - Process ntuples
  # - Make skims
  # - Make slim and skims
  # Note, there is a disk list of commands in step files. 
  # There is also a memory list of commands from this script. 
  # They are compared and missing steps are ran.
  # Step is the step that is running. Will run the next step.
  FIRST_COMMAND = processMc(YEAR=2016, PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME=PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/produce_zgamma_picos.py.'+PRODUCTION_NAME+'.2016.mc.step', LOG_FILENAME= PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/produce_zgamma_picos.py.'+PRODUCTION_NAME+'.2016.mc.log', PICO_DIR=PICO_DIR, NANOAOD_VERSION=NANOAOD_VERSION, FIRST_COMMAND=FIRST_COMMAND, notify_script=notify_script)
  FIRST_COMMAND = processMc(YEAR='2016APV', PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME=PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/produce_zgamma_picos.py.'+PRODUCTION_NAME+'.2016.mc.step', LOG_FILENAME= PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/produce_zgamma_picos.py.'+PRODUCTION_NAME+'.2016.mc.log', PICO_DIR=PICO_DIR, NANOAOD_VERSION=NANOAOD_VERSION, FIRST_COMMAND=FIRST_COMMAND, notify_script=notify_script)
  FIRST_COMMAND = processMc(YEAR=2017, PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME=PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/produce_zgamma_picos.py.'+PRODUCTION_NAME+'.2017.mc.step', LOG_FILENAME= PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/produce_zgamma_picos.py.'+PRODUCTION_NAME+'.2017.mc.log', PICO_DIR=PICO_DIR, NANOAOD_VERSION=NANOAOD_VERSION, FIRST_COMMAND=FIRST_COMMAND, notify_script=notify_script)
  FIRST_COMMAND = processMc(YEAR=2018, PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME=PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/produce_zgamma_picos.py.'+PRODUCTION_NAME+'.2018.mc.step', LOG_FILENAME= PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/produce_zgamma_picos.py.'+PRODUCTION_NAME+'.2018.mc.log', PICO_DIR=PICO_DIR, NANOAOD_VERSION=NANOAOD_VERSION, FIRST_COMMAND=FIRST_COMMAND, notify_script=notify_script)

  FIRST_COMMAND = processData(YEAR=2016, PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME=PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/produce_zgamma_picos.py.'+PRODUCTION_NAME+'.2016.data.step', LOG_FILENAME= PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/produce_zgamma_picos.py.'+PRODUCTION_NAME+'.2016.data.log', PICO_DIR=PICO_DIR, NANOAOD_VERSION=NANOAOD_VERSION, FIRST_COMMAND=FIRST_COMMAND, notify_script=notify_script)
  FIRST_COMMAND = processData(YEAR=2017, PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME=PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/produce_zgamma_picos.py.'+PRODUCTION_NAME+'.2017.data.step', LOG_FILENAME= PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/produce_zgamma_picos.py.'+PRODUCTION_NAME+'.2017.data.log', PICO_DIR=PICO_DIR, NANOAOD_VERSION=NANOAOD_VERSION, FIRST_COMMAND=FIRST_COMMAND, notify_script=notify_script)
  FIRST_COMMAND = processData(YEAR=2018, PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME=PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/produce_zgamma_picos.py.'+PRODUCTION_NAME+'.2018.data.step', LOG_FILENAME= PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/produce_zgamma_picos.py.'+PRODUCTION_NAME+'.2018.data.log', PICO_DIR=PICO_DIR, NANOAOD_VERSION=NANOAOD_VERSION, FIRST_COMMAND=FIRST_COMMAND, notify_script=notify_script)

  # Change permission of directories
  os.system("find "+PICO_DIR+"/"+NANOAOD_VERSION+'/'+PRODUCTION_NAME+" -type d -exec chmod 775 {} \;")
