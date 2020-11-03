#!/usr/bin/env python
import os
import glob
import sys

def processSteps(process_commands, YEAR, PRODUCTION_NAME, STEP_FILEBASENAME, PICO_DIR, NANOAOD_VERSION, FIRST_COMMAND, NO_RUN, tag):
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

  # Backup log file
  fileIndex = 1
  if os.path.isfile(STEP_FILEBASENAME):
    while (1):
      filename = STEP_FILEBASENAME+"."+str(fileIndex)
      if not os.path.isfile(filename):
        os.system('mv '+STEP_FILEBASENAME+' '+STEP_FILEBASENAME+"."+str(fileIndex))
        break
      fileIndex += 1

   # Do steps
  for iStep in xrange(len(process_commands)):
    if didCommandRun[iStep]: continue
    # Run commands
    error = 0
    for command in process_commands[iStep]:
      print('Will run below command')
      print('  '+command)
      if NO_RUN: continue
      error += os.system(command)
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
    # Update log file
    with open(STEP_FILEBASENAME,'w') as step_file:
      for iLine in xrange(len(diskCommandList)):
        step_file.write(diskCommandList[iLine]+'\n')

    # Check if there was an error
    if error != 0:
      print('There was an error')
      os.system(notify_script+' "There was an error on step '+str(iStep)+' for '+tag+'"')
      sys.exit(1)
  return FIRST_COMMAND

def processMc(YEAR, PRODUCTION_NAME, STEP_FILEBASENAME, PICO_DIR, NANOAOD_VERSION, FIRST_COMMAND, notify_script):
  YEAR = str(YEAR)
  mc_tag=PRODUCTION_NAME+'_'+YEAR+'_mc'
  # Add mc commands
  process_commands = [
    #0
    [notify_script+' "Start process nano '+mc_tag+'"',
    './scripts/write_process_nano_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/nano/'+YEAR+'/mc/ --production '+PRODUCTION_NAME+' --dataset_list txt/datasets/'+NANOAOD_VERSION+'_higgsino_'+YEAR+'_mc_dataset_paths --tag '+mc_tag,
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
    [notify_script+' "Start skim met150 '+mc_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/unskimmed/ --skim_name met150 --tag '+mc_tag,
    'auto_submit_jobs.py skim_met150_cmds_'+mc_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim met150 '+mc_tag+'"'],
    
    #4
    [notify_script+' "Start merge met150 '+mc_tag+'"',
    './scripts/write_slim_and_merge_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/skim_met150/ --slim_name higmc --tag '+mc_tag,
    'auto_submit_jobs.py '+mc_tag+'_slim_higmc_met150_cmds.json -f',
    './scripts/confirm_slim.py '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/skim_met150 '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/merged_higmc_met150',
    notify_script+' "Finished merge met150 '+mc_tag+'"'],

    #5
    [notify_script+' "Start skim preselect '+mc_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/merged_higmc_met150/ --skim_name preselect --tag '+mc_tag,
    'auto_submit_jobs.py skim_preselect_cmds_'+mc_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim preselect '+mc_tag+'"'],
    
    #6
    [notify_script+' "Start skim higloose '+mc_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/merged_higmc_met150/ --skim_name higloose --tag '+mc_tag,
    'auto_submit_jobs.py skim_higloose_cmds_'+mc_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim higloose '+mc_tag+'"'],

    #7 1l control sample skim
    [notify_script+' "Start skim higlep1T '+mc_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/unskimmed/ --skim_name higlep1T --tag '+mc_tag,
    'auto_submit_jobs.py skim_higlep1T_cmds_'+mc_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim higlep1T '+mc_tag+'"'],

    #8 1l control sample merge
    [notify_script+' "Start merge higlep1T '+mc_tag+'"',
    './scripts/write_slim_and_merge_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/skim_higlep1T/ --slim_name higmc --tag '+mc_tag,
    'auto_submit_jobs.py '+mc_tag+'_slim_higmc_higlep1T_cmds.json -f',
    './scripts/confirm_slim.py '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/skim_higlep1T '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/merged_higmc_higlep1T',
    notify_script+' "Finished merge higlep1T '+mc_tag+'"'],

    #9 2l control sample skim
    [notify_script+' "Start skim higlep2T '+mc_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/unskimmed/ --skim_name higlep2T --tag '+mc_tag,
    'auto_submit_jobs.py skim_higlep2T_cmds_'+mc_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim higlep2T '+mc_tag+'"'],

    #10 2l control sample merge
    [notify_script+' "Start merge higlep2T '+mc_tag+'"',
    './scripts/write_slim_and_merge_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/skim_higlep2T/ --slim_name higmc --tag '+mc_tag,
    'auto_submit_jobs.py '+mc_tag+'_slim_higmc_higlep2T_cmds.json -f',
    './scripts/confirm_slim.py '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/skim_higlep2T '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/merged_higmc_higlep2T',
    notify_script+' "Finished merge higlep2T '+mc_tag+'"'],

    #11 qcd control sample skim and merge from met150
    [notify_script+' "Start skim higqcd '+mc_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/mc/merged_higmc_met150/ --skim_name higqcd --tag '+mc_tag,
    'auto_submit_jobs.py skim_higqcd_cmds_'+mc_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim higqcd '+mc_tag+'"'],
  ]
  return processSteps(process_commands, YEAR, PRODUCTION_NAME, STEP_FILEBASENAME, PICO_DIR, NANOAOD_VERSION, FIRST_COMMAND, NO_RUN, mc_tag)

def processSignal(YEAR, PRODUCTION_NAME, STEP_FILEBASENAME, PICO_DIR, NANOAOD_VERSION, FIRST_COMMAND, notify_script):
  YEAR = str(YEAR)
  # Add signal commands
  sig_tag=PRODUCTION_NAME+'_'+YEAR+'_sig'
  process_commands = [
    # signal
    #0
    [notify_script+' "Started process nano '+sig_tag+'"',
    './scripts/write_process_nano_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/nano/'+YEAR+'/SMS-TChiHH_2D/ --production '+PRODUCTION_NAME+' --tag '+sig_tag,
    'auto_submit_jobs.py process_nano_cmds_'+sig_tag+'.json -c scripts/check_process_nano_job.py -f',
    notify_script+' "Finished process nano '+sig_tag+'"'],
    
    #1
    [notify_script+' "Started merge corrections '+sig_tag+'"',
    './scripts/merge_corrections.py --wgt_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/wgt_sums/ --corr_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/corrections/',
    notify_script+' "Finished merge corrections '+sig_tag+'"'],
    
    #2
    [notify_script+' "Started applied corrections '+sig_tag+'"',
    './scripts/write_apply_corrections_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/raw_pico/ --tag '+sig_tag,
    'auto_submit_jobs.py '+sig_tag+'_apply_corr_cmds.json -c scripts/check_apply_corrections_job.py -f',
    notify_script+' "Finished applied corrections '+sig_tag+'"'],
    
    #3
    [notify_script+' "Started skim met150 '+sig_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/unskimmed/ --skim_name met150 --tag '+sig_tag,
    'auto_submit_jobs.py skim_met150_cmds_'+sig_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim met150 '+sig_tag+'"'],
    
    #4
    [notify_script+' "Started merge met150 '+sig_tag+'"',
    './scripts/write_slim_and_merge_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/skim_met150/ --slim_name higmc --tag '+sig_tag,
    'auto_submit_jobs.py '+sig_tag+'_slim_higmc_met150_cmds.json -f',
    './scripts/confirm_slim.py '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/skim_met150 '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/merged_higmc_met150',
    notify_script+' "Finished merge met150 '+sig_tag+'"'],
    
    #5
    [notify_script+' "Started skim preselect '+sig_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/merged_higmc_met150 --skim_name preselect --tag '+sig_tag,
    'auto_submit_jobs.py skim_preselect_cmds_'+sig_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim preselect '+sig_tag+'"'],

    #6
    [notify_script+' "Started skim higloose '+sig_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/merged_higmc_met150 --skim_name higloose --tag '+sig_tag,
    'auto_submit_jobs.py skim_higloose_cmds_'+sig_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim higloose '+sig_tag+'"'],

    #7 1l control sample skim
    [notify_script+' "Start skim higlep1T '+sig_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/unskimmed/ --skim_name higlep1T --tag '+sig_tag,
    'auto_submit_jobs.py skim_higlep1T_cmds_'+sig_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim higlep1T '+sig_tag+'"'],

    #8 1l control sample merge
    [notify_script+' "Start merge higlep1T '+sig_tag+'"',
    './scripts/write_slim_and_merge_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/skim_higlep1T/ --slim_name higmc --tag '+sig_tag,
    'auto_submit_jobs.py '+sig_tag+'_slim_higmc_higlep1T_cmds.json -f',
    './scripts/confirm_slim.py '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/skim_higlep1T '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/merged_higmc_higlep1T',
    notify_script+' "Finished merge higlep1T '+sig_tag+'"'],

    #9 2l control sample skim
    [notify_script+' "Start skim higlep2T '+sig_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/unskimmed/ --skim_name higlep2T --tag '+sig_tag,
    'auto_submit_jobs.py skim_higlep2T_cmds_'+sig_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim higlep2T '+sig_tag+'"'],

    #10 2l control sample merge
    [notify_script+' "Start merge higlep2T '+sig_tag+'"',
    './scripts/write_slim_and_merge_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/skim_higlep2T/ --slim_name higmc --tag '+sig_tag,
    'auto_submit_jobs.py '+sig_tag+'_slim_higmc_higlep2T_cmds.json -f',
    './scripts/confirm_slim.py '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/skim_higlep2T '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/merged_higmc_higlep2T',
    notify_script+' "Finished merge higlep2T '+sig_tag+'"'],

    #11 qcd control sample skim and merge from met150
    [notify_script+' "Start skim higqcd '+sig_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/merged_higmc_met150/ --skim_name higqcd --tag '+sig_tag,
    'auto_submit_jobs.py skim_higqcd_cmds_'+sig_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim higqcd '+sig_tag+'"'],

  ]
  return processSteps(process_commands, YEAR, PRODUCTION_NAME, STEP_FILEBASENAME, PICO_DIR, NANOAOD_VERSION, FIRST_COMMAND, NO_RUN, sig_tag)

def processData(YEAR, PRODUCTION_NAME, STEP_FILEBASENAME, PICO_DIR, NANOAOD_VERSION, FIRST_COMMAND, notify_script):
  YEAR = str(YEAR)
  # Add signal commands
  data_tag=PRODUCTION_NAME+'_'+YEAR+'_data'
  process_commands = [
    # signal
    #0
    [notify_script+' "Started process nano '+data_tag+'"',
    './scripts/write_process_nano_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/nano/'+YEAR+'/data/ --production '+PRODUCTION_NAME+' --dataset_list txt/datasets/higgsino_data_infile_list.txt --list_format filename --tag '+data_tag,
    'auto_submit_jobs.py process_nano_cmds_'+data_tag+'.json -c scripts/check_data_process_nano_job.py -f',
    notify_script+' "Finished process nano '+data_tag+'"'],
    
    #1
    [notify_script+' "Started skim met150 '+data_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/raw_pico/ --skim_name met150 --tag '+data_tag,
    'auto_submit_jobs.py skim_met150_cmds_'+data_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim met150 '+data_tag+'"'],
    
    #2
    [notify_script+' "Started merge met150 '+data_tag+'"',
    './scripts/write_slim_and_merge_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/skim_met150/ --slim_name higdata --tag '+data_tag,
    'auto_submit_jobs.py '+data_tag+'_slim_higdata_met150_cmds.json -f',
    './scripts/confirm_slim.py '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/skim_met150 '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/merged_higdata_met150',
    notify_script+' "Finished merge met150 '+data_tag+'"'],
    
    #3
    [notify_script+' "Started skim preselect '+data_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/merged_higdata_met150 --skim_name preselect --tag '+data_tag,
    'auto_submit_jobs.py skim_preselect_cmds_'+data_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim preselect '+data_tag+'"'],

    #4
    [notify_script+' "Started skim higloose '+data_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/merged_higdata_met150 --skim_name higloose --tag '+data_tag,
    'auto_submit_jobs.py skim_higloose_cmds_'+data_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim higloose '+data_tag+'"'],

    #5 1l control sample skim
    [notify_script+' "Start skim higlep1T '+data_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/raw_pico/ --skim_name higlep1T --tag '+data_tag,
    'auto_submit_jobs.py skim_higlep1T_cmds_'+data_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim higlep1T '+data_tag+'"'],

    #6 1l control sample merge
    [notify_script+' "Start merge higlep1T '+data_tag+'"',
    './scripts/write_slim_and_merge_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/skim_higlep1T/ --slim_name higdata --tag '+data_tag,
    'auto_submit_jobs.py '+data_tag+'_slim_higdata_higlep1T_cmds.json -f',
    './scripts/confirm_slim.py '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/skim_higlep1T '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/merged_higdata_higlep1T',
    notify_script+' "Finished merge higlep1T '+data_tag+'"'],

    #7 2l control sample skim
    [notify_script+' "Start skim higlep2T '+data_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/raw_pico/ --skim_name higlep2T --tag '+data_tag,
    'auto_submit_jobs.py skim_higlep2T_cmds_'+data_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim higlep2T '+data_tag+'"'],

    #8 2l control sample merge
    [notify_script+' "Start merge higlep2T '+data_tag+'"',
    './scripts/write_slim_and_merge_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/skim_higlep2T/ --slim_name higdata --tag '+data_tag,
    'auto_submit_jobs.py '+data_tag+'_slim_higdata_higlep2T_cmds.json -f',
    './scripts/confirm_slim.py '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/skim_higlep2T '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/merged_higdata_higlep2T',
    notify_script+' "Finished merge higlep2T '+data_tag+'"'],

    #9 qcd control sample skim and merge from met150
    [notify_script+' "Start skim higqcd '+data_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/'+NANOAOD_VERSION+'/'+PRODUCTION_NAME+'/'+YEAR+'/data/merged_higdata_met150/ --skim_name higqcd --tag '+data_tag,
    'auto_submit_jobs.py skim_higqcd_cmds_'+data_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim higqcd '+data_tag+'"'],

  ]

  return processSteps(process_commands, YEAR, PRODUCTION_NAME, STEP_FILEBASENAME, PICO_DIR, NANOAOD_VERSION, FIRST_COMMAND, NO_RUN, data_tag)

if __name__ == '__main__':
  PRODUCTION_NAME = 'higgsino_inyo'
  PICO_DIR = '/net/cms25/cms25r5/pico'
  NANOAOD_VERSION = 'NanoAODv7'
  notify_script = 'sendTelegramMessage.py'
  # Use sendMail if telegram is not setup
  #notify_script = os.path.dirname(os.path.realpath(__file__))+'/sendEmail.sh'

  # Do not run commands. Used for debugging
  NO_RUN = False
  # To prompt for first command
  FIRST_COMMAND = True

  # Main steps
  # - Process
  # - Skim
  # - Slim and Skim
  # There is a disk list of commands. There is a mem list of commands. Compare and run missing steps.

  # Step is the step that is running. Will run the next step.
  #FIRST_COMMAND = processMc(YEAR=2016, PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME='process_steps_higgsino.py.'+PRODUCTION_NAME+'.2016.mc.step', PICO_DIR=PICO_DIR, NANOAOD_VERSION=NANOAOD_VERSION, FIRST_COMMAND=FIRST_COMMAND, notify_script=notify_script)
  #FIRST_COMMAND = processSignal(YEAR=2016, PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME='process_steps_higgsino.py.'+PRODUCTION_NAME+'.2016.sig.step', PICO_DIR=PICO_DIR, NANOAOD_VERSION=NANOAOD_VERSION, FIRST_COMMAND=FIRST_COMMAND, notify_script=notify_script)
  FIRST_COMMAND = processData(YEAR=2016, PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME='process_steps_higgsino.py.'+PRODUCTION_NAME+'.2016.data.step', PICO_DIR=PICO_DIR, NANOAOD_VERSION=NANOAOD_VERSION, FIRST_COMMAND=FIRST_COMMAND, notify_script=notify_script)

  #FIRST_COMMAND = processMc(YEAR=2017, PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME='process_steps_higgsino.py.'+PRODUCTION_NAME+'.2017.mc.step', PICO_DIR=PICO_DIR, NANOAOD_VERSION=NANOAOD_VERSION, FIRST_COMMAND=FIRST_COMMAND, notify_script=notify_script)
  #FIRST_COMMAND = processSignal(YEAR=2017, PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME='process_steps_higgsino.py.'+PRODUCTION_NAME+'.2017.sig.step', PICO_DIR=PICO_DIR, NANOAOD_VERSION=NANOAOD_VERSION, FIRST_COMMAND=FIRST_COMMAND, notify_script=notify_script)
  FIRST_COMMAND = processData(YEAR=2017, PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME='process_steps_higgsino.py.'+PRODUCTION_NAME+'.2017.data.step', PICO_DIR=PICO_DIR, NANOAOD_VERSION=NANOAOD_VERSION, FIRST_COMMAND=FIRST_COMMAND, notify_script=notify_script)

  #FIRST_COMMAND = processMc(YEAR=2018, PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME='process_steps_higgsino.py.'+PRODUCTION_NAME+'.2018.mc.step', PICO_DIR=PICO_DIR, NANOAOD_VERSION=NANOAOD_VERSION, FIRST_COMMAND=FIRST_COMMAND, notify_script=notify_script)
  #FIRST_COMMAND = processSignal(YEAR=2018, PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME='process_steps_higgsino.py.'+PRODUCTION_NAME+'.2018.sig.step', PICO_DIR=PICO_DIR, NANOAOD_VERSION=NANOAOD_VERSION, FIRST_COMMAND=FIRST_COMMAND, notify_script=notify_script)
  FIRST_COMMAND = processData(YEAR=2018, PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME='process_steps_higgsino.py.'+PRODUCTION_NAME+'.2018.data.step', PICO_DIR=PICO_DIR, NANOAOD_VERSION=NANOAOD_VERSION, FIRST_COMMAND=FIRST_COMMAND, notify_script=notify_script)
