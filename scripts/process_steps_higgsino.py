#!/usr/bin/env python
import os
import glob
import sys
# TODO Check before and after for production with cutflow

def process(YEAR, PRODUCTION_NAME, STEP_FILEBASENAME, PICO_DIR, DO_SIGNAL, FIRST_COMMAND):
  notify_script = 'sendTelegramMessage.py'
  YEAR = str(YEAR)
  mc_tag=PRODUCTION_NAME+'_'+YEAR+'_mc'
  # Add mc commands
  process_commands = [
    #0
    [notify_script+' "Start process nano '+mc_tag+'"',
    './scripts/write_process_nano_cmds.py --in_dir '+PICO_DIR+'/NanoAODv5/nano/'+YEAR+'/mc/ --production '+PRODUCTION_NAME+' --dataset_list datasets/mc_dataset_paths --tag '+mc_tag,
    'auto_submit_jobs.py process_nano_cmds_'+mc_tag+'.json -c scripts/check_process_nano_job.py -f',
    notify_script+' "Finished process nano '+mc_tag+'"'], 
    
    #1
    [notify_script+' "Start merge corrections '+mc_tag+'"',
    './scripts/merge_corrections.py --wgt_dir '+PICO_DIR+'/NanoAODv5/'+PRODUCTION_NAME+'/'+YEAR+'/mc/wgt_sums/ --corr_dir '+PICO_DIR+'/NanoAODv5/'+PRODUCTION_NAME+'/'+YEAR+'/mc/corrections/',
    notify_script+' "Finished merge corrections '+mc_tag+'"'],
    
    #2
    [notify_script+' "Start applied corrections '+mc_tag+'"',
    './scripts/write_apply_corrections_cmds.py --in_dir '+PICO_DIR+'/NanoAODv5/'+PRODUCTION_NAME+'/'+YEAR+'/mc/raw_pico/ --tag '+mc_tag,
    'auto_submit_jobs.py '+mc_tag+'_apply_corr_cmds.json -c scripts/check_apply_corrections_job.py -f',
    notify_script+' "Finished applied corrections '+mc_tag+'"'],

    #3
    [notify_script+' "Start skim met150 '+mc_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/NanoAODv5/'+PRODUCTION_NAME+'/'+YEAR+'/mc/unskimmed/ --skim_name met150 --tag '+mc_tag,
    'auto_submit_jobs.py skim_met150_cmds_'+mc_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim met150 '+mc_tag+'"'],
    
    #4
    [notify_script+' "Start merge met150 '+mc_tag+'"',
    './scripts/write_slim_and_merge_cmds.py --in_dir '+PICO_DIR+'/NanoAODv5/'+PRODUCTION_NAME+'/'+YEAR+'/mc/skim_met150/ --slim_name higmc --tag '+mc_tag,
    'auto_submit_jobs.py '+mc_tag+'_slim_higmc_met150_cmds.json -f',
    './scripts/confirm_slim.py '+PICO_DIR+'/NanoAODv5/'+PRODUCTION_NAME+'/'+YEAR+'/mc/skim_met150 '+PICO_DIR+'/NanoAODv5/'+PRODUCTION_NAME+'/'+YEAR+'/mc/merged_higmc_met150',
    notify_script+' "Finished merge met150 '+mc_tag+'"'],

    #5
    [notify_script+' "Start skim preselect '+mc_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/NanoAODv5/'+PRODUCTION_NAME+'/'+YEAR+'/mc/merged_higmc_met150/ --skim_name preselect --tag '+mc_tag,
    'auto_submit_jobs.py skim_preselect_cmds_'+mc_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim preselect '+mc_tag+'"'],
    
    #6
    [notify_script+' "Start skim higloose '+mc_tag+'"',
    './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/NanoAODv5/'+PRODUCTION_NAME+'/'+YEAR+'/mc/merged_higmc_met150/ --skim_name higloose --tag '+mc_tag,
    'auto_submit_jobs.py skim_higloose_cmds_'+mc_tag+'.json -c scripts/check_skim.py -f',
    notify_script+' "Finished skim higloose '+mc_tag+'"'],
  ]

  # Add signal commands
  if DO_SIGNAL:
    sig_tag=PRODUCTION_NAME+'_'+YEAR+'_sig'
    process_commands_sig = [
      # signal
      #7
      [notify_script+' "Started process nano '+sig_tag+'"',
      './scripts/write_process_nano_cmds.py --in_dir '+PICO_DIR+'/NanoAODv5/nano/'+YEAR+'/SMS-TChiHH_2D/ --production '+PRODUCTION_NAME+' --tag '+sig_tag,
      'auto_submit_jobs.py process_nano_cmds_'+sig_tag+'.json -c scripts/check_process_nano_job.py -f',
      notify_script+' "Finished process nano '+sig_tag+'"'],
      
      #8
      [notify_script+' "Started merge corrections '+sig_tag+'"',
      './scripts/merge_corrections.py --wgt_dir '+PICO_DIR+'/NanoAODv5/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/wgt_sums/ --corr_dir '+PICO_DIR+'/NanoAODv5/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/corrections/',
      notify_script+' "Finished merge corrections '+sig_tag+'"'],
      
      #9
      [notify_script+' "Started applied corrections '+sig_tag+'"',
      './scripts/write_apply_corrections_cmds.py --in_dir '+PICO_DIR+'/NanoAODv5/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/raw_pico/ --tag '+sig_tag,
      'auto_submit_jobs.py '+sig_tag+'_apply_corr_cmds.json -c scripts/check_apply_corrections_job.py -f',
      notify_script+' "Finished applied corrections '+sig_tag+'"'],
      
      #10
      [notify_script+' "Started skim met150 '+sig_tag+'"',
      './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/NanoAODv5/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/unskimmed/ --skim_name met150 --tag '+sig_tag,
      'auto_submit_jobs.py skim_met150_cmds_'+sig_tag+'.json -c scripts/check_skim.py -f',
      notify_script+' "Finished skim met150 '+sig_tag+'"'],
      
      #11
      [notify_script+' "Started merge met150 '+sig_tag+'"',
      './scripts/write_slim_and_merge_cmds.py --in_dir '+PICO_DIR+'/NanoAODv5/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/skim_met150/ --slim_name higmc --tag '+sig_tag,
      'auto_submit_jobs.py '+sig_tag+'_slim_higmc_met150_cmds.json -f',
      './scripts/confirm_slim.py '+PICO_DIR+'/NanoAODv5/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/skim_met150 '+PICO_DIR+'/NanoAODv5/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/merged_higmc_met150',
      notify_script+' "Finished merge met150 '+sig_tag+'"'],
      
      #12
      [notify_script+' "Started skim preselect '+sig_tag+'"',
      './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/NanoAODv5/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/merged_higmc_met150 --skim_name preselect --tag '+sig_tag,
      'auto_submit_jobs.py skim_preselect_cmds_'+sig_tag+'.json -c scripts/check_skim.py -f',
      notify_script+' "Finished skim preselect '+sig_tag+'"'],

      #13
      [notify_script+' "Started skim higloose '+sig_tag+'"',
      './scripts/write_skim_cmds.py --in_dir '+PICO_DIR+'/NanoAODv5/'+PRODUCTION_NAME+'/'+YEAR+'/SMS-TChiHH_2D/merged_higmc_met150 --skim_name higloose --tag '+sig_tag,
      'auto_submit_jobs.py skim_higloose_cmds_'+sig_tag+'.json -c scripts/check_skim.py -f',
      notify_script+' "Finished skim higloose '+sig_tag+'"'],

    ]
    process_commands.extend(process_commands_sig)

  # Search for current step
  search = glob.glob(STEP_FILEBASENAME+"*")
  if len(search) == 1:
    STEP_FILENAME = search[0]
    step = int(STEP_FILENAME.replace(STEP_FILEBASENAME,''))
  else:
    step = -1
  next_step = step + 1
  if next_step >= len(process_commands):
    print('Number of commands is '+str(len(process_commands)))
    print('Log file from index 0: '+STEP_FILENAME)
    print('All commands have run for '+PRODUCTION_NAME+"_"+YEAR)
    return

  # Print commands
  if FIRST_COMMAND:
    if step != -1:
      print('[Step '+str(step)+'] Previous commands')
      print('\n'.join(process_commands[step]))
      print('\n')
    print('[Step '+str(next_step)+'] Next command to run')
    print('\n'.join(process_commands[next_step]))
    raw_input('Press Enter to continue')
    FIRST_COMMAND = False

  # Do steps
  for iRemainStep in range(len(process_commands)-next_step):
    # Run commands
    error = 0
    for command in process_commands[next_step]:
      print('Will run below command')
      print('  '+command)
      error += os.system(command)
    # Log
    if next_step != 0:
      os.system('rm '+STEP_FILEBASENAME+str(next_step-1))
    os.system('touch '+STEP_FILEBASENAME+str(next_step))

    # Check if there was an error
    if error != 0:
      print('There was an error')
      os.system(notify_script+' "There was an error on step '+str(next_step)+' for '+PRODUCTION_NAME+'_'+YEAR+'"')
      sys.exit(1)

    next_step += 1
  return FIRST_COMMAND

if __name__ == '__main__':
  PRODUCTION_NAME = 'higgsino_humboldt'
  PICO_DIR = '/net/cms25/cms25r5/jbkim/pico'

  # To prompt for first command
  FIRST_COMMAND = True

  # Step is the step that is running. Will run the next step.
  FIRST_COMMAND = process(YEAR=2016, PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME='process_steps_higgsino.py.'+PRODUCTION_NAME+'.2016.step', PICO_DIR=PICO_DIR, DO_SIGNAL=True, FIRST_COMMAND=FIRST_COMMAND)
  FIRST_COMMAND = process(YEAR=2017, PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME='process_steps_higgsino.py.'+PRODUCTION_NAME+'.2017.step', PICO_DIR=PICO_DIR, DO_SIGNAL=False, FIRST_COMMAND=FIRST_COMMAND)
  FIRST_COMMAND = process(YEAR=2018, PRODUCTION_NAME=PRODUCTION_NAME, STEP_FILEBASENAME='process_steps_higgsino.py.'+PRODUCTION_NAME+'.2018.step', PICO_DIR=PICO_DIR, DO_SIGNAL=False, FIRST_COMMAND=FIRST_COMMAND)
