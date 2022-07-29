#!/usr/bin/env python
# The ./jobscript_check.py should return 'success' or 'fail' or 'to_submit' or 'submitted' for a job_log_string
# The input is given as sys.argv[1] = queue_system.compress_string(job_log_string) sys.argv[2] = queue_system.compress_string(job_argument_string)
import sys
import queue_system
import ROOT
import argparse
import shlex

# Note that hyphens get removed for output.
def get_args(keys, job_argument_string):
  parser = argparse.ArgumentParser()
  for key in keys:
    parser.add_argument(key)
  args, unknown_args = parser.parse_known_args(shlex.split(job_argument_string))
  return vars(args)

if __name__ == "__main__":
  job_log_string = queue_system.decompress_string(sys.argv[1])
  job_argument_string = queue_system.decompress_string(sys.argv[2])
  
  print(job_argument_string)

  job_args = get_args(['--command'], job_argument_string)

  # command = "slim_and_merge.py -s txt/slim_rules/zgdata.txt -o /merged_zgdata_llg/merged_raw_pico_llg_DoubleEG__Run2016CC_zgdata_llg_nfiles_24.root -i /skim_llg/*raw_pico_llg_DoubleEG__Run2016C*.root"
  command_args = get_args(['-s', '-o', '-i'], job_args['command'])
  print(command_args)
  
  unslimmed_path = command_args['i']
  slimmed_path = command_args['o']
  
  unslimmed = ROOT.TChain("tree")
  unslimmed.Add(unslimmed_path)
  unslimmed_events = unslimmed.GetEntries()
  
  slimmed = ROOT.TChain("tree")
  slimmed.Add(slimmed_path)
  slimmed_events = slimmed.GetEntries()
  
  print(unslimmed_path+' has '+str(unslimmed_events)+' events')
  print(slimmed_path+' has '+str(slimmed_events)+' events')
  
  if unslimmed_events != slimmed_events:
    print('[For queue_system] fail: unslimmed and slimmed do not have same number of events')
    sys.exit(1)
  else:
    print('[For queue_system] success')
