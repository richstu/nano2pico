#!/usr/bin/env python
# The ./jobscript_check.py should return 'success' or 'fail' or 'to_submit' or 'submitted' for a job_log_string
# The input is given as sys.argv[1] = queue_system.compress_string(job_log_string) sys.argv[2] = queue_system.compress_string(job_argument_string)
import sys
import os.path
import queue_system
from ROOT import TFile

# job_argument_string = "/net/top/homes/aovcharova/code/nano2pico/run/process_nano.exe -f DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv5__PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1__250000__233E4358-D599-FE4A-B585-A6B18F4DDEF1.root -i /mnt/hadoop/pico/NanoAODv5/nano/2016/mc/ -o /net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/"
job_log_string = queue_system.decompress_string(sys.argv[1])
job_argument_string = queue_system.decompress_string(sys.argv[2])

print(job_argument_string)

args = job_argument_string.split('--command="')[1].split('"')[0]
tmp = args.split(' ')
infile_path = tmp[4]+'/'+tmp[2]
outfile_path = tmp[6]+'/raw_pico/raw_pico_'+tmp[2]

job_failed = False
if os.path.isfile(outfile_path):
  #file exists, check that is was closed properly
  output_file = TFile(outfile_path)
  if output_file.TestBit(TFile.kRecovered):
    #file had to be recovered, job must have crashed
    fail_message = 'File not closed properly.'
    job_failed = True
  if output_file.IsZombie():
    #file is zombie, job crashed
    fail_message = 'File not closed properly, is now zombie.'
    job_failed = True
else:
  #file doesn't even exist, job must have crashed
  fail_message = 'File doesn\'t exist.'
  job_failed = True

if job_failed:
  print('[For queue_system] fail: ' + fail_message)
else:
  print('[For queue_system] success')
