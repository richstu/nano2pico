#!/usr/bin/env python
# The ./jobscript_check.py should return 'success' or 'fail' or 'to_submit' or 'submitted' for a job_log_string
# The input is given as sys.argv[1] = queue_system.compress_string(job_log_string) sys.argv[2] = queue_system.compress_string(job_argument_string)
import sys
import os.path
import queue_system
from ROOT import TFile

def get_era(run):
  if (run >= 272007 and run <= 275376):
    return 'Run2016B'
  elif (run >= 275657 and run <= 276283):
    return 'Run2016C'
  elif (run >= 276315 and run <= 276811):
    return 'Run2016D'
  elif (run >= 276831 and run <= 277420):
    return 'Run2016E'
  elif (run >= 277772 and run <= 278808):
    return 'Run2016F'
  elif (run >= 278820 and run <= 280385):
    return 'Run2016G'
  elif (run >= 280919 and run <= 284044):
    return 'Run2016H'
  else:
    return 'Unknown'

# job_argument_string = "/net/top/homes/aovcharova/code/nano2pico/run/process_nano.exe -f DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv5__PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1__250000__233E4358-D599-FE4A-B585-A6B18F4DDEF1.root -i /mnt/hadoop/pico/NanoAODv5/nano/2016/mc/ -o /net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/"
# job argument_string = "net/top/homes/oshiro/code/nano2pico/run/combine_datasets.exe -i /net/cms29/cms29r0/pico/NanoAODv5/ttz_cordellbankv2/2016/data/raw_pico/ -o /net/cms29/cms29r0/pico/NanoAODv5/ttz_cordellbankv2/2016/data/raw_pico//../unskimmed -f ./txt/datasets/ttz_2016_data_datasetcombine_list.txt -b 276811 -e 276811 -y Run2016 -n Nano1June2019-v1"
job_log_string = queue_system.decompress_string(sys.argv[1])
job_argument_string = queue_system.decompress_string(sys.argv[2])

print(job_argument_string)

args = job_argument_string.split('--command="')[1].split('"')[0]
tmp = args.split(' ')
infile_path = tmp[4]+'/'+tmp[2]
outfile_path = tmp[6]+'/raw_pico/raw_pico_'+tmp[2]

#load appropriate info from datasets file
datasets_file = open(tmp[6])
datasets_list = datasets_file.read().split("\n")
datasets_list.pop()
sample_names = ""
first_name = True
for dataset_name in datasets_list:
  if first_name:
    sample_names += dataset_name
    first_name = False
  else:
    sample_names += '_'
    sample_names += dataset_name
sample_number = len(datasets_list)

output_filenames = []
era_name = get_era(int(tmp[8]))
if era_name != 'Unknown':
  sample_index = 0
  while sample_index < sample_number:
    output_filenames.append(tmp[4] + '/pico_' + era_name + '_' + str(sample_index) + '_' + sample_names + '_' + tmp[14] + '_runs' + tmp[8] + '.root')
    sample_index += 1
else:
  sample_index = 0
  while sample_index < sample_number:
    output_filenames.append(tmp[4] + '/pico_' + tmp[12] + '_' + str(sample_index) + '_' + sample_names + '_' + tmp[14] + '_runs' + tmp[8] + '.root')
    sample_index += 1

job_failed = False
fail_message = ''
for output_filename in output_filenames:
  if os.path.isfile(output_filename):
    #file exists, check that is was closed properly
    output_file = TFile(output_filename)
    if output_file.TestBit(TFile.kRecovered):
      #file had to be recovered, job must have crashed
      fail_message = 'File not closed properly.'
      job_failed = True
      break
  else:
    #file doesn't even exist, job must have crashed
    fail_message = 'File doesn\'t exist.'
    job_failed = True
    break

if job_failed:
  print('[For queue_system] fail: ' + fail_message)
else:
  print('[For queue_system] success')
