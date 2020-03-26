#!/usr/bin/env python
# The ./jobscript_check.py should return 'success' or 'fail' or 'to_submit' or 'submitted' for a job_log_string
# The input is given as sys.argv[1] = queue_system.compress_string(job_log_string) sys.argv[2] = queue_system.compress_string(job_argument_string)
import sys
import queue_system
from ROOT import TChain

# job_argument_string = "/net/top/homes/aovcharova/code/nano2pico/run/process_nano.exe -f DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv5__PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1__250000__233E4358-D599-FE4A-B585-A6B18F4DDEF1.root -i /mnt/hadoop/pico/NanoAODv5/nano/2016/mc/ -o /net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/"
job_log_string = queue_system.decompress_string(sys.argv[1])
job_argument_string = queue_system.decompress_string(sys.argv[2])

print(job_argument_string)

args = job_argument_string.split('--command="')[1].split('"')[0]
tmp = args.split(' ')
infile_path = tmp[4]+'/'+tmp[2]
outfile_path = tmp[6]+'/raw_pico/raw_pico_'+tmp[2]

infile = TChain("Events");
infile.Add(infile_path);
in_nent = infile.GetEntries()

outfile = TChain("tree");
outfile.Add(outfile_path);
out_nent = outfile.GetEntries()

if outfile.GetNbranches() == 0:
  print('[For queue_system] fail: output ({}) has no branches.'.format(outfile_path))
elif in_nent == out_nent:
  print('[For queue_system] success')
else:
  print('[For queue_system] fail: Input ({}) has {} entries, while output ({}) has {} entries.'.format(infile_path, in_nent, outfile_path, out_nent))
