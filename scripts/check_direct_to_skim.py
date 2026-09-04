#!/usr/bin/env python3
# The ./jobscript_check.py should return 'success' or 'fail' or 'to_submit' or 'submitted' for a job_log_string
# The input is given as: sys.argv[1] = queue_system.compress_string(job_log_string) sys.argv[2] = queue_system.compress_string(job_argument_string)

import sys
import os
import re
import zlib
from ROOT import TChain


# Decompress string for passed things through posix
# this is normally in queue system, but it breaks with python3 so here is a
# reimplementation for now
def decompress_string(compressed_string):
  return zlib.decompress(bytes.fromhex(compressed_string)).decode('utf-8')


job_log_string = decompress_string(sys.argv[1])
job_argument_string = decompress_string(sys.argv[2])

print(job_argument_string)

args = job_argument_string.split('--command="')[1].split('"')[0]
tmp = args.split(' ')

infile_name = tmp[2]
in_dir = tmp[4]
out_dir = tmp[6]

infile_path = os.path.join(in_dir, infile_name)
outfile_path = os.path.join(out_dir,'skim_'+tmp[10],'pico_'+tmp[10]+'_'+infile_name)
print("infile: ", infile_path)
print("outfile: ", outfile_path)

if ('error' in job_log_string.lower()) or ('segmentation fault' in job_log_string.lower()):
  print('[For queue_system] fail: Error in job_log')
  sys.exit(0)

# Extract expected number of skimmed events from log
match = re.search(r'SKIM_PASS_EVENTS:\s*(\d+)', job_log_string)

if match is None:
  print('[For queue_system] fail: Could not find SKIM_PASS_EVENTS in job log.')
  sys.exit(0)

in_nent = int(match.group(1))


# Check output file
if os.path.exists(outfile_path):
  outfile = TChain("tree")
  outfile.Add(outfile_path)
  out_nent = outfile.GetEntries()
  nbranches = outfile.GetNbranches()
else:
    out_nent = 0
    nbranches = 0

print('DEBUG: expected entries = {}, output entries = {}'.format(in_nent, out_nent))

if os.path.exists(outfile_path) and nbranches == 0:
  print('[For queue_system] fail: output ({}) has no branches.'.format(outfile_path))
elif in_nent == out_nent:
    print('[For queue_system] success')
else:
  print('[For queue_system] fail: Expected {} entries from log, while output ({}) has {} entries.'.format(in_nent, outfile_path, out_nent))
