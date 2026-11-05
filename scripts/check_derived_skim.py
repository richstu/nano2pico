#!/usr/bin/env python3
# Check script for skim_file.py jobs of the form
#   skim_file.py -k <skim> -i <infile> -o <outdir>/ [--from_skim <skim>]
# Like check_direct_to_skim.py, the expected number of entries is taken from the
# job log (SKIM_PASS_EVENTS: N, printed by skim_file.py), so the input file is not re-read.
# argv[1] = compressed job log, argv[2] = compressed job argument string.

import sys
import os
import re
import zlib
from pico_count import count_entries


def decompress_string(compressed_string):
  return zlib.decompress(bytes.fromhex(compressed_string)).decode('utf-8')


job_log_string = decompress_string(sys.argv[1])
job_argument_string = decompress_string(sys.argv[2])

args = job_argument_string.split('--command="')[1].split('"')[0]
tmp = args.split(' ')
# tmp: [0]=skim_file.py [1]=-k [2]=skim_name [3]=-i [4]=infile [5]=-o [6]=outdir [7:]=--from_skim X
skim_name = tmp[2]
infile_path = tmp[4]
out_dir = tmp[6]
in_name = os.path.basename(os.path.realpath(infile_path))
if '--from_skim' in tmp:
  from_prefix = 'pico_'+tmp[tmp.index('--from_skim')+1]+'_'
  if not in_name.startswith(from_prefix):
    print('[For queue_system] fail: input name {} does not start with {}'.format(in_name, from_prefix))
    sys.exit(0)
  out_name = 'pico_'+skim_name+'_'+in_name[len(from_prefix):]
else:
  out_name = in_name.replace('pico_', 'pico_'+skim_name+'_')
outfile_path = os.path.join(out_dir, out_name)
print('outfile: ', outfile_path)

low = job_log_string.lower()
if 'traceback (most recent call last)' in low or 'segmentation fault' in low or 'error in <' in low:
  print('[For queue_system] fail: Error in job_log')
  sys.exit(0)

match = re.search(r'SKIM_PASS_EVENTS:\s*(\d+)', job_log_string)
if match is None:
  print('[For queue_system] fail: Could not find SKIM_PASS_EVENTS in job log.')
  sys.exit(0)
in_nent = int(match.group(1))

if not os.path.exists(outfile_path):
  print('[For queue_system] fail: output ({}) does not exist.'.format(outfile_path))
  sys.exit(0)
out_nent, nbranches = count_entries(outfile_path, 'tree')
print('DEBUG: expected entries = {}, output entries = {}'.format(in_nent, out_nent))
if nbranches == 0:
  print('[For queue_system] fail: output ({}) has no branches.'.format(outfile_path))
elif in_nent == out_nent:
  print('[For queue_system] success')
else:
  print('[For queue_system] fail: Expected {} entries from log, while output ({}) has {} entries.'.format(in_nent, outfile_path, out_nent))

