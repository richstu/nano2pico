#!/usr/bin/env python3

import os, argparse
from glob import glob

def getTag(path):
  tag = path.split("/")[-1]
  tag = tag.split("RunIISummer16NanoAODv4")[0]
  tag = tag.split("RunIISummer16NanoAODv5")[0]
  tag = tag.split("_ext")[0]
  tag = tag.replace("raw_pico_","")
  tag = tag.strip("_")
  return tag

parser = argparse.ArgumentParser(description="Submits batch jobs to apply new SFs and compute sum-of-weights",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i","--in_dir", default="/net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/raw_pico/",
                    help="Directory where the raw pico files are")
parser.add_argument("-o","--out_cmd_file", default="cmds.py",
                    help="File with list of commands for batch system.")
parser.add_argument("--overwrite", default=False,
                    help="Process all input files regardless whether output exists.")
args = parser.parse_args()

in_dir = args.in_dir
if 'raw_pico' not in in_dir:
  print('Input directory expected to contain string "/raw_pico/"!')
  sys.exit(1)
out_dir = args.in_dir.replace('/raw_pico/','/unskimmed/')
if not os.path.exists(out_dir): 
  os.mkdir(out_dir)

in_file_paths = glob(os.path.join(in_dir,'*.root'))

nfiles = 0
cmdfile = open(args.out_cmd_file,'w')
cmdfile.write('#!/bin/env python\n')
for ifile_path in in_file_paths:
  outfile_path = ifile_path.replace('/raw_pico/','/unskimmed/').replace('/raw_pico_','/pico_')
  if not args.overwrite and os.path.exists(outfile_path):
    continue
  ifile = os.path.basename(os.path.realpath(ifile_path))
  corr_file = 'corr_'+getTag(ifile)+'.root'
  cmd = '{}/run/apply_corrections.exe -f {} -i {} -c {}'.format(os.getcwd(), ifile, in_dir, corr_file)
  cmdfile.write('print(\"'+cmd+'\")\n')
  nfiles +=1

cmdfile.close()
os.chmod(args.out_cmd_file, 0o755)

print("Found {} input files.\n".format(nfiles))

print("To generate job json and submit jobs do: ")
print('convert_cl_to_jobs_info.py '+args.out_cmd_file+' apply_corrections.json')
print('auto_submit_jobs.py apply_corrections.json -c scripts/check_apply_corrections_job.py')
  
