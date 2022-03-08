#!/usr/bin/env python3

import os, argparse
from glob import glob

def getTag(path):
  tag = path.split("/")[-1]
  tag = tag.split("RunIISummer16NanoAODv4")[0]

  tag = tag.split("RunIISummer16NanoAODv5")[0]
  tag = tag.split("RunIIFall17NanoAODv5")[0]
  tag = tag.split("RunIIAutumn18NanoAODv5")[0]

  tag = tag.split("RunIISummer16NanoAODv7")[0]
  tag = tag.split("RunIIFall17NanoAODv7")[0]
  tag = tag.split("RunIIAutumn18NanoAODv7")[0]

  tag = tag.split("RunIISummer20UL16NanoAODv9")[0]
  tag = tag.split("RunIISummer20UL17NanoAODv9")[0]
  tag = tag.split("RunIISummer20UL18NanoAODv9")[0]

  tag = tag.split("_ext")[0]
  tag = tag.replace("raw_pico_","")
  tag = tag.strip("_")
  return tag

parser = argparse.ArgumentParser(description="Submits batch jobs to apply new SFs and compute sum-of-weights",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i","--in_dir", default="/net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/raw_pico/",
                    help="Directory where the raw pico files are")
parser.add_argument("--overwrite", action='store_true',
                    help="Process all input files regardless whether output exists.")
parser.add_argument('-t', '--tag', default='',
                    help='Optionally specify a tag to be used to differentiate helper files for batch submission.')
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
cmdfile_name = 'apply_corr_cmds.py'
if (args.tag!=''): cmdfile_name = args.tag+'_'+cmdfile_name
cmdfile = open(cmdfile_name,'w')
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
os.chmod(cmdfile_name, 0o755)

print("Found {} input files.\n".format(nfiles))

json_name = cmdfile_name.replace('.py','.json')
print('To print a sample command:')
print('cat '+cmdfile_name+' | tail -n 1\n')
print("To generate job json and submit jobs do: ")
os.system('convert_cl_to_jobs_info.py -f '+cmdfile_name+' '+json_name)
print('auto_submit_jobs.py '+json_name+' -c scripts/check_apply_corrections_job.py')
  
