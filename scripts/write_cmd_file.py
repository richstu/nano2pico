#!/usr/bin/env python3

import os, argparse
from glob import glob

parser = argparse.ArgumentParser(description="Submits batch jobs to apply new SFs and compute sum-of-weights",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i","--in_dir", default="/mnt/hadoop/pico/NanoAODv5/nano/2016/mc/",
                    help="Directory where the NanoAOD files are")
parser.add_argument("-p","--production", default="higgsino_angeles",
                    help="Determines the output folder.")
parser.add_argument("-d","--datasets_file", default="datasets/higgsino_2016_mc_dataset_paths.txt",
                    help="File with the list of dataset names as they appear in DAS")
parser.add_argument("-o","--out_cmd_file", default="cmds.py",
                    help="File with list of commands for batch system.")
args = parser.parse_args()

in_dir = args.in_dir
out_base_dir = args.in_dir.replace('/mnt/hadoop/','/net/cms29/cms29r0/').replace('/nano/','/'+args.production+'/')

if not os.path.exists(out_base_dir): os.mkdir(out_base_dir)

datasets = []
with open(args.datasets_file) as f:
  datasets = f.readlines()


cmdfile = open(args.out_cmd_file,'w')
cmdfile.write('#!/bin/env python\n')

for ds in datasets:
  tmp_ = ds.split("/")
  if len(tmp) > 1: # in case of empty lines in file
    ds = tmp_[1]
  in_file_paths = glob(os.path.join(in_dir,'*'+ds+'*.root'))
  for ifile_path in in_file_paths:
    ifile = os.path.basename(os.path.realpath(ifile_path))
    cmd = '{}/run/process_nano.exe -f {} -i {} -o {}'.format(os.getcwd(), ifile, in_dir, out_base_dir)
    cmdfile.write('print(\"'+cmd+'\")\n')

cmdfile.close()
os.chmod(args.out_cmd_file, 0o755)

print("To generate job json and submit jobs do: ")
print('convert_cl_to_jobs_info.py '+cmdfile+' '+production+'.json')
print('auto_submit_jobs.py '+production+'.json')
  
