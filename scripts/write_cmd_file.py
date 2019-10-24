#!/usr/bin/env python3

import os, argparse
from glob import glob

parser = argparse.ArgumentParser(description="Submits batch jobs to apply new SFs and compute sum-of-weights",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i","--in_dir", default="/mnt/hadoop/pico/NanoAODv5/nano/2016/mc/",
                    help="Directory where the NanoAOD files are")
parser.add_argument("-p","--production", default="higgsino_angeles",
                    help="Determines the output folder.")
parser.add_argument("-d","--datasets_file", default="2016_mc_dataset_paths",
                    help="File with the list of dataset names as they appear in DAS")
parser.add_argument("-o","--out_cmd_file", default="cmds.py",
                    help="File with list of commands for batch system.")
args = parser.parse_args()

in_dir = args.in_dir
wgt_dir = os.path.join(args.in_dir.replace('/nano/','/'+args.production+'/'),'wgt_sums')
out_dir = os.path.join(args.in_dir.replace('/nano/','/'+args.production+'/'),'reweighted')

if not os.path.exists(wgt_dir): os.mkdir(wgt_dir)
if not os.path.exists(out_dir): os.mkdir(out_dir)

datasets = []
with open(args.datasets_file) as f:
  datasets = f.readlines().split()


cmdfile = open(args.out_cmd_file,'w')
cmdfile.write('#!/bin/env python')

for ds in datasets:
  ds = ds.split("/")[1]
  infile_paths = glob(os.path.join(in_dir,'*'+ds+'*.root'))
  for infile_path in infile_paths:
    infile = os.path.basename(os.path.realpath(infile))
    wfile_path = os.path.join(wgt_dir,'wgt_sums_'+infile)
    outfile_path = os.path.join(out_dir,'pico_'+infile)
    cmd = '{}/run/process_nano.exe -i {} -w {} -o {}'.format(os.getcwd(), infile_path, wfile_path, outfile_path)
    cmdfile.write('print(\"'+cmd+'\")\n')

cmdfile.close()
os.chmod(cmdfile, 0o755)
os.cmd('convert_cl_to_jobs_info.py '+cmdfile+' '+production+'.json')

print("To submit jobs do: ")
print('auto_submit_jobs.py '+production+'.json')
  
