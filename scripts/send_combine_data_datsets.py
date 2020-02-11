#!/usr/bin/env python

# Script to send jobs to merge ntuples
# Heavily modified from version in richstu/babymaker

import os, sys, subprocess, argparse
from glob import glob
import json
import re
import string

parser = argparse.ArgumentParser(description="Submits batch jobs to remove overlapping events from pico data datasets.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i","--in_dir", default="",
                    help="Directory where the pico files are, e.g. /net/cms29/cms29r0/pico/NanoAODv5/ttz_cordellbank/2016/data/raw_pico")
parser.add_argument("-o","--out_dir", default="",
                    help="Determines the output folder. If left empty, will use in_dir/../unskimmed.")
parser.add_argument("-f","--out_cmd_file", default="cmds.py",
                    help="File with list of commands for batch system.")
parser.add_argument("-d","--dataset_list", default="txt/datasets/dataset_list.txt",
                    help="File with the list of dataset names")
args = parser.parse_args()
in_dir = args.in_dir
out_dir = args.out_dir
dataset_list = args.dataset_list
if len(out_dir)==0:
    out_dir = args.in_dir+"/../unskimmed"
if not os.path.exists(out_dir): 
    os.makedirs(out_dir)
if not os.path.isfile(dataset_list):
    print("Dataset file invalid.")
    sys.exit(0)

#figure out year from input folder and get appropriate golden json
year_match = re.search('20\d\d', in_dir)
if not year_match:
    print("Could not determine year from input directory")
    sys.exit(0)
else:
    year = year_match.group(0)
jsonfile = ''
if year=="2016":
    jsonfile = 'txt/json/golden_Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16.json'
elif year=="2017":
    jsonfile = 'txt/json/golden_Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17.json'
elif year=="2018":
    jsonfile = 'txt/json/golden_Cert_314472-325175_13TeV_PromptReco_Collisions18.json'
else:
    print("No golden json available for specified year.")
    sys.exit(0)

#get NanoAOD production from file names
all_input_files = glob(os.path.join(in_dir,'*.root'))
print("Found {} input files.\n".format(len(all_input_files)))
nano_version = ""
if len(all_input_files) > 0:
    fname_parts = (all_input_files[0].split("/"))[-1].split("__")
    nano_version = ""
    for fname_part in fname_parts:
        if (fname_part[0:4] == "Nano"):
            nano_version = fname_part
    if (len(fname_part) == 0):
        print("Unknown Nano/pico naming convention, this script may need to be modified.")
        sys.exit(0)
    nano_version = string.replace(nano_version,"_ver1","")
    nano_version = string.replace(nano_version,"_ver2","")
    nano_version = string.replace(nano_version,"_ver3","")
    nano_version = string.replace(nano_version,"_ver4","")
else:
    print("No files in input directory. Exiting.")
    sys.exit(0)

# Number of runs in each ntuple
runs_file = 1 
runs=[]
with open(jsonfile) as jfile:
  for line in jfile:
    for word in line.split():
      if '"' in word: 
        word = word.split('"')[1]
        runs.append(word)
# Dividing runs into sets of "runs_file" elements
runs = [runs[i:i+runs_file] for i in xrange(0, len(runs), runs_file)]
print('Number of runs: {}\n'.format(len(runs)))

# write python file for running commands
cmdfile = open(args.out_cmd_file,'w')
cmdfile.write('#!/bin/env python\n')
for run in runs:
  cmd = '{}/run/combine_datasets.exe -i {} -o {} -f {} -b {} -e {} -y Run{} -n {}'.format(os.getcwd(), in_dir, out_dir, dataset_list, run[0], run[-1], year, nano_version)
  cmdfile.write('print(\"'+cmd+'\")\n')
cmdfile.close()
os.chmod(args.out_cmd_file, 0o755)

print("To generate job json and submit jobs do: ")
print('convert_cl_to_jobs_info.py '+args.out_cmd_file+' '+'<filename>.json')
print('auto_submit_jobs.py <filename>.json')

sys.exit(0)
