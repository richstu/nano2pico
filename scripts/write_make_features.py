#!/usr/bin/env python3

import os, sys, argparse
from utilities import bcolors
from glob import glob

parser = argparse.ArgumentParser(description="Submits batch jobs to make feature ntuples for NN training.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i","--in_dir", default="",
                    help="Directory where the pico files are")
parser.add_argument('-t', '--tag', default='',
                  help='Optionally specify a tag to be used to differentiate helper files for batch submission.')
args = vars(parser.parse_args())

in_dir = args['in_dir']
if (in_dir[-1]!='/'): in_dir = in_dir + '/'
in_file_paths = glob(os.path.join(in_dir,'*.root'))
print('Found {} input files.\n'.format(len(in_file_paths)))
feats = 'zg'
if 'higgsino' in in_dir:
  feats = 'hig'

out_dir = in_dir.replace(in_dir.split('/')[-2],'features_'+in_dir.split('/')[-2])
if not os.path.exists(out_dir): 
  os.makedirs(out_dir)

cmdfile_name = 'feature_cmds.py'
if (args['tag']!=''): cmdfile_name = args['tag']+'_'+cmdfile_name
cmdfile = open(cmdfile_name,'w')
cmdfile.write('#!/bin/env python\n')
cmd =''
for ifile_path in in_file_paths:
  ifile = os.path.basename(os.path.realpath(ifile_path))
  cmd = '{}/run/make_{}feats.exe -f {} -i {} -o {}'.format(os.getcwd(), feats, ifile, in_dir, out_dir)
  cmdfile.write('print(\"'+cmd+'\")\n')

cmdfile.close()
os.chmod(cmdfile_name, 0o755)

json_name = cmdfile_name.replace('.py','.json')
os.system('convert_cl_to_jobs_info.py '+cmdfile_name+' '+json_name)

print('Example command:')
print(bcolors.WARNING+cmd+bcolors.ENDC)

print("To submit jobs do: ")
print('auto_submit_jobs.py '+json_name)
  
