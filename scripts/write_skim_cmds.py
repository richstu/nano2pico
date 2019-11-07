#!/usr/bin/env python3

import os, argparse
from glob import glob

if __name__ == '__main__':

  parser = argparse.ArgumentParser(description='Submits batch jobs to apply new SFs and compute sum-of-weights',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('-i','--in_dir', required=True, 
                      default='/net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/unskimmed/',
                      help='Directory where the raw pico files are')
  parser.add_argument('-k','--skim_name', required=True, default='',
                      help='Plain text name for the skim. Output folder name will be named according to this.')
  parser.add_argument("--overwrite", default=False,
                    help="Process all input files regardless whether output exists.")
  parser.add_argument('-t', '--tag', default='',
                    help='Optionally specify a tag to be used to differentiate helper files for batch submission.')
  args = vars(parser.parse_args())

  skim_name = args['skim_name']
  in_dir = args['in_dir']
  if (in_dir[-1]!='/'): in_dir = in_dir + '/'
  in_file_paths = glob(os.path.join(in_dir,'*.root'))
  print('Found {} input files.\n'.format(len(in_file_paths)))

  enclosing_dir = os.path.dirname(os.path.dirname(in_dir))
  out_dir = os.path.join(enclosing_dir,'skim_'+skim_name)
  # in case running on top of a slim
  if '/merged_' in in_dir: 
    slim_folder_split = in_dir.rstrip('/').split('/')[-1].split('_')
    slim_name = slim_folder_split[0]+'_'+slim_folder_split[1]+'_'+skim_name # keep the slim rules name
    out_dir = out_dir.replace('skim_'+skim_name, slim_name)
  if not os.path.exists(out_dir): 
    os.mkdir(out_dir)

  cmdfile_name = 'skim_'+skim_name+'_cmds.py'
  if (args['tag']!=''): cmdfile_name = args['tag']+'_'+cmdfile_name
  cmdfile = open(cmdfile_name,'w')
  cmdfile.write('#!/bin/env python\n')
  for ifile_path in in_file_paths:
    out_file_path = ifile_path.replace(in_dir,out_dir).replace('pico_','pico_'+skim_name+'_')
    if not args['overwrite'] and os.path.exists(out_file_path):
      continue
    cmd = '{}/scripts/skim_file.py -k {} -i {} -o {}'.format(os.getcwd(), skim_name, ifile_path, out_dir)
    cmdfile.write('print(\''+cmd+'\')\n')

  cmdfile.close()
  os.chmod(cmdfile_name, 0o755)

  json_name = cmdfile_name.replace('.py','.json')
  print('To print a sample command:')
  print('cat '+cmdfile_name+' | tail -n 1\n')
  print('To generate job json and submit jobs:')
  print('convert_cl_to_jobs_info.py '+cmdfile_name+' '+json_name)
  print('auto_submit_jobs.py '+json_name)



