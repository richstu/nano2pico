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
  parser.add_argument("--overwrite", action='store_true',
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
  out_dir = os.path.join(enclosing_dir,'skim_'+skim_name+'/')
  # in case running on top of a slim
  if '/merged_' in in_dir: 
    slim_folder_split = in_dir.rstrip('/').split('/')[-1].split('_')
    slim_name = slim_folder_split[0]+'_'+slim_folder_split[1]+'_'+skim_name # keep the slim rules name
    out_dir = out_dir.replace('skim_'+skim_name, slim_name)
  if not os.path.exists(out_dir): 
    os.mkdir(out_dir)

  cmdfile_name = 'skim_'+skim_name+'_cmds.py'
  if (args['tag']!=''): cmdfile_name = cmdfile_name.replace('.py', '_'+args['tag']+'.py')
  cmdfile = open(cmdfile_name,'w')
  cmdfile.write('#!/bin/env python\n')
  nexisting=0
  for ifile_path in in_file_paths:
    if "raw_pico" in ifile_path:
      out_file_path = ifile_path.replace(in_dir,out_dir).replace('raw_pico_','pico_'+skim_name+'_')
    else: 
      out_file_path = ifile_path.replace(in_dir,out_dir).replace('pico_','pico_'+skim_name+'_')
    if os.path.exists(out_file_path):
      nexisting +=1
      if not args['overwrite']: 
        continue
    cmd = '{}/scripts/skim_file.py -k {} -i {} -o {}'.format(os.getcwd(), skim_name, ifile_path, out_dir)
    cmdfile.write('print(\''+cmd+'\')\n')

  cmdfile.close()
  os.chmod(cmdfile_name, 0o755)
  if not args['overwrite']:
    print('Found existing {} output files, so will submit  {} jobs. Use --overwrite to ignore the existing output and run all jobs.\n'.format(nexisting, len(in_file_paths)-nexisting))
  else:
    print('Found existing {} output files, which will be overwritten.\n'.format(nexisting))

  json_name = cmdfile_name.replace('.py','.json')
  os.system('convert_cl_to_jobs_info.py -f '+cmdfile_name+' '+json_name)

  print('To print a sample command:')
  print('cat '+cmdfile_name+' | tail -n 1\n')

  print('To generate job json and submit jobs:')
  print('auto_submit_jobs.py '+json_name+' -c scripts/check_skim.py')



