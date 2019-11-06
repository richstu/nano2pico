#!/usr/bin/env python3

import os, argparse
from glob import glob

def findBaseSampleNames(folder):
  infiles = set() # to remove duplicates
  for file in glob(folder+'/*.root'):
    tag = file.split('/')[-1]
    tag = tag.split('__RunIISummer16NanoAODv5__')[0]
    infiles.add(tag)
    sortedfiles = list()
  for file in infiles:
    sortedfiles.append(file)
  sortedfiles = sorted(sortedfiles)

  return sortedfiles

if __name__ == '__main__':

  parser = argparse.ArgumentParser(description='Submits batch jobs to apply new SFs and compute sum-of-weights',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('-i','--in_dir', required=True, 
                      default='/net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/skim_met150/',
                      help='Directory where the skim is.')
  parser.add_argument('-l','--slim_name', required=True, default='higmc',
                      help='Name of the slim, e.g. "higmc". N.B. code then assumes that txt/slim_rules/higmc.txt exists.')
  parser.add_argument('--overwrite', default=False,
                    help='Process all input files regardless whether output exists.')
  args = vars(parser.parse_args())

  in_dir = args['in_dir']
  skim_name = in_dir.rstrip('/').split('/')[-1].replace('skim_','')

  slim_name = args['slim_name']
  slim_rules = os.path.join(os.getcwd(),'txt/slim_rules/'+slim_name+'.txt')

  enclosing_dir = os.path.dirname(os.path.dirname(in_dir))
  out_dir = os.path.join(enclosing_dir, 'merged_'+slim_name+'_'+skim_name)
  if not os.path.exists(out_dir): 
    os.mkdir(out_dir)

  tags = findBaseSampleNames(in_dir)
  print('Found {} tags.\n'.format(len(tags)))

  cmdfile_name = 'slim_'+slim_name+'_'+skim_name+'_cmds.py'
  cmdfile = open(cmdfile_name,'w')
  cmdfile.write('#!/bin/env python\n')
  for tag in tags:
    #if 'TTJets_SingleLeptFromT_' not in tag: continue
    in_file_paths = os.path.join(in_dir,'*'+tag+'*.root')
    out_file_name = 'merged_'+tag+'_'+slim_name+'_'+skim_name+'_nfiles_'+str(len(glob(in_file_paths)))
    out_file_path = os.path.join(out_dir,out_file_name+'.root')

    if not args['overwrite'] and os.path.exists(out_file_path):
      continue

    cmd = '{}/scripts/slim_and_merge.py -s {} -o {} -i {}'.format(os.getcwd(), slim_rules, out_file_path, in_file_paths)
    cmdfile.write('print(\''+cmd+'\')\n')

  cmdfile.close()
  os.chmod(cmdfile_name, 0o755)

  json_name = 'slim_'+slim_name+'_'+skim_name+'_json.py'
  print('To print a sample command:')
  print('cat '+cmdfile_name+' | tail -n 1\n')
  os.system('convert_cl_to_jobs_info.py '+cmdfile_name+' '+json_name)
  print('\nTo generate job json and submit jobs:')
  print('auto_submit_jobs.py '+json_name)



