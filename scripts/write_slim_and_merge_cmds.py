#!/usr/bin/env python3

import os, argparse
from glob import glob

def findBaseSampleNames(folder):
  infiles = set() # to remove duplicates
  for file in glob(folder+'/*.root'):
    dataset_tag = file.split('/')[-1]
    dataset_tag = dataset_tag.split('__RunIISummer16NanoAODv5__')[0]
    dataset_tag = dataset_tag.split('_ext')[0]
    infiles.add(dataset_tag)
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
  parser.add_argument('--overwrite', action='store_true',
                    help='Process all input files regardless whether output exists.')
  parser.add_argument('-t', '--tag', default='',
                    help='Optionally specify a tag to be used to differentiate helper files for batch submission.')
  args = vars(parser.parse_args())

  in_dir = args['in_dir']
  skim_name = in_dir.rstrip('/').split('/')[-1].replace('skim_','')

  slim_name = args['slim_name']
  slim_rules = os.path.join(os.getcwd(),'txt/slim_rules/'+slim_name+'.txt')

  enclosing_dir = os.path.dirname(os.path.dirname(in_dir))
  out_dir = os.path.join(enclosing_dir, 'merged_'+slim_name+'_'+skim_name)
  if not os.path.exists(out_dir): 
    os.mkdir(out_dir)

  dataset_tags = findBaseSampleNames(in_dir)
  print('Found {} dataset_tags.\n'.format(len(dataset_tags)))

  cmdfile_name = 'slim_'+slim_name+'_'+skim_name+'_cmds.py'
  if (args['tag']!=''): cmdfile_name = args['tag']+'_'+cmdfile_name
  cmdfile = open(cmdfile_name,'w')
  cmdfile.write('#!/bin/env python\n')
  nexisting = 0
  for dstag in dataset_tags:
    #if 'TTJets_SingleLeptFromT_' not in dstag: continue
    in_file_paths = os.path.join(in_dir,'*'+dstag+'*.root')
    out_file_name = 'merged_'+dstag+'_'+slim_name+'_'+skim_name+'_nfiles_'+str(len(glob(in_file_paths)))
    out_file_path = os.path.join(out_dir,out_file_name+'.root')

    if os.path.exists(out_file_path):
      nexisting +=1
      if not args['overwrite']: 
        continue
    cmd = '{}/scripts/slim_and_merge.py -s {} -o {} -i {}'.format(os.getcwd(), slim_rules, out_file_path, in_file_paths)
    cmdfile.write('print(\''+cmd+'\')\n')

  cmdfile.close()
  os.chmod(cmdfile_name, 0o755)

  if not args['overwrite']:
    print('Found existing {} output files, will submit jobs only for remainder. Use --overwrite to ignore the existing output and run all jobs.\n'.format(nexisting))
  else:
    print('Found existing {} output files, which will be overwritten.\n'.format(nexisting))

  json_name = cmdfile_name.replace('.py','.json')
  print('To print a sample command:')
  print('cat '+cmdfile_name+' | tail -n 1\n')
  os.system('convert_cl_to_jobs_info.py '+cmdfile_name+' '+json_name)
  print('\nTo generate job json and submit jobs:')
  print('auto_submit_jobs.py '+json_name)



