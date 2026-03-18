"""Python script that produces NanoAOD dataset names from (pico) files
01/26 M Oshiro
"""

import re
from glob import glob
from argparse import ArgumentParser

if __name__=='__main__':
  parser = ArgumentParser()
  parser.add_argument('-i','--input_dir', default='')
  parser.add_argument('-o','--output', default='dataset_list')
  args = parser.parse_args()
  filename_list = glob(args.input_dir+'/*.root')
  dataset_names = []
  for filename in filename_list:
    #transform into dataset name
    filename = filename[len(args.input_dir):]
    filename = filename.replace('pico_','/').replace('pythia8__','pythia8/')
    match = re.search('v\d+(_ext\d+)?-v\d+',filename)
    if not match:
      raise RuntimeError(
          f'Could not find sample version in filename {filename}')
    filename = filename[:match.end()]+'/NANOAODSIM'
    match = re.search('__(1\d\dX)',filename)
    if not match:
      match = re.search('__(forPOG_1\d\dX)',filename)
      if not match:
        match = re.search('__(without_JHUGEN_1\d\dX)',filename)
        if not match:
          raise RuntimeError(
              f'Could not find CMSSW version in filename {filename}')
    filename = filename[:match.start()]+'-'+match[1]+filename[match.end():]
    if not (filename in dataset_names):
      dataset_names.append(filename)
  dataset_names.sort()
  with open(args.output, 'w') as output_file:
    first = True
    for dataset_name in dataset_names:
      if first:
        first = False
      else:
        output_file.write('\n')
      output_file.write(dataset_name)
