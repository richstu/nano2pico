#!/usr/bin/env python3
import argparse

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='''Organize cross section log''', formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('-l','--log', required=True, help='Cross section log filename.')
  args = parser.parse_args()

  #cross_section_filename = 'unit_test_cross_section.log'
  cross_section_filename = args.log

  missing_datasets = set()

  with open(cross_section_filename) as cross_section_file:
    for line in cross_section_file:
      if '-999999.0 pb' in line:
        dataset = line.strip().split('/')[-1].split('__')[0]
        missing_datasets.add(dataset)

  for dataset in sorted(missing_datasets):
    print(dataset)
