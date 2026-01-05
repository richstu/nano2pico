#!/usr/bin/env python

import subprocess
from argparse import ArgumentParser

def print_and_run(command):
  '''
  Given string command, print and run it through bash
  '''
  print(command)
  subprocess.call(command.split())

if __name__=='__main__':
  #parse arguments
  argument_parser = ArgumentParser(prog='test_process_nano',
      description='Runs nano2pico on a test file')
  argument_parser.add_argument('-i','--input_filename')
  argument_parser.add_argument('-n','--nevents')
  argument_parser.add_argument('-p','--process_nano_only',action='store_true')
  argument_parser.add_argument('-c','--clean',action='store_true')
  args = argument_parser.parse_args()

  if args.input_filename != None:
    input_file = args.input_filename
    path_pos = input_file.rfind('/')
    indir = input_file[:path_pos]
    infile_name = input_file[path_pos+1:]

    #process nano
    cmd = './run/process_nano.exe --in_file '+infile_name+' --in_dir '+indir+' --out_dir out/zgamma/'
    if args.nevents != None:
      cmd += ' --nent '+args.nevents
    print_and_run(cmd)

    if not args.process_nano_only:
      #merge corrections
      cmd = './run/merge_corrections.exe out/zgamma/corrections/corr_'+infile_name
      cmd += ' out/zgamma/wgt_sums/wgt_sums_'+infile_name
      print_and_run(cmd)

      #apply corrections
      cmd = './run/apply_corrections.exe --in_file raw_pico_'+infile_name+' --in_dir out/zgamma/raw_pico/ --corr_file corr_'+infile_name
      print_and_run(cmd)

  if args.clean:
    subprocess.call('rm out/zgamma/unskimmed/*.root',shell=True)
    subprocess.call('rm out/zgamma/corrections/*.root',shell=True)
    subprocess.call('rm out/zgamma/wgt_sums/*.root',shell=True)
    subprocess.call('rm out/zgamma/raw_pico/*.root',shell=True)
