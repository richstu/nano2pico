#!/usr/bin/env python
#script used to quickly run tests after nano2pico changes

from argparse import ArgumentParser

import subprocess

def print_and_run(command):
  '''
  Given string command, print and run it through bash
  '''
  print(command)
  subprocess.call(command.split())

if __name__=='__main__':
  argument_parser = ArgumentParser(prog='processnano')
  argument_parser.add_argument('-i','--input_file')
  argument_parser.add_argument('-n','--nent')
  args = argument_parser.parse_args()

  path_pos = args.input_file.rfind('/')
  indir = args.input_file[:path_pos]
  infile = args.input_file[path_pos+1:]
  cmd = './run/process_nano.exe --in_file '+infile+' --in_dir '+indir+' --out_dir out/zgamma/'
  if args.nent != None:
    cmd += ' --nent {}'.format(args.nent)
  print_and_run(cmd)

