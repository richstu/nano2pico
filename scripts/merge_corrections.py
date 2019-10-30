#! /usr/bin/env python

from __future__ import print_function

import argparse
import subprocess
import os

def fullPath(path):
  return os.path.realpath(os.path.abspath(os.path.expanduser(path)))

def ensureDir(path):
  try:
    os.makedirs(path)
  except OSError:
    if not os.path.isdir(path):
      raise

def getTag(path):
  tag = path.split("/")[-1]
  tag = tag.split("RunIISummer16NanoAODv5")[0]
  tag = tag.replace("wgt_sums_","")
  tag = tag.strip("_")
  return tag

def mergeCorrections(wgt_dir, corr_dir, year):
  wgt_dir = fullPath(wgt_dir)
  corr_dir = fullPath(corr_dir)

  ensureDir(corr_dir)
  
  input_files = [os.path.join(wgt_dir,f) for f in os.listdir(wgt_dir)
                 if os.path.isfile(os.path.join(wgt_dir, f))
                 and os.path.splitext(f)[1] == ".root"]

  tags = list(set([getTag(f) for f in input_files]))

  for i in range(len(tags)):
    tag = tags[i]
    output_file = os.path.join(corr_dir, "corr_"+tag+".root")
    if os.path.exists(output_file):
      print("Processing tag {} of {}: Output file already exists. Continue.".format(i+1,len(tags)))
      continue
    command = ["run/merge_corrections.exe", str(year), output_file]
    for f in input_files:
      if tag in f:
        command.append(f)
    print("Processing tag {} of {}: {}".format(i+1,len(tags),tag))
    subprocess.call(command)

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Merges multiple sum-of-weights files into one corrections file per tag.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("-w", "--wgt_dir", default=os.getcwd()+"/wgt_sums/",
                      help="Directory from which to read sum-of-weights files")
  parser.add_argument("-c", "--corr_dir", default=os.getcwd()+"/corrections/",
                      help="Directory in which to store corrections files")
  parser.add_argument("-y","--year", type=int, default=2016, help="Sample year.")
  args = parser.parse_args()

  mergeCorrections(args.wgt_dir, args.corr_dir, args.year)
