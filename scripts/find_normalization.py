#! /usr/bin/env python
"""@package docstring 
Script that aggregates the sum-of-weights for NanoAOD files and stores them in 
a json file
"""

import argparse
from glob import glob
import json
import os
import re
import ROOT

def full_path(path: str) -> str:
  """Expands path including symlinks

  Args:
    path: path to expand, possibly with symlinks

  Returns:
    expanded path
  """
  return os.path.realpath(os.path.abspath(os.path.expanduser(path)))

def get_tag(path: str) -> str:
  """Gets name ("tag") of dataset from NanoAOD by ignoring extensions and hash

  Args:
    path: NanoAOD filename

  Returns:
    NanoAOD dataset name without directory, extension, or hash
  """
  tag = path.split("/")[-1]
  tag = re.findall('(.*?)__',tag)[0]
  return tag

def find_normalization(nano_dir: str, output_name: str, overwrite: bool, 
                       dataset_list: str, verbose: bool):
  """Aggregates sum-of-weights for NanoAOD files in nano_dir and writes them to
  txt/wgt_sums/<output_name>

  Args:
    nano_dir: Directory to parse
    output_name: output json filename
    overwrite: if false, will throw an error if nano_dir has been parsed before
               in output file
    dataset_list: file containing list of data sets to process
    verbose: print debug info
  """
  if dataset_list != "":
    if not os.path.isfile(dataset_list):
      raise RuntimeError("Data set list {dataset_list} does not exist.")

  # load data from file if it already exists and check for overwrite
  wgt_sum_dir = os.path.join(os.path.join(os.getcwd(),"txt"),"wgt_sums")
  if not os.path.isdir(wgt_sum_dir):
    os.mkdir(wgt_sum_dir)

  output_filename = os.path.join(wgt_sum_dir,output_name)
  wgt_sums = {}
  if os.path.isfile(output_filename):
    with open(output_filename,"r") as output_file:
      wgt_sums = json.loads(output_file.read())

  nano_dir = full_path(nano_dir)
  if (nano_dir in wgt_sums) and (not overwrite):
    raise RuntimeError("Directory wgt_sums exists in output file. "
                       "Use --overwrite to overwrite existing content.")
  if not (nano_dir in wgt_sums):
    wgt_sums[nano_dir] = {}

  # get tags for all NanoAOD data sets 
  if dataset_list != "":
    tags = set()
    tag_filelist = dict()
    with open(dataset_list, "r") as dataset_file:
      lines = dataset_file.read().split("\n")
      for line in lines:
        if len(line) > 0:
          dataset_info = line.split("/")
          tag = dataset_info[1]
          wildcard = (dataset_info[1]+"__"+dataset_info[2].replace("-","__",1)
                      +"*.root")
          if tag in tags:
            tag_filelist[tag] += glob(os.path.join(nano_dir, wildcard))
          else:
            tags.add(tag)
            tag_filelist[tag] = glob(os.path.join(nano_dir, wildcard))
  else:
    input_files = [os.path.join(nano_dir,f) for f in os.listdir(nano_dir)
                   if os.path.isfile(os.path.join(nano_dir, f))
                   and os.path.splitext(f)[1] == ".root"]
    tags = list(set([get_tag(f) for f in input_files]))
  if verbose:
    print(f"Found {len(tags)} tags")

  # sum weights for each tag
  for tag in tags:
    if verbose:
      print(f"Processing tag {tag}")
    wgt_sums[nano_dir][tag] = {}
    if dataset_list != "":
      tag_files = tag_filelist[tag]
    else:
      tag_files = glob(os.path.join(nano_dir, tag)+"*.root")
    gen_event_sumw = 0.0
    lhe_scale_sumw = [0.0 for imurf in range(9)]
    for tag_file in tag_files:
      with ROOT.TFile(tag_file, "READ") as nano_file:
        for event in nano_file.Runs: #only 1 event in metadata tree
          try:
            gen_event_sumw += event.genEventSumw
          except:
            print(f"WARNING: genEventSumw not set in {tag_file}")
          try:
            if len(event.LHEScaleSumw) == 9:
              for imurf in range(9):
                lhe_scale_sumw[imurf] += event.LHEScaleSumw[imurf]
          except:
            print(f"WARNING: LHEScaleSumw not set in {tag_file}")
    wgt_sums[nano_dir][tag]['genEventSumw'] = gen_event_sumw
    for imurf in range(9):
      wgt_sums[nano_dir][tag][f'LHEScaleSumw{imurf}'] = lhe_scale_sumw[imurf]
  
  # save output
  with open(output_filename,"w") as output_file:
    output_file.write(json.dumps(wgt_sums, indent=2))

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Script that aggregates the "
      "sum-of-weights for NanoAOD files and stores them in a json file",
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("-n", "--nano_dir")
  parser.add_argument("-d", "--dataset_list", default="")
  parser.add_argument("-o", "--output_name", default="wgt_sums.json")
  parser.add_argument("-w", "--overwrite", action="store_true")
  parser.add_argument("-v", "--verbose", action="store_true")
  args = parser.parse_args()

  find_normalization(args.nano_dir, args.output_name, args.overwrite, 
                     args.dataset_list, args.verbose)
