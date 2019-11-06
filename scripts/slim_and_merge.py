#! /usr/bin/env python

from __future__ import print_function

import sys
import argparse
import fnmatch
import os

from ROOT import TChain, TFile

def getRules(slim_rules_file_name):
    rules = [ line.strip().split() for line in open(slim_rules_file_name) ]
    good_rules = [ rule for rule in rules
                   if len(rule)==0
                   or (len(rule)>0 and rule[0].startswith("#"))
                   or (len(rule)>=2 and (rule[0]=="keep" or rule[0]=="drop")) ]
    bad_rules = [ rule for rule in rules if rule not in good_rules ]
    good_rules = [ rule for rule in good_rules if len(rule)>=2 ]
    for rule in bad_rules:
        utilities.ePrint("Invalid rule:",rule,"\n")
    return good_rules

def passRules(branch, rules):
    matched_rules =  [ rule for rule in rules if fnmatch.fnmatch(branch, rule[1]) ]
    return len(matched_rules)==0 or matched_rules[-1][0] == "keep"

def sortInputFilesBySize(input_file_paths):
    input_file_paths = [ (f, os.path.getsize(f)) for f in input_file_paths ]

    input_file_paths.sort(key=lambda f: f[1], reverse=True)

    input_file_paths = [ f[0] for f in input_file_paths ]

    return input_file_paths

def slimFile(slim_rules_file_name, out_file_path, input_file_paths, test_mode):
    print("     INPUT FILES:",input_file_paths,"\n")
    print("     OUTPUT FILE:",out_file_path,"\n")
    print("      RULES FILE:",slim_rules_file_name,"\n")

    in_tree = TChain("tree", "tree")

    input_file_paths = sortInputFilesBySize(input_file_paths)

    for input_file_name in input_file_paths:
        in_tree.Add(input_file_name)

    branch_names = [ branch.GetName() for branch in in_tree.GetListOfBranches() ]
    rules = getRules(slim_rules_file_name)
    kept_branches = [ branch for branch in branch_names if passRules(branch, rules) ]
    kept_branches.sort()
    dropped_branches = [ branch for branch in branch_names if branch not in kept_branches ]
    dropped_branches.sort()

    print("DROPPED BRANCHES:",dropped_branches,"\n")
    print("   KEPT BRANCHES:",kept_branches,"\n")
    if test_mode: return

    for branch in branch_names:
        if branch in kept_branches: in_tree.SetBranchStatus(branch, True)
        else:                       in_tree.SetBranchStatus(branch, False)

    out_file = TFile(out_file_path, "recreate")
    in_tree.Merge(out_file, 0, "fast keep")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prunes branches from an ntuple",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-t", "--test", action="store_true",
                        help="Run in test mode, quickly diplaying the list of kept and dropped branchs without actually copying the trees.")
    parser.add_argument("-s", "--slim_rules_file",
                        help="File containing rules for pruning branches (one rule per line). Rules are are the form \"keep XXX\" or \"drop YYY\". Unix shell-style wildcards (e.g., '*') allow pattern matching. Branches are kept by default if no matching rule is found for the branch. If multiple rules match, the last takes precedence.")
    parser.add_argument("-o", "--out_file_path",
                        help="File in which to save the slimmed and merged ntuple.")
    parser.add_argument("-i", "--input_file_paths", nargs="+",
                        help="Files containing ntuples to be slimmed and merged.")
    args = parser.parse_args()

    slimFile(args.slim_rules_file, args.out_file_path, args.input_file_paths, args.test)
