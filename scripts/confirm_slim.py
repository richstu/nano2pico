#!/usr/bin/env python
import ROOT
import sys

unslimmed_dir = sys.argv[1]
slimmed_dir = sys.argv[2]

unslimmed = ROOT.TChain("tree")
unslimmed.Add(unslimmed_dir+'/*.root')
unslimmed_events = unslimmed.GetEntries()

slimmed = ROOT.TChain("tree")
slimmed.Add(slimmed_dir+'/*.root')
slimmed_events = slimmed.GetEntries()

print(unslimmed_dir+' has '+str(unslimmed_events)+' events')
print(slimmed_dir+' has '+str(slimmed_events)+' events')

if unslimmed_events != slimmed_events:
  print('[Error] unslimmed and slimmed do not have same number of events')
  sys.exit(1)
else:
  print('[Success] unslimmed and slimmed have same number of events')
