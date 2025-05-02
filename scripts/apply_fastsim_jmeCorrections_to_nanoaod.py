#!/usr/bin/env python3
# Applies fastsim jmeCorrections to NanoAOD using CMSSW
# Example: python apply_fastsim_jmeCorrections_to_nanoaod.py -o FOLDER -y YEAR -i FILENAME...
import os, sys
import argparse
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *

parser = argparse.ArgumentParser()
parser.add_argument('-o','--output', required=True)
parser.add_argument('-y','--year', required=True, default="2016")
parser.add_argument('-i','--input', required=True, nargs='+')
parser.add_argument('-c','--cut_string', default=None)
parser.add_argument('-m','--max_events', default=None)
args = parser.parse_args()

metBranchName = "MET" if args.year != "2017" else "METFixEE2017"
#metBranchName = "MET"

jmeCorrections = createJMECorrector(isMC=True, dataYear=args.year, runPeriod="A", jesUncert="Total", jetType="AK4PFchs", noGroom=False, metBranchName=metBranchName, applySmearing=True, isFastSim=True)
print(args.max_events)
p=PostProcessor(outputDir=args.output,inputFiles=args.input,cut=args.cut_string,modules=[jmeCorrections()],provenance=True, postfix="", maxEntries=args.max_events)
p.run()
