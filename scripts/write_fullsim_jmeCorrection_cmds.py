#!/usr/bin/env python
# Writes commands to apply JECs on nanoAODs using apply_fastsim_jmeCorrections_to_nanoaod.py
# Example: ./scripts/write_fastsim_jmeCorrection_cmds.py 
#  -i /net/cms25/cms25r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D_unsplit 
#  -o /net/cms25/cms25r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D_unsplit_fastSimJmeCorrection 
#  -c apply_fastsim_jmeCorrection_2016.py

import os, sys
import argparse
import glob
parser = argparse.ArgumentParser()
parser.add_argument('-i','--inputNanoAodFolder', required=True, default="/net/cms25/cms25r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D_unsplit")
parser.add_argument('-o','--outputNanoAodFolder', required=True, default="/net/cms25/cms25r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D_unsplit_fastSimJmeCorrection")
parser.add_argument('-c', '--commandFilename', required=True, default="apply_fastsim_jmeCorrection_2016.py")
parser.add_argument('-g', '--glob', default='*.root')
args = parser.parse_args()

year = ""
if args.inputNanoAodFolder.find("/2016/")!= -1: year = "2016"
elif args.inputNanoAodFolder.find("/2017/")!= -1: year = "2017"
elif args.inputNanoAodFolder.find("/2018/")!= -1: year = "2018"

filenames = glob.glob(args.inputNanoAodFolder+"/"+args.glob)
commandFile = open(args.commandFilename, 'w')
commandFile.write("#!/bin/env python\n")
for filename in filenames:
  command = "print ('./scripts/apply_fullsim_jmeCorrections_to_nanoaod.py -o "+args.outputNanoAodFolder+" -y "+year+" -i "+filename+"')\n"
  commandFile.write(command)
commandFile.close()

if not os.path.exists(args.outputNanoAodFolder):
  os.makedirs(args.outputNanoAodFolder)

os.chmod(args.commandFilename, 0o755)
os.system('convert_cl_to_jobs_info.py '+args.commandFilename+' '+ args.commandFilename.replace('.py','.json'))

print('Last line in '+args.commandFilename+' is:')
os.system('cat '+args.commandFilename+' | tail -n 1\n')

print("To submit jobs do: ")
print('auto_submit_jobs.py '+args.commandFilename.replace('.py','.json')+' -c scripts/check_jmeCorrection.py')
