#!/usr/bin/env python3

import os, sys, argparse
from utilities import bcolors
from glob import glob

parser = argparse.ArgumentParser(description="Update pico tree with DNN output.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-p","--pico_dir", 
                    default="/net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/SMS-TChiHH_2D/mergednn_higmc_higloose/",
                    help="Directory where the input files are.")
parser.add_argument("-d","--dnnout_dir", 
                    default="/net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/SMS-TChiHH_2D/dnnout_higloose/",
                    help="Directory where the output files should go.")
args = vars(parser.parse_args())

pico_dir = args['pico_dir']
if (pico_dir[-1]!='/'): pico_dir = pico_dir + '/'
pico_paths = glob(os.path.join(pico_dir,'*.root'))
print('Found {} input files.\n'.format(len(pico_paths)))

for ipico_path in pico_paths:
  ipico_file = os.path.basename(ipico_path)
  idnnout_path = os.path.join(args['dnnout_dir'], 'dnnout_'+ipico_file)
  if os.path.exists(idnnout_path):
    cmd = './run/update_pico.exe -p {} -d {}'.format(ipico_path, idnnout_path)
    print(cmd)
    os.system(cmd)
  else:
    print('Could not find dnn output file:',idnnout_path)


print(bcolors.OKGREEN+'\nYou are finally done with processing and it\'s time for some fun!!!\n'+bcolors.ENDC)

