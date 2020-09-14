#!/usr/bin/env python3
import subprocess
import argparse
import glob
import re

# The inputs needed to run this are at: /afs/cern.ch/user/a/amete/public/EWKGauginoCrossSections_13TeV
# N.B. Initiliaze the vectors connected to the branches to 0 to avoid seg faults in get_gaugino.C
# then run this inside the folder to avoid the need for further modifications

if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser.add_argument('-i','--signal_path', required=True, default = "")
  parser.add_argument('-m','--model', default = "CN")
  args = parser.parse_args()

  doRound =True

  # Find mass of higgsino
  signal_files = glob.glob(args.signal_path+'/*.root')
  chi_mass_list = set()
  for signal_file in signal_files:
    mChi = re.findall(r"mChi-\d+",signal_file)[0].replace('mChi-','')
    mChi = int(mChi)
    if (doRound):
      if (mChi%25!=0 and mChi!=127): 
        mChi = int(round(mChi/25))*25
    chi_mass_list.add(mChi)
  print("Mass points")
  print(sorted(chi_mass_list))

  # Print cross sections
  for mass in sorted(chi_mass_list):
    if mass == 1500: continue
    result = subprocess.check_output('root -l -q \'get_gaugino.C("'+args.model+'","hino",'+str(mass)+')\'', shell=True)
    # print(result)
    result = result.decode()
    xsec = float(result.split(' is ')[-1].split(' [pb] ')[0])
    xsec_unc = float(result.split(' +/- ')[-1].split(' [rel')[0])
    if (doRound):
      print('else if(hig_mass >{}-12 && hig_mass <= {}+12) {{ xsec = .5824*.5824*{}; xsec_unc = {}; return;}}'.format(mass, mass, xsec, xsec_unc))
    else:
      print('else if(hig_mass =={}) {{ xsec = .5824*.5824*{}; xsec_unc = {}; return;}}'.format(mass, xsec, xsec_unc))

  #for mass in range(125,1501,25):
  #  if mass==125: mass=127
  #  result = subprocess.check_output('root -l -q \'get_gaugino.C("N1N2","hino",'+str(mass)+')\'', shell=True)
  #  # print(result)
  #  result = result.decode()
  #  xsec = float(result.split(' is ')[-1].split(' [pb] ')[0])
  #  xsec_unc = float(result.split(' +/- ')[-1].split(' [rel')[0])
  #  print('else if(hig_mass =={}) {{ xsec = .5824*.5824*{}; xsec_unc = {}; return;}}'.format(mass, xsec, xsec_unc))
