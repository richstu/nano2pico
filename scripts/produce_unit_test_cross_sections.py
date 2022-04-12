#!/bin/sh
''''exec python3 -u -- "$0" ${1+"$@"} # '''
import os
import re
import subprocess
import threading
import time
import multiprocessing
import argparse

def output_reader(process, commandOutput):
  for line in iter(process.stdout.readline, b''):
    output = line.decode('utf-8')
    print(output, end='')
    commandOutput[0] += output

def runCommand(command):
  print("\n[Info] Running command: "+command)
  process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)

  commandOutput = ['']
  thread = threading.Thread(target=output_reader, args=(process,commandOutput))
  thread.start()
  thread.join()
  # Try to get poll
  for iTime in range(60):
    if (process.poll() == None): time.sleep(1)
    else: break
  return process.poll(), commandOutput[0]

def get_nanoaod_files(folders):
  nanoaod_files = []
  for folder in folders:
    for obj in os.listdir(folder):
      if ".root" not in obj: continue
      path = os.path.join(folder,obj)
      year = re.findall("nano/(\d\d\d\d)/", folder)[0]
      nanoaod_files.append([path, year])
  return nanoaod_files
  
def get_cross_sections(files):
  # cross_sections = [(filename, cross_section in pb)]
  cross_sections = []
  # Make commands
  command_list = []
  for file_info in files:
    file_path = file_info[0]
    year = file_info[1]
    command = './run/print_cross_sections.exe -f '+file_path+' -y '+year
    command_list.append(command)
  # Run commands
  pool = multiprocessing.Pool()
  command_results = pool.map(runCommand, command_list)
  # Organize results
  for command_result in command_results:
    return_code, output = command_result
    cross_section = float(re.findall('is ([-]?\d.*)pb', output)[0])
    cross_sections.append([file_path, cross_section])
  return cross_sections

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='''Gets cross sections that will be used for NanoAOD files.''', formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('-l','--output_log', required=True, help='Ouptut log filename.')
  args = parser.parse_args()

  #cross_section_log_name = "unit_test_cross_section.log"
  cross_section_log_name = args.output_log

  # Print cross sections for NanoAODv7 mc files
  # Get all NanoAODv7 files
  nanoaodv7_folders = ["/net/cms17/cms17r0/pico/NanoAODv7/nano/2016/mc", 
                       "/net/cms24/cms24r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D_fastSimJmeCorrection",
                       "/net/cms17/cms17r0/pico/NanoAODv7/nano/2017/mc",
                       "/net/cms24/cms24r0/pico/NanoAODv7/nano/2017/SMS-TChiHH_2D_fastSimJmeCorrection",
                       "/net/cms17/cms17r0/pico/NanoAODv7/nano/2018/mc",
                       "/net/cms24/cms24r0/pico/NanoAODv7/nano/2018/SMS-TChiHH_2D_fastSimJmeCorrection",
                       ]

  #nanoaodv7_files = [(path, year)]
  nanoaodv7_files = get_nanoaod_files(nanoaodv7_folders)
  #nanoaodv7_cross_sections = [(path, cross section in pb)]
  nanoaodv7_cross_sections = get_cross_sections(nanoaodv7_files)
  # Organize cross sections
  #map_nanoaodv7_cross_sections[filename][year] = (path, cross section in pb)
  map_nanoaodv7_cross_sections = {}
  for ifile, file_info in enumerate(nanoaodv7_files):
    filename = file_info[0]
    year = file_info[1]
    cross_section = nanoaodv7_cross_sections[ifile][1]
    map_nanoaodv7_cross_sections[os.path.basename(filename)] = {year:[filename, cross_section]}
  # Write to log
  cross_section_log = open(cross_section_log_name,'w')
  for filename in sorted(map_nanoaodv7_cross_sections):
    for year in sorted(map_nanoaodv7_cross_sections[filename]):
      path = map_nanoaodv7_cross_sections[filename][year][0]
      cross_section = map_nanoaodv7_cross_sections[filename][year][1]
      cross_section_log.write("filepath: "+path+" year: "+str(year)+" cross-section: "+str(cross_section)+" pb\n")

  # Print cross sections for NanoAODv9 mc files
  # Get all NanoAODv9 filenames
  nanoaodv9_folders = ["/net/cms17/cms17r0/pico/NanoAODv9/nano/2016/mc", 
                       "/net/cms17/cms17r0/pico/NanoAODv9/nano/2017/mc",
                       "/net/cms17/cms17r0/pico/NanoAODv9/nano/2018/mc",
                       ]
  nanoaodv9_files = get_nanoaod_files(nanoaodv9_folders)
  #nanoaodv9_cross_sections = [(path, cross section in pb)]
  nanoaodv9_cross_sections = get_cross_sections(nanoaodv9_files)
  # Organize cross sections
  #map_nanoaodv9_cross_sections[filename][year] = (path, cross section in pb)
  map_nanoaodv9_cross_sections = {}
  for ifile, file_info in enumerate(nanoaodv9_files):
    filename = file_info[0]
    year = file_info[1]
    cross_section = nanoaodv9_cross_sections[ifile][1]
    map_nanoaodv9_cross_sections[os.path.basename(filename)] = {year:[filename, cross_section]}
  # Write to log
  for filename in sorted(map_nanoaodv9_cross_sections):
    for year in sorted(map_nanoaodv9_cross_sections[filename]):
      path = map_nanoaodv9_cross_sections[filename][year][0]
      cross_section = map_nanoaodv9_cross_sections[filename][year][1]
      cross_section_log.write("filepath: "+path+" year: "+str(year)+" cross-section: "+str(cross_section)+" pb\n")

  cross_section_log.close()
