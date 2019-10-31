#!/bin/env python
import os, argparse

if __name__ == '__main__':

  parser = argparse.ArgumentParser(description="Submits batch jobs to split mass points of signal NanoAOD files",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("-i","--in_dir", default="/mnt/hadoop/pico/NanoAODv5/nano/2016/TChiHH_unsplit/",
                      help="Directory where the unsplit NanoAOD files are")
  parser.add_argument("-p","--target_dir", default="/mnt/hadoop/pico/NanoAODv5/nano/2016/TChiHH",
                      help="Determines the output folder.")
  parser.add_argument("-d","--dataset_filenames", default="SMS-TChiHH_HToBB_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16NanoAODv5_PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1_*.root",
                      help="File with list of commands for batch system.")
  parser.add_argument("-o","--out_cmd_filename", default="cmds.py",
                      help="File with list of commands for batch system.")
  args = parser.parse_args()
 
  source_directory = args.in_dir
  target_directory = args.target_dir

  if not os.path.exists(target_directory):
    os.makedirs(target_directory)
  
  out_string = '''#!/bin/env python
source_directory = "'''+source_directory+'''"
target_directory = "'''+target_directory+'''"
mass_points = [127, 150, 175, 
               200, 225, 250, 275, 
               300, 325, 350, 375, 
               400, 425, 450, 475, 
               500, 525, 550, 575, 
               600, 625, 650, 675, 
               700, 725, 750, 775, 
               800, 825, 850, 875, 
               900, 925, 950, 975, 
               1000, 1025, 1050, 1075,
               1100, 1125, 1150, 1175,
               1200, 1225, 1250, 1275,
               1300, 1325, 1350, 1375,
               1400, 1425, 1450, 1475,
               1500
               ]

for mass in mass_points:
  print("'''+os.getcwd()+'''/scripts/skim_file.py -m "+str(mass)+" -i \\""+source_directory+"'''+args.dataset_filenames+'''\\" -o "+target_directory)
'''

  with open(args.out_cmd_filename, 'w') as out_cmd_file:
    out_cmd_file.write(out_string)

  os.chmod(args.out_cmd_filename, 0o755)
  print("To generate job json and submit jobs do: ")
  print('convert_cl_to_jobs_info.py '+args.out_cmd_filename+' split_mass_points.json')
  print('auto_submit_jobs.py split_mass_points.json -c scripts/check_apply_corrections_job.py')
