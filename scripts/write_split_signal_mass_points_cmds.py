#!/bin/env python
import os, argparse
import ROOT

def get_2d_mass_points(signal_chain, pdgId_1, pdgId_2):
  signal_chain.SetEstimate(signal_chain.GetEntries()+1)
  signal_chain.Draw("MinIf$(GenPart_mass,GenPart_pdgId=="+str(pdgId_1)+"):MinIf$(GenPart_mass,GenPart_pdgId=="+str(pdgId_2)+")","", "goff")
  number_variables = signal_chain.GetSelectedRows()
  mass_array_1 = signal_chain.GetV1()
  mass_array_2 = signal_chain.GetV2()
  mass_points = set()
  '''
  #print(number_variables, signal_chain.GetEntries())
  for iVar in range(number_variables):
    #print (mass_array_1[iVar], mass_array_2[iVar])
    if (mass_array_1[iVar] == 127):
      mass_points.add((mass_array_1[iVar], int(round(mass_array_2[iVar]/25)*25)))
    else:
      mass_points.add((int(round(mass_array_1[iVar]/25)*25), int(round(mass_array_2[iVar]/25)*25)))
  print(mass_points)'''

  
  # for 4b samples:
  nlsp_min, nlsp_max, nlsp_step = 200, 1300, 100
  lsp_min, lsp_max, lsp_step = 100, 1100, 100

  cols = []        
  cols.append([150,1])
  for mx in range(nlsp_min, nlsp_max+1, nlsp_step):
      cols.append([mx,1])
      for my in range(lsp_min, lsp_max+1, lsp_step):
          if mx - my < 127:
              continue
          cols.append([mx, my])
      if mx>250:
          cols.append([mx, mx-250])
      if mx>150:
          cols.append([mx, mx-150])
  mpoints = cols

  for point in mpoints:
    mchi, mlsp = point[0], point[1]
    mass_points.add((mchi, mlsp))


  return sorted(mass_points)





if __name__ == '__main__':

  parser = argparse.ArgumentParser(description="Submits batch jobs to split mass points of signal NanoAOD files",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("-2","--two_dim", action="store_true")
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

  if args.two_dim:
    chain = ROOT.TChain('Events')
    print(args.in_dir+"/"+args.dataset_filenames)
    chain.Add(args.in_dir+"/"+args.dataset_filenames)
    # signal_id=1000023, lsp_id=1000022
    # mass_points = [(mass of signal, mass of lsp)]
    mass_points = get_2d_mass_points(chain, 1000023, 1000022)
    out_string = '''#!/bin/env python
source_directory = "'''+source_directory+'''/"
target_directory = "'''+target_directory+'''/"
mass_points = '''+str(mass_points)+'''
for mass_point in mass_points:
  print("'''+os.getcwd()+'''/scripts/skim_file.py -m "+str(mass_point[0])+" -l "+str(mass_point[1])+" -i '"+source_directory+"'''+args.dataset_filenames+'''' -o "+target_directory)
'''
  else:
    out_string = '''#!/bin/env python
source_directory = "'''+source_directory+'''/"
target_directory = "'''+target_directory+'''/"
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
  print("'''+os.getcwd()+'''/scripts/skim_file.py -m "+str(mass)+" -i '"+source_directory+"'''+args.dataset_filenames+'''' -o "+target_directory)
'''

  with open(args.out_cmd_filename, 'w') as out_cmd_file:
    out_cmd_file.write(out_string)

  os.chmod(args.out_cmd_filename, 0o755)
  print("To generate job json and submit jobs do: ")
  print('convert_cl_to_jobs_info.py '+args.out_cmd_filename+' '+os.path.splitext(args.out_cmd_filename)[0]+'.json')
  print('auto_submit_jobs.py '+os.path.splitext(args.out_cmd_filename)[0]+'.json -c jobscript_check.py -n cms1')
