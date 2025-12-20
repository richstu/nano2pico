#!/usr/bin/env python
#script used to quickly run tests after nano2pico changes

import subprocess

def print_and_run(command):
  '''
  Given string command, print and run it through bash
  '''
  print(command)
  subprocess.call(command.split())

#do a data and MC file file of each era and data/mc type
input_files = [
               '/net/cms11/cms11r0/pico/NanoAODv9/nano/2016APV/data/DoubleEG__Run2016B__ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v3__2520000__B45617BE-6146-6546-8D17-72C1998A19A0.root',
               '/net/cms11/cms11r0/pico/NanoAODv9/nano/2016APV/mc/ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v1__30000__F995D7AB-9F03-D246-AE5C-3B084629CEE8.root',
               '/net/cms11/cms11r0/pico/NanoAODv9/nano/2016/data/DoubleMuon__Run2016G__UL2016_MiniAODv2_NanoAODv9-v2__2430000__D1989538-A8AE-854E-9F7B-6638D5D45817.root',
               '/net/cms11/cms11r0/pico/NanoAODv9/nano/2016/mc/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16NanoAODv9__106X_mcRun2_asymptotic_v17-v1__30000__E94099CB-E3C8-EB43-B182-4EA89C7C5411.root',
               '/net/cms11/cms11r0/pico/NanoAODv9/nano/2017/data/DoubleEG__Run2017C__UL2017_MiniAODv2_NanoAODv9-v1__270000__528C245A-DC66-7140-9F91-23C029849439.root',
               '/net/cms17/cms17r0/pico/NanoAODv9/nano/2017/mc/GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1__70000__93B37B9F-16D5-6F4F-B920-3B9C682CA8A8.root',
               '/net/cms11/cms11r0/pico/NanoAODv9/nano/2018/data/EGamma__Run2018D__UL2018_MiniAODv2_NanoAODv9-v3__250000__973B1E8F-106B-1F44-92C8-3A2926E6761C.root',
               '/net/cms11/cms11r0/pico/NanoAODv9/nano/2018/mc/ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018_realistic_v16_L1v1-v1__2520000__3363228C-76B3-9245-AEE9-8644875EE6D0.root',
               '/net/cms11/cms11r0/pico/NanoAODv12/nano/2022/data/DoubleMuon__Run2022C__22Sep2023-v1__50000__bbfa5e3e-9080-43c4-ad61-f335ac9171fc.root',
               '/net/cms11/cms11r0/pico/NanoAODv12/nano/2022/mc/DYGto2LG-1Jets_MLL-50_PTG-50to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8__Run3Summer22NanoAODv12__130X_mcRun3_2022_realistic_v5-v3__2550000__e2050779-e2bb-4702-88a0-976643ec1618.root',
               '/net/cms11/cms11r0/pico/NanoAODv12/nano/2022EE/data/EGamma__Run2022F__22Sep2023-v1__2540000__26e202f6-027d-42c1-ba51-295dc6dcbfdf.root',
               '/net/cms11/cms11r0/pico/NanoAODv12/nano/2022EE/mc/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8__Run3Summer22EENanoAODv12__130X_mcRun3_2022_realistic_postEE_v6-v2__30000__3e44937a-8f96-41bc-abb5-c50ea9eb6fdf.root',
               '/net/cms11/cms11r0/pico/NanoAODv12/nano/2023/data/Muon1__Run2023C__22Sep2023_v1-v1__50000__a927c73e-4a29-4630-8814-7151d5d98cb9.root',
               '/net/cms11/cms11r0/pico/NanoAODv12/nano/2023/mc/GluGluHtoZG_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8__Run3Summer23NanoAODv12__130X_mcRun3_2023_realistic_v15-v3__50000__00f0d80c-0a3a-4243-8b81-0321d5f9c36e.root',
               '/net/cms11/cms11r0/pico/NanoAODv12/nano/2023BPix/data/EGamma0__Run2023D__22Sep2023_v2-v1__40000__9b38e671-508e-43b8-8486-9b341bfc2cac.root',
               '/net/cms11/cms11r0/pico/NanoAODv12/nano/2023BPix/mc/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8__Run3Summer23BPixNanoAODv12__130X_mcRun3_2023_realistic_postBPix_v2-v3__2550000__844b724c-9c84-4ae3-8b91-04d6d1d84212.root'
               ]

#input_files = ['/net/cms17/cms17r0/pico/NanoAODv9/nano/2017/mc/GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1__70000__93B37B9F-16D5-6F4F-B920-3B9C682CA8A8.root']

#input_files = ['/net/cms11/cms11r0/pico/NanoAODv9/nano/2018/mc/ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018_realistic_v16_L1v1-v1__2520000__3363228C-76B3-9245-AEE9-8644875EE6D0.root']

#input_files = ['/net/cms11/cms11r0/pico/NanoAODv9/nano/2016APV/mc/ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v1__30000__F995D7AB-9F03-D246-AE5C-3B084629CEE8.root']

#input_files = ['/net/cms11/cms11r0/pico/NanoAODv12/nano/2023/mc/DYGto2LG-1Jets_MLL-50_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8__Run3Summer23NanoAODv12__130X_mcRun3_2023_realistic_v15-v4__120000__14592dd6-3aa5-4e8c-b756-984c8d35392c.root']

#input_files = ['/net/cms11/cms11r0/pico/NanoAODv12/nano/2022EE/mc/DYGto2LG-1Jets_MLL-50_PTG-200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8__Run3Summer22EENanoAODv12__130X_mcRun3_2022_realistic_postEE_v6-v4__30000__f5929cb0-6abf-4d15-a539-82dd8a25e07e.root']

for input_file in input_files:
  path_pos = input_file.rfind('/')
  indir = input_file[:path_pos]
  infiles = [input_file[path_pos+1:]]

  #process nano
  for infile in infiles:
    #cmd = './run/process_nano.exe --in_file '+infile+' --in_dir '+indir+' --out_dir out/zgamma/'
    cmd = './run/process_nano.exe --in_file '+infile+' --in_dir '+indir+' --out_dir out/zgamma/ --nent 1000'
    print_and_run(cmd)
  if not 'data' in input_file:
    #merge corrections
    for infile in infiles:
      cmd = './run/merge_corrections.exe out/zgamma/corrections/corr_'+infile
      cmd += ' out/zgamma/wgt_sums/wgt_sums_'+infile
      print_and_run(cmd)
    #apply corrections
    for infile in infiles:
      cmd = './run/apply_corrections.exe --in_file raw_pico_'+infile+' --in_dir out/zgamma/raw_pico/ --corr_file corr_'+infile
      print_and_run(cmd)

##clean output
#for dirname in ['raw_pico','unskimmed','corrections','wgt_sums']:
#  print_and_run('pwd')
#  cmd = 'rm out/zgamma/'+dirname+'/*.root'
#  print_and_run(cmd)
