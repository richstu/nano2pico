#!/usr/bin/env python
#script used to quickly run tests after nano2pico changes

import subprocess

def print_and_run(command):
  '''
  Given string command, print and run it through bash
  '''
  print(command)
  subprocess.call(command.split())

#dirs = ['2016_singleele','2016_dieleleg1','2016_dieleleg2','2016apv_singleele','2016apv_dieleleg1','2016apv_dieleleg2','2017_singleele','2017_dieleleg1','2017_dieleleg2','2018_singleele','2018_dieleleg1','2018_dieleleg2']
#folders = ['2016preVFP_UL','2016preVFP_UL','2016preVFP_UL','2016postVFP_UL','2016postVFP_UL','2016postVFP_UL','2017_UL','2017_UL','2017_UL','2018_UL','2018_UL','2018_UL']
#names = ['ele27','ele23','ele12','ele27','ele23','ele12','ele32','ele23','ele12','ele32','ele23','ele12']
#for i in range(len(dirs)):
#  cmd = ('./scripts/convert_correctionlib.py -i eg_'+dirs[i]+'.txt -f egamma_text -o data/zgamma/'+folders[i]+'/trigeff_'+names[i]).split()
#  subprocess.call(cmd)

#muon_dirs_names =  [
#        ['2016pre_singlemu',  '2016preVFP_UL', 'mu24','NUM_IsoMu24_or_IsoTkMu24_DEN_HToZGamma_SignalMuons_abseta_pt'],
#        ['2016post_singlemu', '2016postVFP_UL','mu24','NUM_IsoMu24_or_IsoTkMu24_DEN_HToZGamma_SignalMuons_abseta_pt'],
#        ['2017_singlemu',     '2017_UL',       'mu27','NUM_IsoMu27_DEN_HToZGamma_SignalMuons_abseta_pt'],
#        ['2018_singlemu',     '2018_UL',       'mu24','NUM_IsoMu24_DEN_HToZGamma_SignalMuons_abseta_pt'],
#        ['2016pre_dimuleg8',  '2016preVFP_UL', 'mu8','NUM_Mu8leg_DEN_HToZGamma_SignalMuons_abseta_pt'],
#        ['2016post_dimuleg8', '2016postVFP_UL','mu8','NUM_Mu8leg_DEN_HToZGamma_SignalMuons_abseta_pt'],
#        ['2017_dimuleg8',     '2017_UL',       'mu8','NUM_Mu8leg_DEN_HToZGamma_SignalMuons_abseta_pt'],
#        ['2018_dimuleg8',     '2018_UL',       'mu8','NUM_Mu8leg_DEN_HToZGamma_SignalMuons_abseta_pt'],
#        ['2016pre_dimuleg17', '2016preVFP_UL', 'mu17', 'NUM_Mu17leg_DEN_HToZGamma_SignalMuons_abseta_pt'],
#        ['2016post_dimuleg17','2016postVFP_UL','mu17', 'NUM_Mu17leg_DEN_HToZGamma_SignalMuons_abseta_pt'],
#        ['2017_dimuleg17',    '2017_UL',       'mu17', 'NUM_Mu17leg_DEN_HToZGamma_SignalMuons_abseta_pt'],
#        ['2018_dimuleg17',    '2018_UL',       'mu17', 'NUM_Mu17leg_DEN_HToZGamma_SignalMuons_abseta_pt']
#        ]
#for muon_name in muon_dirs_names:
#  cmd = ('./scripts/convert_correctionlib.py -i mu_'+muon_name[0]+'.root -f muon_root -o data/zgamma/'+muon_name[1]+'/trigeff_'+muon_name[2]+' -n '+muon_name[3]).split()
#  print(cmd)
#  subprocess.call(cmd)

#indir = '/net/cms17/cms17r0/pico/NanoAODv9/nano/2016APV/mc/'
#infiles = ['GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v1__40000__2948F584-1149-B14C-BCC7-ADFFAADCA068.root']

#indir = '/net/cms17/cms17r0/pico/NanoAODv9/nano/2016/signal/' 
#infiles = ['GluGluHToZG_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL16NanoAODv9__900.root']

#indir = '/net/cms17/cms17r0/pico/NanoAODv9/nano/2017/signal/'
#infile = 'GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1__70000__93B37B9F-16D5-6F4F-B920-3B9C682CA8A8.root'

#indir = '/net/cms17/cms17r0/pico/NanoAODv9/nano/2017/signal/'
#infiles = ['GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1__70000__93B37B9F-16D5-6F4F-B920-3B9C682CA8A8.root']
#           'GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1__70000__AED2B035-3C41-E044-92B7-5A908898EE0F.root',
#           'GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1__70000__E3D30224-C63D-CC48-9259-B0BE9FED9BB1.root',
#           'GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1__70000__E794DDC9-6560-064E-B36C-4E0D9DEF8556.root']

#indir = '/net/cms17/cms17r0/pico/NanoAODv9/nano/2017/mc/'
#infiles = ['DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v2__250000__06626381-A2F9-9E42-85B7-7A4B229C3FF2.root']

#indir = '/net/cms17/cms17r0/pico/NanoAODv9/nano/2017/mc/'
#infiles = ['TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1__50000__DF68DFBA-FC52-7248-8837-1073A83553A7.root']

#indir = '/net/cms11/cms11r0/pico/NanoAODv9/nano/2017/data/'
#infiles = ['DoubleMuon__Run2017F__UL2017_MiniAODv2_NanoAODv9-v1__120000__D4865BB4-A044-D84B-923B-9F78A8F669CD.root']

#indir = '/net/cms17/cms17r0/pico/NanoAODv9/nano/2018/mc/'
#infiles = ['GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018_realistic_v16_L1v1-v1__30000__743CEC74-C563-BB45-8248-8D571110CB02.root']
#infiles = ['DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018_realistic_v16_L1v1-v2__230000__917DE446-21B3-F248-AD70-6DA45CEC8914.root']

#indir = '/net/cms11/cms11r0/pico/NanoAODv11/nano/2022/mc/'
#infiles = ['DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8__Run3Summer22NanoAODv11__126X_mcRun3_2022_realistic_v2-v1__2560000__f2cc3f7c-3b6c-4035-b9ce-ea71117214c9.root']

#indir = '/net/cms11/cms11r0/pico/NanoAODv11/nano/2022/data/'
#infiles = ['EGamma__Run2022C__ReRecoNanoAODv11-v1__2540000__21ad7fb9-3fbd-46a8-b0c2-1bb6ab16aaa3.root']

#indir = '/net/cms11/cms11r0/pico/NanoAODV11p9/nano/2023/data/'
#infiles = ['EGamma0__Run2023C__PromptNanoAODv11p9_v1-v1__70000__674bf7c8-cca2-4f82-a9ec-5b60384b6f07.root']

#input_files = ['/net/cms11/cms11r0/pico/NanoAODv9/nano/2017/mc/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1__2510000__AF3EC52D-5109-A94D-876A-AFCC224E3944.root']
#input_files = ['/net/cms11/cms11r0/pico/NanoAODv9/nano/2018/mc/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018_realistic_v16_L1v1_ext1-v1__50000__2BFD0669-8DD8-2041-8262-08BC1EEB10F0.root']
#input_files = ['/net/cms11/cms11r0/pico/NanoAODv9/nano/2016APV/mc/WWW_4F_DiLeptonFilter_TuneCP5_13TeV-amcatnlo-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v2__50000__E8090E21-BF9C-AF4F-BB7C-582AFEBAA1BE.root']
#input_files = ['/net/cms11/cms11r0/pico/NanoAODv9/nano/2016/mc/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8__RunIISummer20UL16NanoAODv9__106X_mcRun2_asymptotic_v17_ext1-v1__120000__6A9C1237-C169-424F-B8AE-C91F18D3F040.root']
#input_files = ['/net/cms11/cms11r0/pico/NanoAODv9/nano/2017/mc/WZ_TuneCP5_13TeV-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1__130000__4563AB48-B968-8E40-A2B6-CC2267F666FD.root']
#input_files = ['/net/cms11/cms11r0/pico/NanoAODv9/nano/2017/mc/ZZ_TuneCP5_13TeV-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1__70000__FD2A8987-A9D3-A049-8ECE-CE7E52D922BA.root']
#input_files = ['/net/cms11/cms11r0/pico/NanoAODv9/nano/2018/mc/EWKZ2Jets_ZToLL_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018_realistic_v16_L1v1-v2__110000__81CD1433-385A-F547-94F0-A9C1AD2B1F54.root']
#input_files = ['/net/cms11/cms11r0/pico/NanoAODv9UCSB1/nano/2017/mc/ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1__045.root']

#do a data and MC file file of each era and data/mc type
input_files = [
               '/net/cms11/cms11r0/pico/NanoAODv9/nano/2016APV/data/DoubleEG__Run2016B__ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v3__2520000__B45617BE-6146-6546-8D17-72C1998A19A0.root',
               '/net/cms11/cms11r0/pico/NanoAODv9/nano/2016APV/mc/ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v1__30000__F995D7AB-9F03-D246-AE5C-3B084629CEE8.root',
               '/net/cms11/cms11r0/pico/NanoAODv9/nano/2016/data/DoubleMuon__Run2016G__UL2016_MiniAODv2_NanoAODv9-v2__2430000__D1989538-A8AE-854E-9F7B-6638D5D45817.root',
               '/net/cms11/cms11r0/pico/NanoAODv9/nano/2016/mc/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16NanoAODv9__106X_mcRun2_asymptotic_v17-v1__30000__E94099CB-E3C8-EB43-B182-4EA89C7C5411.root',
               '/net/cms11/cms11r0/pico/NanoAODv9/nano/2017/data/DoubleEG__Run2017C__UL2017_MiniAODv2_NanoAODv9-v1__270000__528C245A-DC66-7140-9F91-23C029849439.root',
               '/net/cms11/cms11r0/pico/NanoAODv9/nano/2017/mc/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1__2510000__AF3EC52D-5109-A94D-876A-AFCC224E3944.root',
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

#input_files = ['/net/cms17/cms17r0/pico/NanoAODv9/nano/2017/signal/GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1__70000__93B37B9F-16D5-6F4F-B920-3B9C682CA8A8.root']

#input_files = ['/net/cms11/cms11r0/pico/NanoAODv9/nano/2018/mc/ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018_realistic_v16_L1v1-v1__2520000__3363228C-76B3-9245-AEE9-8644875EE6D0.root']

#input_files = ['/net/cms11/cms11r0/pico/NanoAODv9/nano/2016APV/mc/ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16NanoAODAPVv9__106X_mcRun2_asymptotic_preVFP_v11-v1__30000__F995D7AB-9F03-D246-AE5C-3B084629CEE8.root']

#input_files = ['/net/cms11/cms11r0/pico/NanoAODv12/nano/2023/mc/DYGto2LG-1Jets_MLL-50_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8__Run3Summer23NanoAODv12__130X_mcRun3_2023_realistic_v15-v4__120000__14592dd6-3aa5-4e8c-b756-984c8d35392c.root']

for input_file in input_files:
  path_pos = input_file.rfind('/')
  indir = input_file[:path_pos]
  infiles = [input_file[path_pos+1:]]

  #process nano
  for infile in infiles:
    #cmd = './run/process_nano.exe --in_file '+infile+' --in_dir '+indir+' --out_dir out/zgamma/'
    cmd = './run/process_nano.exe --in_file '+infile+' --in_dir '+indir+' --out_dir out/zgamma/ --nent 1000'
    print_and_run(cmd)
  ##merge corrections
  #for infile in infiles:
  #  cmd = './run/merge_corrections.exe out/zgamma/corrections/corr_'+infile
  #  cmd += ' out/zgamma/wgt_sums/wgt_sums_'+infile
  #  print_and_run(cmd)
  ##apply corrections
  #for infile in infiles:
  #  cmd = './run/apply_corrections.exe --in_file raw_pico_'+infile+' --in_dir out/zgamma/raw_pico/ --corr_file corr_'+infile
  #  print_and_run(cmd)

##clean output
#for dirname in ['raw_pico','unskimmed','corrections','wgt_sums']:
#  print_and_run('pwd')
#  cmd = 'rm out/zgamma/'+dirname+'/*.root'
#  print_and_run(cmd)
