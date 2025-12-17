from subprocess import run

BASE_DIR = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_redwood_v0_systsignal/'
BASE_DIR_R3 = '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0_systsignal/'
#years = ['2016APV','2016','2017','2018']
years = ['2018']
subdirs = ['raw_pico','wgt_sums','corrections','unskimmed','skim_ll',
           'skim_llg','merged_zgmc_ll','merged_zgmc_llg']

if __name__=='__main__':
  for year in years:
    for subdir in subdirs:
      #run((f'rm {BASE_DIR}{year}/mc/{subdir}/*HToMuMu*.root').split())
      #run((f'rm {BASE_DIR}{year}/mc/{subdir}/*WplusH*M-125*.root').split())
      #print(f'rm {BASE_DIR}{year}/mc/{subdir}/*WminusH*M-125*.root')
      print(f'rm {BASE_DIR}{year}/mc/{subdir}/*ZH_HToZG*ZToLL*M-125_TuneCP5*.root')
      #print(f'rm {BASE_DIR}{year}/mc/{subdir}/*GluGluHToZG_ZToLL_M-125_TuneCP5*.root')
      #print(f'rm {BASE_DIR_R3}{year}/mc/{subdir}/*GluGluHtoZG_Zto2L_M-125_TuneCP5*.root')
