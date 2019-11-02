# nano2pico

Utility package for converting CMS NanoAOD to analysis-ready ntuples called "pico".

## Environment setting

Use one of the servers supporting CMSSW, e.g. cms1,cms3,cms4,cms5...

~~~~bash
    . /cvmfs/cms.cern.ch/cmsset_default.sh
    cd /net/cms29/cms29r0/pico/CMSSW_10_2_11_patch1/src
    eval `scramv1 runtime -sh`
    cd -
    git clone --recurse-submodules https://github.com/richstu/nano2pico
    source set_env.sh
~~~~

## Latest production

Variables stored in the pico can be seen in [variables/pico](variables/pico). For an overview of the available branches, see the dedicated section at the bottom of this README.

The latest Higgsino production nicknamed "Angeles" can be found here:

~~~~bash
  /net/cms2/cms2r0/pico/NanoAODv5/higgsino_angeles/2016/mc/
  /net/cms2/cms2r0/pico/NanoAODv5/higgsino_angeles/2016/TChiHH/
~~~~

To see the sizes and number of files, one can do:

~~~~bash
  ./scripts/count_root_files.py -f /net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/
~~~~

## What is nano2pico?

This package is used to do the Nano -> pico conversion in three steps in order to allow parallelizing the production at the sub-dataset level:
  1. All variables and event weights (except normalization) are calculated using [src/process_nano.cxx](src/process_nano.cxx). This step also keeps a tally of the weights for all events in the file being run over as input to the next step.
  2. The sums of weights from step 1 are further aggregate to get the total per dataset. A correction is then calculated to ensure that the weights do not change the total expected number of events for the dataset. This is done in [src/merge_corrections.cxx](src/merge_corrections.cxx). The luminosity normalization weight `w_lumi` to be applied to get the yield in 1fb-1 is also calculated for each dataset in this step.
  3. The `raw_pico` files from step 1 are corrected by the per-dataset correction factors derived in step 2 and written to the `unskimmed` folder.

At this point, various skims can be made as defined in [scripts/skim_file.py](scripts/skim_file.py).

## Interactive test usage

Define some paths, e.g.:

~~~~bash
export INDIR=/net/cms29/cms29r0/pico/NanoAODv5/nano/2016/TChiHH/
export INFILE=SMS-TChiHH_mChi-1000_mLSP-1_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv5__PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1.root
~~~~

Step 1. Produce raw pico ntuple from a nano input file, adding `--isFastsim` and `--isData` if applicable:

~~~~bash
./compile.sh && ./run/process_nano.exe --in_file $INFILE --in_dir $INDIR --out_dir out/ --nent 10000
~~~~

Note that for interactive jobs, you need to ensure that the output directory subfolders `wgt_sums` and `raw_pico` exist.

:bangbang: Code functionality relies on the input NanoAOD filename! Specifically, `INFILE` is parsed for:

* flag `isData = infile.Contains("Run201") ? true : false;`
* flag `isFastsim = infile.Contains("Fast") ? true : false;`
* variable `year = infile.Contains("RunIISummer16") ? 2016 : (infile.Contains("RunIIFall17") ? 2017 : 2018)`
* output branch `type` is set based on the presence of dataset name substrings (see event_tools.cpp)
* branches related on ISR also depend on the presence of dataset name substrings

Step 2. For each dataset, add up the sums of weights obtained for each file in step 1 and calculate the corrections needed to normalize each individual weight as well as the total weight. Note that the order of options is fixed with the arguments after the first being the input files. This is to allow arbitrary number of input files. Note that again functionality depends on the naming, e.g. correction file name is used to decide what cross-section to use.

~~~~bash
./compile.sh && ./run/merge_corrections.exe out/corrections/corr_$INFILE out/wgt_sums/wgt_sums_$INFILE
~~~~

Step 3. Using the pico file from step 1 and the corrections file from step 2 as input, we can renormalize the weight branches as follows:

~~~~bash
./compile.sh && ./run/apply_corrections.exe --in_file raw_pico_$INFILE --in_dir out/raw_pico/ --corr_file corr_$INFILE
~~~~

## Batch system usage

### Step 0. Setup environment

  source set_env.sh

### Step 1. Converting Nano to Pico:
Generate a python file that prints the commands to be run in the batch (input for the queue system):

~~~~bash 
./scripts/write_process_nano_cmds.py --in_dir /mnt/hadoop/pico/NanoAODv5/nano/2016/mc/ \
                                      --production higgsino_angeles \
                                      --dataset_list datasets/higgsino_2016_mc_dataset_list.txt
~~~~

or for signal, just specify the appropriate input folder and omit the `--dataset_list` argument to run on all files in the input folder.

This produces the commands in `cmds.py`. You can perform a last check by running one of the commands interactively. Next, submit the jobs to the batch system. Note the -c option which allows to attach a script that compares the input and output number of entries when each job is done. Note the check can be performed later if one needs to detach the session. Alternatively, this command can be started in screen:

~~~~bash 
convert_cl_to_jobs_info.py cmds.py higgsino_angeles.json
auto_submit_jobs.py higgsino_angeles.json -c scripts/check_process_nano_job.py
~~~~

To check whether the jobs were successful later on do something like, where `auto_higgsino_angeles.json` is generated earlier by the `auto_submit_jobs.py` command:

~~~~bash 
check_jobs.py auto_higgsino_angeles.json -c scripts/check_process_nano_job.py
~~~~

This command will result in `checked_auto_higgsino_angeles.json`, which can then be used to resubmit failed jobs if any:

~~~~bash 
select_resubmit_jobs.py checked_auto_higgsino_angeles.json -c scripts/check_process_nano_job.py 
auto_submit_jobs.py resubmit_checked_auto_higgsino_angeles.json -c scripts/check_process_nano_job.py 
~~~~

### Step 2. Merge sums of weights

For example:

~~~~bash 
./scripts/merge_corrections.py --wgt_dir /net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/wgt_sums/ \
                               --corr_dir /net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/corrections/ 
~~~~

### Step 3. Submit the weight correction jobs

To generate the commands use:

~~~~bash 
./scripts/write_apply_corrections_cmds.py --in_dir /net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/raw_pico/
~~~~

Follow similar process as in Step 1 to submit the commands as batch jobs. 

## Calculating b-tagging efficiencies

Use `parameterize_efficiency.cxx`, giving the directory with all the MC files and the year as arguments. Below is an example run for 2016 MC.

~~~~bash
./compile.sh && ./run/parameterize_efficiency.exe -i /mnt/hadoop/jbkim/2019_09_30/2016/mc/ -y 2016
~~~~


## Description of pico branches

These refer to the branches obtained with `ZGamma = false`, i.e. higgsino production!

:blue_book: Documentation of the Nano variables used as input throughout the code can be found [here](https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc102X_doc.html).

####   Global
* `run, lumiblock, event` - as expected
* `type` - integer encoding the physics process, see [here](https://github.com/richstu/nano2pico/blob/5e62c553fb306f6c1f27bccafee037fb939c0f75/src/event_tools.cpp#L124).
* `stitch` - include this variable in order to run on an inclusive sample together with an overlapping slice in a different dataset, e.g. stitch = false for events with GenMET > 150 in the inclusive TTJets sample, in order to remove them when using the inclusive sample together with the deidicated *genMET-150* samples, see [here](https://github.com/richstu/nano2pico/blob/f4b99417bd65c134796b703552522a7de5429f19/src/event_tools.cpp#L31-L45).
* `npv` - number of reconstructed PV
* `ht` - sum of pt of jets not associated with a lepton
* `met, met_phi, met_calo, met_tru, met_tru_phi` - as expected
* `mt` - transverse mass, only calculated for nlep==1
* `mt_tru` - transverse mass at truth level, only calculated for ntrulep==1

####   Higgsino variables
Using the 4-jet with highest DeepCSV, calculate the higgsino variables for the three possible pairings. The 0th index stores the pairing with smalled Delta m
* `hig_cand_dm` - Mass difference between the two Higgs candidates
* `hig_cand_am` - Average mass between the two Higgs candidates
* `hig_cand_drmax` - Max opening angle between the two b jets out of the two Higgs candidates.

Same variables using the 4-jet with highest DeepFlavour discriminant value are stored in:
* `hig_df_cand_*` 

* `low_dphi` - require dPhi(jet, MET) be less than 0.5 for jets 1,2 and less than 0.3 for jets 3,4

####   Jets

Filled in `jet_producer`:
* `nbl, nbm, nbt` - number of loose, medium and tight tagged jets according to DeepCSV tagger
* `nbdfl, nbdfm, nbdft` - number of loose, medium and tight tagged jets according to DeepFlavour tagger
* `njet` - number of jets that pass the pt and eta cuts and do not overlap with a _signal_ lepton
* `jet_*` - basic jet related variables and also:
  * `jet_h1d, jet_h2d` - booleans indicating whether this jet is one of the two jets in Higgs 1 or Higgs 2 of the 0th pair of Higgs candidates stored in `hig_cand_*`
  * `jet_fjet_idx` - index of any fat jets within 0.8

* `int nfjet` - number of AK8 jets that pass the pt and eta cuts
* `fjet_*` - basic fat jet related variables and also:
  * `fjet_deep_md_hbb_btv` - Mass-decorrelated Deep Double B, H->bb vs QCD discriminator, endorsed by BTV
  * `fjet_deep_md_hbb_jme` - Mass-decorrelated DeepAk8, H->bb vs QCD discriminator, endorsed by JME

#### Leptons 

* `nlep = nel + nmu`
* `nvlep = nvel + nvmu`

Calculated in [mu_producer](src/mu_producer.cpp):
* `nmu` - number of muons satisfying all signal muon requirements for resolved Higgsino analysis
* `nvmu` - number of muons satisfying all veto muon requirements for resolved Higgsino analysis
* `mu_*` - variables of interest for all muons satisfying the veto id, eta and pt requirements, but no isolation requirements

Calculated in [el_producer](src/el_producer.cpp):
* `nel` - number of electrons satisfying all signal electron requirements for resolved Higgsino analysis
* `nvel` - number of electrons satisfying all veto electron requirements for resolved Higgsino analysis
* `el_*` - variables of interest for all electrons satisfying the veto id, eta and pt requirements, but no isolation 

Calculated in [dilep_producer](src/dilep_producer.cpp):
* `elel_*, mumu_*` - variables relating to the dilepton system (all combinations stored if more than 2 leptons)

#### Photons 

These are not really used in Higgsino, but just in case...Calculated in [photon_producer](src/photon_producer.cpp):
* `nphoton` - number of signal photons
* `photon_*` - photon variables

#### Tracks 

Calculated in [tk_producer](src/tk_producer.cpp):
* `ntk` - number of tracks passing criteria for resolved Higgsino analysis
* `tk_*` - track variables

#### Quality 

* `pass_*` - recommended MET filters
* `pass_jets` - set to false if any of the jets fails loose ID
* `pass` - combination of all required filters and `pass_jets`

#### Truth 

* `mc_*` - information for a set of the generator particles in the hard process
* `ntrumu,ntruel,ntrutauh,ntrutaul` - # of true leptons of particular type, where tauh is hadronically decaying taus and taul is leptonically decaying taus
* `ntrulep = ntrumu + ntruel + ntrutaul` 
* `mprod, mlsp` - higgsino and lsp mass, with lsp mass always equal to one for the higgsino model

####   ISR

* `nisr` - number of ISR jets according to matching to truth, used for ISR reweighting used by the SUS PAG for strong production
* `isr_tru_*` - MC truth, hadronic recoil, used for ISR reweighting used by the SUS PAG for weak production

* `jetsys_*` - hadronic recoil, i.e. vector sum of all jets, used in V+jets ISR studies
* `jetsys_nob_*` - hadronic recoil, i.e. vector sum of all jets that are not b-tagged, used in 2L tt+jets ISR studies

#### Weights 

Calculated in [process_nano](src/process_nano.cxx) and then re-normalized in subsequent production steps:
* `weight` - product of all the individual weights below
* `w_lumi` - weight to be applied to get the expected yield in 1 fb-1. 
* `w_lep` - product of fullsim lepton SFs for leptons with > 20 GeV
* `w_fs_lep` - product of fastsim lepton SFs for leptons with > 20 GeV
* `w_btag` - product of fullsim and fastsim b-tag SFs if counting _medium tags only_, DeepCSV tagger
* `w_btag_df` - product of fullsim and fastsim b-tag SFs if counting _medium tags only_, DeepFlavour tagger
* `w_bhig` - product of fullsim and fastsim b-tag SFs accounting for _all three WPs_ for the DeepCSV tagger
* `w_bhig_df` - product of fullsim and fastsim b-tag SFs accounting for _all three WPs_ for the DeepFlavour tagger
* `w_isr` - 1., except for TTJets 2016 and signal, SUSY ISR reweighting
* `w_pu` - currently just set 1.
* `w_prefire` - currently just set 1.

#### Other

* `HLT_*` - trigger decisions
* `sys_*` - systematic variations of weights up=0, down=1
