# nano2pico

Utility package for converting CMS NanoAOD to analysis-ready ntuples called "pico". These have been used in the HH+MET and H to Zgamma analyses.

## Installation, environment, and compilation

If you are not on one of the UCSB servers, prerequisites that must be installed include CERN ROOT and optionally Scons for building nano2pico. Assuming your ROOT environment is active, you can build nano2pico after cloning this repository by running `scons` if you have Scons installed, or by running `compile.sh` otherwise.

If you are using one of the UCSB servers that supports CMSSW (e.g. cms1,cms3,cms4,cms5...), you can use the following commands to install nano2pico and set up the environment:

~~~~bash
# Setup git version for SL6.5
. /cvmfs/cms.cern.ch/cmsset_default.sh;cd /net/cms29/cms29r0/pico/CMSSW_10_2_11_patch1/src;eval `scramv1 runtime -sh`;cd -
# Clone git
git clone --recurse-submodules git@github.com:richstu/nano2pico.git
# If did not use recurse at clone, use following command: git submodule update --init --remote --recursive
# Setup environemnt
source set_env.sh
~~~~

You can then compile via `scons` or `compile.sh`.

## Productions

Variables stored in the pico can be seen in [variables/pico](variables/pico). For an overview of the available branches, see the dedicated section at the bottom of this README.

To see the sizes and number of files in all available productions do:

~~~~bash
  ./scripts/count_root_files.py -f /net/cms29/cms29r0/pico/NanoAODv5/
~~~~

## How does nano2pico work?

This package is used to do the Nano -> pico conversion in three steps in order to allow parallelizing the production at the sub-dataset level:

  1. All variables and event weights (except normalization) are calculated using [src/process_nano.cxx](src/process_nano.cxx). For Monte Carlo events, this step also keeps a tally of the weights for all events in the file being run over as input to the next step. For data events, only this step is needed.
  2. The sums of weights from step 1 are further aggregate to get the total per dataset. A correction is then calculated to ensure that the weights do not change the total expected number of events for the dataset. This is done in [src/merge_corrections.cxx](src/merge_corrections.cxx). The luminosity normalization weight `w_lumi` to be applied to get the yield in 1fb-1 is also calculated for each dataset in this step.
  3. The `raw_pico` files from step 1 are corrected by the per-dataset correction factors derived in step 2 and written to the `unskimmed` folder.

At this point, various skims can be made as defined in [scripts/skim_file.py](scripts/skim_file.py).

*Note:* The input path in which the input NanoAOD files are stored as well as the output path are analyzed to determine the behavior of nano2pico. To run with settings for the Higgs to Z gamma analysis, the 'out_dir' should contian "zgamma" in its name. To run on custom NanoAODv9 files, the input directory must contain "NanoAODv9UCSB" in its name.

## Interactive test usage

Define some paths, e.g.:

~~~~bash
export INDIR=/net/cms29/cms29r0/pico/NanoAODv5/nano/2016/TChiHH/
export INFILE=SMS-TChiHH_mChi-1000_mLSP-1_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv5__PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1.root
~~~~

Step 1. Make an output directory out/ with subdirectories `wgt_sums` and `raw_pico`. Produce raw pico ntuple from a nano input file:

~~~~bash
./compile.sh && ./run/process_nano.exe --in_file $INFILE --in_dir $INDIR --out_dir out/ --nent 10000
~~~~

:bangbang: Code functionality relies on the input NanoAOD filename! Specifically, `INFILE` is parsed for:

* flag `isData = infile.Contains("Run201") ? true : false;`
* flag `isFastsim = infile.Contains("Fast") ? true : false;`
* flag `isSignal = Contains(in_file, "TChiHH") || Contains(in_file, "T5qqqqZH") ? true : false;`
* variable `year = infile.Contains("RunIISummer16") ? 2016 : (infile.Contains("RunIIFall17") ? 2017 : 2018)`
* output branch `type` is set based on the presence of dataset name substrings (see event_tools.cpp)
* branches related on ISR also depend on the presence of dataset name substrings

Step 2. If you are using data, you are done! If you are using MC, for each dataset, add up the sums of weights obtained for each file in step 1 and calculate the corrections needed to normalize each individual weight as well as the total weight. Note that the order of options is fixed with the arguments after the first being the input files. This is to allow arbitrary number of input files. Note that again functionality depends on the naming, e.g. correction file name is used to decide what cross-section to use.
Make subdirectory `corrections` in `out`.

~~~~bash
./compile.sh && ./run/merge_corrections.exe out/corrections/corr_$INFILE out/wgt_sums/wgt_sums_$INFILE
~~~~

Step 3. Make subdirectory `unskimmed in `out`. Using the pico file from step 1 and the corrections file from step 2 as input, we can renormalize the weight branches as follows:

~~~~bash
./compile.sh && ./run/apply_corrections.exe --in_file raw_pico_$INFILE --in_dir out/raw_pico/ --corr_file corr_$INFILE
~~~~

## Batch system usage

nano2pico supports batch system usage to process many datasets in parallel. Currently, this system is only configured for use at UCSB, and must be modified if nano2pico is being run on another system.

While previous versions of nano2pico used the custom UCSB batch system, the current version now uses HTCondor. For this reason, jobs must be run on the cms11 server. To get permission to run jobs on HTCondor, please contact Jaebak.

### Step 0. Setup environment

  source set_env.sh

### Step 1. Converting Nano to Pico:

First, generate a text file containing the datasets in DAS format (this is produced by copy\_dataset) or the filenames to be processed, one per line. If you use filenames, you must add the argument `--list_format filename` when invoking `scripts/write_process_nano_cmds.py`.

Next, generate a python file that prints the commands to be run in the batch (input for the queue system):

~~~~bash 
./scripts/write_process_nano_cmds.py --in_dir /mnt/hadoop/pico/NanoAODv5/nano/2016/mc/ \
                                      --production higgsino_angeles \
                                      --dataset_list datasets/higgsino_2016_mc_dataset_list.txt
~~~~

or for signal, just specify the appropriate input folder and omit the `--dataset_list` argument to run on all files in the input folder.

To run on data, use `--list_format filename` in order to interpret the lines in the file passes to `--dataset_list` as a list of filenames with wildcards. For an example file, see [txt/datasets/higgsino_data_infile_list.txt](txt/datasets/higgsino_data_infile_list.txt). For example:

~~~~bash 
./scripts/write_process_nano_cmds.py --in_dir /net/cms29/cms29r0/pico/NanoAODv5/nano/2016/data/ \
                                     --production higgsino_humboldt \
                                     --dataset_list txt/datasets/higgsino_data_infile_list.txt \
                                     --list_format filename
~~~~

This produces the commands in `cmds.py`. You can perform a last check by running one of the commands interactively. Next, submit the jobs to the batch system. Note the -c option which allows to attach a script that compares the input and output number of entries when each job is done. Note the check can be performed later if one needs to detach the session. Alternatively, this command can be started in screen:

~~~~bash 
auto_submit_jobs.py process_nano_cmds.json -c scripts/check_process_nano_job.py
~~~~

For data use the scripts/check_data_process_nano_job.py like below

~~~~bash 
auto_submit_jobs.py process_nano_cmds.json -c scripts/check_data_process_nano_job.py
~~~~

If the above script is interrupted, one can check whether the jobs were successful later on by passing the json produced by auto_submit_job.py to check_jobs.py:

~~~~bash 
check_jobs.py auto_higgsino_angeles.json -c scripts/check_process_nano_job.py
~~~~

This command will result in `checked_auto_higgsino_angeles.json`, which can then be used to resubmit failed jobs if any:

~~~~bash 
select_resubmit_jobs.py checked_auto_higgsino_angeles.json -c scripts/check_process_nano_job.py 
auto_submit_jobs.py resubmit_checked_auto_higgsino_angeles.json -c scripts/check_process_nano_job.py 
~~~~

One can also resume the auto_submit_job.py like below

~~~~bash 
auto_submit_jobs.py auto_higgsino_angeles.json -c scripts/check_process_nano_job.py -o auto_higgsino_angeles.json
~~~~

If processing Monte Carlo, then proceed to steps 2 and 3 (MC). If processing data, proceed directly to step 4.

### Step 2 (MC). Merge sums of weights

For example:

~~~~bash 
./scripts/merge_corrections.py --wgt_dir /net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/wgt_sums/ \
                               --corr_dir /net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/corrections/ 
~~~~

### Step 3 (MC). Submit the weight correction jobs

To generate the commands use:

~~~~bash 
./scripts/write_apply_corrections_cmds.py --in_dir /net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/raw_pico/
~~~~

Follow similar process as in Step 1 to submit the commands as batch jobs. 

### Step 4. Making skims

It's recommended to start with a relatively inclusive skim which would then serve as the starting point for tighter skims to minimize total time spent on skimming. For example:

~~~~bash 
./scripts/write_skim_cmds.py --in_dir /net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/TChiHH/skim_higloose/ \
                             --skim_name higtight \
                             --tag apples
auto_submit_jobs.py skim_hightight_cmds_apples.json -c scripts/check_skim.py
~~~~

The skim names are defined in [scripts/skim_file.py](scripts/skim_file.py). If defining a new skim, please commit the definition!! This eliminates confusion of what is in various folders on disk later on.

The argument `--tag` is optional. It is used to differentiate the JSON files created by the queue system in case of running multiple skims of the same type. It will not affect the folder structure.

Use `--overwrite` to run over all files even if output already exists. Otherwise, restarting the process of batch submission will skip files that have already been processed. Note that if you just re-issue the `auto_submit_jobs.py` with the original json file WILL overwrite. To omit files with existing output re-start from this step.

Note that this step works also on slims produced by Step 5.

### Step 5. Making slims

Finally, one can remove branches that are not commonly used and merge all files pertaining to one dataset into a single file to further reduce size and speed up making plots. For example:

~~~~bash 
./scripts/write_slim_and_merge_cmds.py --in_dir /net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/skim_met150/ \
                                       --slim_name higmc
~~~~

Here the slim name must correspond to a txt file in the slim_rules folder, so in this example `txt/slim_rules/higmc.txt`. The file contains the list of branches to be dropped/kept.

Similarly to above, one can optionally use `--overwrite` or `--tag`.

### Step 6. Prior to DNN training: Prepare tree with DNN inputs (deprecated)

For the higgsino analysis, one can prepare a tree with all the necessary DNN inputs for either training or inference using the executable `make_higfeats.exe`, and in the batch system, e.g.:

~~~~bash 
./scripts/write_generic_cmds.py ./scripts/write_generic_cmds.py \
           -i /net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado/2017/mc/merged_higmc_higloose/ \
           -o /net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado/2017/mc/higfeats_higloose/ \
           -e ./run/make_higfeats.exe -t mc2017
~~~~

As usual, the tag is optional and only relevant for the filename of the resulting cmd file.

### Step 7. After DNN evaluation: Merge pico with DNN output (deprecated)

After training the DNN and evaluating its output for all samples of interest using the `diboson_ml` package, one can update the corresponding pico trees to add a new branch containing the DNN output. This relies on having the events in the same order, so one has to update the pico ntuples used as input to higfeats! Given it is rather quick, it's done interactively.

For now, copy the input folder just in case...

~~~~bash 
cp -r /net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/mc/merged_higmc_higloose/ \
      /net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/mc/mergednn_higmc_higloose/ 
./scripts/run_update_pico.py \
     --pico_dir /net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/mc/mergednn_higmc_higloose/ \
     --dnnout_dir /net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/mc/dnnout_higloose/
~~~~

## Calculating b-tagging efficiencies (deprecated)

Use `parameterize_efficiency.cxx`, giving the directory with all the MC files and the year as arguments. Below is an example run for 2016 MC.

~~~~bash
./compile.sh && ./run/parameterize_efficiency.exe -i /mnt/hadoop/jbkim/2019_09_30/2016/mc/ -y 2016
~~~~

## Description of pico branches

This section describes the event content (branches) of pico n-tuples.

:blue_book: Documentation of the Nano variables used as input throughout the code can be found [here](https://cms-nanoaod-integration.web.cern.ch/autoDoc).

####   Global

Filled by various files including [event_tools](src/event_tools.cpp), [jetmet_producer](src/jetmet_producer.cpp), and [process_nano](src/process_nano.cxx):
* `run, lumiblock, event` - as expected
* `type` - integer encoding the physics process, see `src/event_tools.cpp`.
* `stitch, stitch_*` - include `stitch` in order to run on an inclusive sample together with an overlapping slice in a different dataset, e.g. stitch = false for events with GenMET > 150 in the inclusive TTJets sample, in order to remove them when using the inclusive sample together with the deidicated *genMET-150* samples, see [here](https://github.com/richstu/nano2pico/blob/f4b99417bd65c134796b703552522a7de5429f19/src/event_tools.cpp#L31-L45). The full stitch variable is determined by the `stich_*` pieces.
* `is_overlap, use_event` - overlap removal similar to stitch, but as implemented for the H to Z gamma analysis. Mostly deals with overlap between inclusive samples (ex. Z to ll) and those with photons generated at matrix level. `use_event` will be true for events in one overlapping sample, but not the other
* `npv, npv_good` - number of (good) reconstructed PVs
* `rho` - pileup energy density calculated from all PF candidates without foreground removal
* `ht, ht5` - sum of pt of jets not associated with a lepton, including jets with |eta|<2.4 (|eta|<5)
* `mht, mht_phi` - negative vector sums of jets not associated with a lepton
* `met, met_phi, met_calo` - PF/calorimeter missing transverse momentum, as expected
* `mt` - transverse mass, only calculated for nlep==1

####   Jets

Filled in `jetmet_producer`:
* `nbl, nbm, nbt` - number of loose, medium and tight tagged jets according to DeepCSV tagger
* `nbdfl, nbdfm, nbdft` - number of loose, medium and tight tagged jets according to DeepFlavour tagger
* `njet` - number of jets that pass the pt and eta cuts and do not overlap with a _signal_ lepton
* `jet_*` - basic jet related variables and also:
  * `jet_h1d, jet_h2d` - booleans indicating whether this jet is one of the two jets in Higgs 1 or Higgs 2 of the 0th pair of Higgs candidates stored in `hig_cand_*`
  * `jet_fjet_idx` - index of any fat jets within 0.8

* `nfjet` - number of AK8 jets that pass the pt and eta cuts
* `fjet_*` - basic fat jet related variables and also:
  * `fjet_deep_md_hbb_btv` - Mass-decorrelated Deep Double B, H->bb vs QCD discriminator, endorsed by BTV
  * `fjet_deep_md_hbb_jme` - Mass-decorrelated DeepAk8, H->bb vs QCD discriminator, endorsed by JME

* `nsubfjet` - number of AK8 jet subjets
* `subfjet_*` - AK8 jet subjet variables

#### Leptons 

* `nlep = nel + nmu`
* `nvlep = nvel + nvmu`
* `lep_*` - signal leptons (electron and muon) variables

Calculated in [mu_producer](src/mu_producer.cpp):
* `nmu` - number of muons satisfying all signal muon requirements for resolved Higgsino analysis
* `nvmu` - number of muons satisfying all veto muon requirements for resolved Higgsino analysis
* `mu_*` - variables of interest for all muons satisfying the veto id, eta and pt requirements, but no isolation requirements

Calculated in [el_producer](src/el_producer.cpp):
* `nel` - number of electrons satisfying all signal electron requirements for resolved Higgsino analysis
* `nvel` - number of electrons satisfying all veto electron requirements for resolved Higgsino analysis
* `el_*` - variables of interest for all electrons satisfying the veto id, eta and pt requirements, but no isolation 

#### Photons 

Calculated in [photon_producer](src/photon_producer.cpp):
* `nphoton` - number of signal photons
* `photon_*` - photon variables

* `nfsrphoton` - number of FSR photons
* `fsrphoton_*` - FSR photon variables

#### Tracks 

Calculated in [tk_producer](src/tk_producer.cpp):
* `ntk` - number of tracks passing criteria for resolved Higgsino analysis
* `tk_*` - track variables

####   Higgsino (HH->4b+MET) candidates
Using the 4-jet with highest DeepCSV, calculate the higgsino variables for the three possible pairings. The 0th index stores the pairing with smalled Delta m
* `hig_cand_dm` - Mass difference between the two Higgs candidates
* `hig_cand_am` - Average mass between the two Higgs candidates
* `hig_cand_drmax` - Max opening angle between the two b jets out of the two Higgs candidates.

Same variables using the 4-jet with highest DeepFlavour discriminant value are stored in:
* `hig_df_cand_*` 

* `low_dphi_*` - require dPhi(jet, MET) be less than 0.5 for jets 1,2 and less than 0.3 for jets 3,4

#### Higgs to Z gamma-related variables

Calculated in [dilep_producer](src/dilep_producer.cpp) and [zgamma_producer](src/zgamma_producer.cpp).

* `dijet_*` - variables for pairs of two highest pT signal jets
* `zg_cutBitMap` - stores whether events pass a particular version of the H to Z gamma baseline selection
* `nllphoton` - number of Higgs to Z gamma candidates
* `llphoton_*` - Higgs candidate variables

#### Higgsino (HH->bbgammagamma+MET) candidates

* `nphotonphoton` -  number of diphoton Higgs candidates
* `photonphoton_*` - diphoton Higgs candidate variables
* `nbb` -  number of b bbar Higgs candidates using deepCSV
* `bb_*` - b bbar Higgs candidate variables
* `nbb_df` -  number of b bbar Higgs candidates using deep flavor
* `bb_df_*` - b bbar Higgs candidate variables
* `nbbphotonphoton` - number of bb gammagamma Higgs candidate pairs using deepCSV
* `bbphotonphoton_*` - bb gammagamma Higgs candidate pair variables
* `nbbphotonphoton_df` - number of bb gammagamma Higgs candidate pairs using deep flavor
* `bbphotonphoton_df_*` - bb gammagamma Higgs candidate pair variables

#### Quality 

* `pass_*` - recommended MET filters
* `pass_jets` - set to false if any of the jets fails loose ID
* `pass` - logical AND of recommended and optional filters (currently not including RA2b for UL)

#### Truth 

* `mc_*` - information for a set of the generator particles in the hard process
* `ntrumu,ntruel,ntrutauh,ntrutaul` - # of true leptons of particular type, where tauh is hadronically decaying taus and taul is leptonically decaying taus
* `ntrulep = ntrumu + ntruel + ntrutaul` 
* `isr_tru_*` - MC truth, hadronic recoil, used for ISR reweighting used by the SUS PAG for weak production
* `mprod, mlsp` - higgsino and lsp mass, with lsp mass always equal to one for the higgsino model
* `npu_tru, npu_tru_mean` - true number of pileup vertices
* `mt_tru` - transverse mass at truth level, only calculated for ntrulep==1
* `met_tru, met_tru_phi` - true missing transverse momentum
* `ht_isr_me` - scalar sum of parton pT at matrix element level
* `ngenjet` - number of truth jets
* `genjet_*` - truth jet variables

####   ISR

* `nisr` - number of ISR jets according to matching to truth, used for ISR reweighting used by the SUS PAG for strong production

* `jetsys_*` - hadronic recoil, i.e. vector sum of all jets, used in V+jets ISR studies
* `jetsys_nob_*` - hadronic recoil, i.e. vector sum of all jets that are not b-tagged, used in 2L tt+jets ISR studies

#### Weights 

Previously calculated in many dedicated files, but for UL, now calculated in [event_weighter](src/event_weighter.cpp) and [trigger_weighter](src/trigger_weighter.cpp) and then re-normalized in subsequent production steps:
* `weight` - product of some of the individual weights below
* `w_lumi` - weight to be applied to get the expected yield in 1 fb-1. 
* `w_lep` - weights to correct lepton ID efficiency, product of `w_el` and `w_mu`
* `w_el` - weights to correct electron ID efficiency
* `w_mu` - weights to correct muon ID efficiency
* `w_fs_lep` - weights to correct FastSim lepton ID efficiency
* `w_photon` - weights to correct photon ID and electron veto efficiency _currently not correct_
* `w_photon_id` - weights to correct photon ID efficiency _currently not correct_
* `w_photon_csev` - weights to correct photon electron-veto efficiency
* `w_btag` - weight to correct _medium WP only_ of deepCSV b-jet ID efficiency
* `w_btag_df` - weight to correct _medium WP only_ of deepFlavor b-jet ID efficiency
* `w_bhig` - weight to correct _all WPs_ of deepCSV b-jet ID efficiency
* `w_bhig_df` - weight to correct _all WPs_ of deepFlavor b-jet ID efficiency. This should be used as the primary b-tag weight
* `w_isr` - 1., except for TTJets 2016 and signal, SUSY ISR reweighting
* `w_pu` - weight to correct pileup distribution
* `w_prefire` - weight to correct for inefficiency caused by run 2 trigger prefiring
* `w_trig` - weight to correct for trigger efficiency

#### Trigger

The final trigger menus from each year from a post on the [Trigger HN](https://hypernews.cern.ch/HyperNews/CMS/get/trigger-prim-datasets/52/1/2.html): [2016](https://docs.google.com/spreadsheets/d/1bII_92pCrgk20A9FMIIHsOsYYj3f-lLjsoRkP_ZNQW4/edit?usp=sharing), [2017](https://docs.google.com/spreadsheets/d/1SqSKATiHcZdSB-t7rMb0SRSUL4Ot1jr2TExjnMrDUpA/edit?usp=sharing), [2018](https://docs.google.com/spreadsheets/d/1D_og1_J6Hp4uALUg-R4Hkq3CF0VN6IK5komHPHv5R-o/edit?usp=sharing). Trigger menus for run 3: [2022?](https://docs.google.com/spreadsheets/d/1ThboqmcVpMGdK8uGi2PZu4zllzMwU88nCD9t_ZhO_d0/edit#gid=0), [2023](https://docs.google.com/spreadsheets/d/1yaj9Za8pcXknuO25O5Jpa565FiC425g4g6ACqI32gkU/edit#gid=1253663040).

Calcalated in [event_tools](src/event_tools.cpp).
* `trig_single_el, trig_double_el, trig_single_mu, trig_double_mu` - flag indicating events passed the lowest unprescaled single/double electron/muon trigger for the particular data taking era considered
* `HLT_*` - trigger decisions

#### Systematic uncertainties

* `sys_*` - systematic variations of weights up=0, down=1
* `sys_bchig, sys_udsghig, sys_fs_bchig, sys_fs_udsghig` - variations in (FastSim) b-tagging weights, split by heavy flavor (b/c) and light flavor (u/d/s/g) jets
* `sys_murf` - variations in renormalization and factorization scales

* `sys_jet_*` - systematic variations in jet pt/mass 
* `sys_photon_*` - systematic variations in photon pt
* `sys_el_*` - systematic variations in electron pt

For the following uncertainties the indices correspond to 0=JER up, 1=JER down, 2=JES up, 3=JES down.
* `sys_njet` - `njet` as jet systematics are varied
* `sys_nb*` - `nb*` as jet systematics are varied
* `sys_hig_cand_*` - `hig_cand_*` as jet systematics are varied
* `sys_low_dphi_met` - `low_dphi_met` as jet systematics are varied
* `sys_ht` - `ht` as jet systematics are varied

For the following uncertainties the indices correspond to 0=JER up, 1=JER down, 2=JES up, 3=JES down, 4=unclustered energy up, 5=unclustered energy down, 6=lepton/photon up, 7=lepton/photon down
* `sys_met` - `met` as jet systematics are varied
* `sys_met_phi` - `met_phi` as jet systematics are varied

## Spliting Higgsino signal

### Step 0. Setup environment

  source set_env.sh

### Step . Correct FastSim JEC on NanoAODs / Add variations on FullSIM signal

~~~~bash 
./scripts/write_fastsim_jmeCorrection_cmds.py --inputNanoAodFolder /net/cms25/cms25r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D_unsplit
                                              --outputNanoAodFolder /net/cms24/cms24r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D_unsplit_fastSimJmeCorrection
                                              --commandFilename apply_fastsim_jmeCorrection_2016.py

./scripts/write_fastsim_jmeCorrection_cmds.py --inputNanoAodFolder /net/cms25/cms25r0/pico/NanoAODv7/nano/2017/SMS-TChiHH_2D_unsplit
                                              --outputNanoAodFolder /net/cms24/cms24r0/pico/NanoAODv7/nano/2017/SMS-TChiHH_2D_unsplit_fastSimJmeCorrection
                                              --commandFilename apply_fastsim_jmeCorrection_2017.py

./scripts/write_fastsim_jmeCorrection_cmds.py --inputNanoAodFolder /net/cms25/cms25r0/pico/NanoAODv7/nano/2018/SMS-TChiHH_2D_unsplit
                                              --outputNanoAodFolder /net/cms24/cms24r0/pico/NanoAODv7/nano/2018/SMS-TChiHH_2D_unsplit_fastSimJmeCorrection
                                              --commandFilename apply_fastsim_jmeCorrection_2018.py
~~~~

~~~~bash 
./scripts/write_fastsim_jmeCorrection_cmds.py --inputNanoAodFolder /net/cms24/cms24r0/pico/NanoAODv7/nano/2016/SMS-T5qqqqZH_unsplit
                                              --outputNanoAodFolder /net/cms24/cms24r0/pico/NanoAODv7/nano/2016/SMS-T5qqqqZH_unsplit_fastSimJmeCorrection
                                              --commandFilename gluino_apply_fastsim_jmeCorrection_2016.py

./scripts/write_fastsim_jmeCorrection_cmds.py --inputNanoAodFolder /net/cms24/cms24r0/pico/NanoAODv7/nano/2017/SMS-T5qqqqZH_unsplit
                                              --outputNanoAodFolder /net/cms24/cms24r0/pico/NanoAODv7/nano/2017/SMS-T5qqqqZH_unsplit_fastSimJmeCorrection
                                              --commandFilename gluino_apply_fastsim_jmeCorrection_2017.py

./scripts/write_fastsim_jmeCorrection_cmds.py --inputNanoAodFolder /net/cms24/cms24r0/pico/NanoAODv7/nano/2018/SMS-T5qqqqZH_unsplit
                                              --outputNanoAodFolder /net/cms24/cms24r0/pico/NanoAODv7/nano/2018/SMS-T5qqqqZH_unsplit_fastSimJmeCorrection
                                              --commandFilename gluino_apply_fastsim_jmeCorrection_2018.py
~~~~

~~~~bash 
auto_submit_jobs.py apply_fastsim_jmeCorrection_2016.json -c scripts/check_jmeCorrection.py
auto_submit_jobs.py apply_fastsim_jmeCorrection_2017.json -c scripts/check_jmeCorrection.py
auto_submit_jobs.py apply_fastsim_jmeCorrection_2018.json -c scripts/check_jmeCorrection.py
~~~~

For FullSIM, 

~~~~bash
./scripts/write_fullsim_jmeCorrection_cmds.py --inputNanoAodFolder /net/cms24/cms24r0/pico/NanoAODv7/nano/2016/SMS-T5qqqqZH_FullSim
                                              --outputNanoAodFolder /net/cms24/cms24r0/pico/NanoAODv7/nano/2016/SMS-T5qqqqZH_FullSimJmeVariations
                                              --commandFilename gluino_apply_fullsim_jmeCorrection_2016.py
./scripts/write_fullsim_jmeCorrection_cmds.py --inputNanoAodFolder /net/cms24/cms24r0/pico/NanoAODv7/nano/2017/SMS-T5qqqqZH_FullSim
                                              --outputNanoAodFolder /net/cms24/cms24r0/pico/NanoAODv7/nano/2017/SMS-T5qqqqZH_FullSimJmeVariations
                                              --commandFilename gluino_apply_fullsim_jmeCorrection_2017.py
./scripts/write_fullsim_jmeCorrection_cmds.py --inputNanoAodFolder /net/cms24/cms24r0/pico/NanoAODv7/nano/2018/SMS-T5qqqqZH_FullSim
                                              --outputNanoAodFolder /net/cms24/cms24r0/pico/NanoAODv7/nano/2018/SMS-T5qqqqZH_FullSimJmeVariations
                                              --commandFilename gluino_apply_fullsim_jmeCorrection_2018.py
~~~~

~~~~bash 
auto_submit_jobs.py apply_fullsim_jmecorrection_2016.json -c scripts/check_jmeCorrection.py
auto_submit_jobs.py apply_fullsim_jmecorrection_2017.json -c scripts/check_jmeCorrection.py
auto_submit_jobs.py apply_fullsim_jmecorrection_2018.json -c scripts/check_jmeCorrection.py
~~~~

### Step 1. Split scan
Generate a python file that prints the commands to be run in the batch (input for the queue system):

~~~~bash 
./scripts/write_split_signal_mass_points_cmds.py --two_dim
                                                 --in_dir /net/cms29/cms29r0/pico/NanoAODv5/nano/2016/SMS-TChiHH_2D_unsplit 
                                                 --target_dir /net/cms29/cms29r0/pico/NanoAODv5/nano/2016/SMS-TChiHH_2D 
                                                 --dataset_filenames SMS-TChiHH_*_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv5__PUSummer16v3Fast_Nano1June2019_102X_mcRun2_asymptotic_v7-v1_*.root
                                                 --out_cmd_filename cmds_split.py
~~~~

~~~~bash 
./scripts/write_split_signal_mass_points_cmds.py --in_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D_unsplit_fastSimJmeCorrection 
                                                 --target_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D_fastSimJmeCorrection 
                                                 --dataset_filenames SMS-TChiHH_HToBB_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv7__PUSummer16v3Fast_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1*.root
                                                 --out_cmd_filename cmds_split_2016_1D.py
./scripts/write_split_signal_mass_points_cmds.py --two_dim --in_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D_unsplit_fastSimJmeCorrection 
                                                 --target_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D_fastSimJmeCorrection 
                                                 --dataset_filenames SMS-TChiHH_HToBB_HToBB_2D_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv7__PUSummer16v3Fast_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1*.root
                                                 --out_cmd_filename cmds_split_2016_2D.py
./scripts/write_split_signal_mass_points_cmds.py --in_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2017/SMS-TChiHH_2D_unsplit_fastSimJmeCorrection 
                                                 --target_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2017/SMS-TChiHH_2D_fastSimJmeCorrection 
                                                 --dataset_filenames SMS-TChiHH_HToBB_HToBB_TuneCP2_13TeV-madgraphMLM-pythia8__RunIIFall17NanoAODv7__PUFall17Fast_Nano02Apr2020_102X_mc2017_realistic_v8-v1*.root
                                                 --out_cmd_filename cmds_split_2017_1D.py
./scripts/write_split_signal_mass_points_cmds.py --two_dim --in_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2017/SMS-TChiHH_2D_unsplit_fastSimJmeCorrection 
                                                 --target_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2017/SMS-TChiHH_2D_fastSimJmeCorrection 
                                                 --dataset_filenames SMS-TChiHH_HToBB_HToBB_2D_TuneCP2_13TeV-madgraphMLM-pythia8__RunIIFall17NanoAODv7__PUFall17Fast_Nano02Apr2020_102X_mc2017_realistic_v8-v1*.root
                                                 --out_cmd_filename cmds_split_2017_2D.py
./scripts/write_split_signal_mass_points_cmds.py --in_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2018/SMS-TChiHH_2D_unsplit_fastSimJmeCorrection 
                                                 --target_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2018/SMS-TChiHH_2D_fastSimJmeCorrection 
                                                 --dataset_filenames SMS-TChiHH_HToBB_HToBB_TuneCP2_13TeV-madgraphMLM-pythia8__RunIIAutumn18NanoAODv7__PUFall18Fast_Nano02Apr2020_102X_upgrade2018_realistic_v21-v1*.root
                                                 --out_cmd_filename cmds_split_2018_1D.py
./scripts/write_split_signal_mass_points_cmds.py --two_dim --in_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2018/SMS-TChiHH_2D_unsplit_fastSimJmeCorrection 
                                                 --target_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2018/SMS-TChiHH_2D_fastSimJmeCorrection 
                                                 --dataset_filenames SMS-TChiHH_HToBB_HToBB_2D_TuneCP2_13TeV-madgraphMLM-pythia8__RunIIAutumn18NanoAODv7__PUFall18Fast_Nano02Apr2020_102X_upgrade2018_realistic_v21-v1*.root
                                                 --out_cmd_filename cmds_split_2018_2D.py

convert_cl_to_jobs_info.py cmds_split.py cmds_split.py.json
auto_submit_jobs.py cmds_split.py.json -c jobscript_check.py -n cms1
~~~~


model in the below commands are use for globbing files in the input directory.

~~~~bash 
./scripts/write_split_gluino_mass_points_cmds.py --in_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2016/SMS-T5qqqqZH_unsplit_fastSimJmeCorrection 
                                                 --out_cmd_filename split_gluino_2016.py 
                                                 --target_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2016/SMS-T5qqqqZH_fastSimJmeCorrection 
                                                 --model "SMS-T5qqqqZH_HToBB-mGluino"
./scripts/write_split_gluino_mass_points_cmds.py --in_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2016/SMS-T5qqqqZH_unsplit_fastSimJmeCorrection 
                                                 --out_cmd_filename split_gluino_mN2_2016.py 
                                                 --target_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2016/SMS-T5qqqqZH_fastSimJmeCorrection 
                                                 --model "SMS-T5qqqqZH_HToBB-mN2"

convert_cl_to_jobs_info.py split_gluino_2016.py split_gluino_2016.json
auto_submit_jobs.py split_gluino_2016.json -c jobscript_check.py -n cms1
convert_cl_to_jobs_info.py split_gluino_mN2_2016.py split_gluino_mN2_2016.json
auto_submit_jobs.py split_gluino_mN2_2016.json -c jobscript_check.py -n cms1

./scripts/write_split_gluino_mass_points_cmds.py --in_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2017/SMS-T5qqqqZH_unsplit_fastSimJmeCorrection 
                                                 --out_cmd_filename split_gluino_2017.py 
                                                 --target_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2017/SMS-T5qqqqZH_fastSimJmeCorrection 
                                                 --model "SMS-T5qqqqZH_HToBB-mGluino"
./scripts/write_split_gluino_mass_points_cmds.py --in_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2017/SMS-T5qqqqZH_unsplit_fastSimJmeCorrection 
                                                 --out_cmd_filename split_gluino_mN2_2017.py 
                                                 --target_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2017/SMS-T5qqqqZH_fastSimJmeCorrection 
                                                 --model "SMS-T5qqqqZH_HToBB-mN2"

./scripts/write_split_gluino_mass_points_cmds.py --in_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2018/SMS-T5qqqqZH_unsplit_fastSimJmeCorrection 
                                                 --out_cmd_filename split_gluino_2018.py 
                                                 --target_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2018/SMS-T5qqqqZH_fastSimJmeCorrection 
                                                 --model "SMS-T5qqqqZH_HToBB-mGluino"
./scripts/write_split_gluino_mass_points_cmds.py --in_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2018/SMS-T5qqqqZH_unsplit_fastSimJmeCorrection 
                                                 --out_cmd_filename split_gluino_mN2_2018.py 
                                                 --target_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2018/SMS-T5qqqqZH_fastSimJmeCorrection 
                                                 --model "SMS-T5qqqqZH_HToBB-mN2"
~~~~

~~~~bash 
./scripts/write_split_signal_mass_points_cmds.py --in_dir /net/cms25/cms25r5/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D_unsplit 
                                                 --target_dir /net/cms25/cms25r5/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D 
                                                 --dataset_filenames SMS-TChiHH_HToBB_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv7__PUSummer16v3Fast_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1*.root
                                                 --out_cmd_filename cmds_split_2016_1D.py
./scripts/write_split_signal_mass_points_cmds.py --two_dim --in_dir /net/cms25/cms25r5/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D_unsplit 
                                                 --target_dir /net/cms25/cms25r5/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D 
                                                 --dataset_filenames SMS-TChiHH_HToBB_HToBB_2D_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv7__PUSummer16v3Fast_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1*.root
                                                 --out_cmd_filename cmds_split_2016_2D.py
./scripts/write_split_signal_mass_points_cmds.py --in_dir /net/cms25/cms25r5/pico/NanoAODv7/nano/2017/SMS-TChiHH_2D_unsplit 
                                                 --target_dir /net/cms25/cms25r5/pico/NanoAODv7/nano/2017/SMS-TChiHH_2D 
                                                 --dataset_filenames SMS-TChiHH_HToBB_HToBB_TuneCP2_13TeV-madgraphMLM-pythia8__RunIIFall17NanoAODv7__PUFall17Fast_Nano02Apr2020_102X_mc2017_realistic_v8-v1*.root
                                                 --out_cmd_filename cmds_split_2017_1D.py
./scripts/write_split_signal_mass_points_cmds.py --two_dim --in_dir /net/cms25/cms25r5/pico/NanoAODv7/nano/2017/SMS-TChiHH_2D_unsplit 
                                                 --target_dir /net/cms25/cms25r5/pico/NanoAODv7/nano/2017/SMS-TChiHH_2D 
                                                 --dataset_filenames SMS-TChiHH_HToBB_HToBB_2D_TuneCP2_13TeV-madgraphMLM-pythia8__RunIIFall17NanoAODv7__PUFall17Fast_Nano02Apr2020_102X_mc2017_realistic_v8-v1*.root
                                                 --out_cmd_filename cmds_split_2017_2D.py
./scripts/write_split_signal_mass_points_cmds.py --in_dir /net/cms25/cms25r5/pico/NanoAODv7/nano/2018/SMS-TChiHH_2D_unsplit 
                                                 --target_dir /net/cms25/cms25r5/pico/NanoAODv7/nano/2018/SMS-TChiHH_2D 
                                                 --dataset_filenames SMS-TChiHH_HToBB_HToBB_TuneCP2_13TeV-madgraphMLM-pythia8__RunIIAutumn18NanoAODv7__PUFall18Fast_Nano02Apr2020_102X_upgrade2018_realistic_v21-v1*.root
                                                 --out_cmd_filename cmds_split_2018_1D.py
./scripts/write_split_signal_mass_points_cmds.py --two_dim --in_dir /net/cms25/cms25r5/pico/NanoAODv7/nano/2018/SMS-TChiHH_2D_unsplit 
                                                 --target_dir /net/cms25/cms25r5/pico/NanoAODv7/nano/2018/SMS-TChiHH_2D 
                                                 --dataset_filenames SMS-TChiHH_HToBB_HToBB_2D_TuneCP2_13TeV-madgraphMLM-pythia8__RunIIAutumn18NanoAODv7__PUFall18Fast_Nano02Apr2020_102X_upgrade2018_realistic_v21-v1*.root
                                                 --out_cmd_filename cmds_split_2018_2D.py

convert_cl_to_jobs_info.py cmds_split.py cmds_split.py.json
auto_submit_jobs.py cmds_split.py.json -c jobscript_check.py

or 
[cms25] ./scripts/run_commands.py cmds_split_2018_2D.py 
Should check log file for error, segmentation..
~~~~

Dataset name needs to be writen out like above.
This produces the commands in `cmds_split.py`. You can perform a last check by running one of the commands interactively. Next, submit the jobs to the batch system. Note the -c option which allows to attach a script that compares the input and output number of entries when each job is done. Note the check can be performed later if one needs to detach the session. Alternatively, this command can be started in screen:

~~~~bash 
convert_cl_to_jobs_info.py cmds_split.py cmds_split.py.json
auto_submit_jobs.py cmds_split.py.json -c scripts/check_skim.py
~~~~

This command will result in `checked_auto_cmds_split.py.json`, which can then be used to resubmit failed jobs if any:

~~~~bash 
select_resubmit_jobs.py checked_auto_cmds_split.py.json -c scripts/check_skim.py
auto_submit_jobs.py resubmit_checked_auto_cmds_split.py.json -c scripts/check_skim.py
~~~~

~~~~bash
./scripts/write_split_gluino_mass_points_cmds.py --in_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2016/SMS-T5qqqqZH_FullSim_unsplit 
                                                 --out_cmd_filename split_fullsim_gluino_2016.py 
                                                 --target_dir /net/cms24/cms24r0/pico/NanoAODv7/nano/2016/SMS-T5qqqqZH_FullSim 
                                                 --model "SMS-T5qqqqZH-mGluino"
~~~~

## Getting Higgsino cross-section

### Step 1. Download cross-section files

~~~~bash
scp -r lxplus:/afs/cern.ch/user/a/amete/public/EWKGauginoCrossSections_13TeV cross_section
~~~~

### Step 2. Fix bug in script

Set masses, xsecs, xsecUncs to 0 when initializing.

~~~~bash
cd cross_section
sed -i 's/std::vector<double>\* masses;/std::vector<double>\* masses=0;/' get_gaugino.C
sed -i 's/std::vector<double>\* xsecs;/std::vector<double>\* xsecs=0;/' get_gaugino.C
sed -i 's/std::vector<double>\* xsecUncs;/std::vector<double>\* xsecUncs=0;/' get_gaugino.C
~~~~

## Getting Gluino cross-section

wget https://raw.githubusercontent.com/fuenfundachtzig/xsec/master/json/pp13_gluino_NNLO%2BNNLL.json
scripts/get_gluino_cross_sections.py

### Step 3. Run script for all mass points.

Use model "CN" (mixing) or "N1N2" (no mixing).

~~~~bash
cd cross_section
../scripts/get_higgsino_cross_sections.py -i /net/cms29/cms29r0/pico/NanoAODv5/nano/2016/SMS-TChiHH_2D --model "N1N2"
~~~~

# Tagging code with git

Confirm tags: `git tag`

Add a lightweight tag: `git tag <tagname>`

Pushing tag: `git push origin <tagname>`

Deleting tag: `git tag -d <tagname>` and `git push origin --delete <tagname>`

Checking out tag: `git checkout <tagname>`

# Update submodule

In case didn't checkout submodule: `git submodule init` and `git submodule update`

Update all submodules: `git submodule update --recursive --remote --merge` and then commit.

## Validation

Validation is done by running `./script/produce_unit_test*` on nano2pico code versions and then using `./script/validate_unit_test*` to compare between results.

Below are examples

~~~bash
# Compares picos file between old code and new code. Also compares production time.
[In old code folder] ./scripts/produce_unit_test_htozgamma_NanoAODv9.py --output_folder unit_test_htozgamma_nanoaodv9 --output_log unit_test_htozgamma_nanoaodv9.log
[In new code folder] ./scripts/produce_unit_test_htozgamma_NanoAODv9.py --output_folder unit_test_htozgamma_nanoaodv9 --output_log unit_test_htozgamma_nanoaodv9.log
./scripts/validate_unit_test_picos.py --output_log_filename validate_unit_test_htozgamma_nanoaodv9.log --unit_test_log_filename unit_test_htozgamma_nanoaodv9.log --golden_base_folder OLD_CODE/unit_test_htozgamma_nanoaodv9 --validate_base_folder NEW_CODE/unit_test_htozgamma_nanoaodv9

# Compares cross section between old code and new code.
[In old code folder] ./scripts/produce_unit_test_cross_sections.py --output_log unit_test_cross_section.log
[In new code folder] ./scripts/produce_unit_test_cross_sections.py --output_log unit_test_cross_section.log
./scripts/validate_unit_test_cross_section.py --output_filename validate_unit_test_cross_section.log --golden_cross_section_log OLD_CODE/unit_test_cross_section.log --validate_cross_section_log NEW_CODE/unit_test_cross_section.log
~~~
