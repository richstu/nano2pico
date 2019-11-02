# nano2pico

Utility package for converting CMS NanoAOD to analysis-ready ntuples called "pico".

## Environment setting

    git clone --recurse-submodules https://github.com/richstu/nano2pico
    git submodule update --init --recursive
    source set_env.sh

Also setup CMSSW and the UCSB job environment variables (JOBBIN, JOBS, LOG, PATH).

## Intro

Variables stored in the pico can be seen in [variables/pico](variables/pico). For an overview of the available branches, see the dedicated section at the bottom of this README.

The Nano -> pico conversion is done in three steps in order to allow parallelizing the production at the sub-dataset level:
  1. All variables and event weights (except normalization) are calculated using [src/process_nano.cxx](src/process_nano.cxx). This step also keeps a tally of the weights for all events in the file being run over as input to the next step.
  2. The sums of weights from step 1 are further aggregate to get the total per dataset. A correction is then calculated to ensure that the weights do not change the total expected number of events for the dataset. This is done in [src/merge_corrections.cxx](src/merge_corrections.cxx). The luminosity normalization weight `w_lumi` to be applied to get the yield in 1fb-1 is also calculated for each dataset in this step.
  3. The `raw_pico` files from step 1 are corrected by the per-dataset correction factors derived in step 2 and written to the `unskimmed` folder.

At this point, various skims can be made as defined in [scripts/skim_file.py](scripts/skim_file.py).

The current Higgsino production nicknamed "Angeles" can be found at:
  /net/cms2/cms2r0/pico/NanoAODv5/higgsino_angeles/2016/mc/
  /net/cms2/cms2r0/pico/NanoAODv5/higgsino_angeles/2016/TChiHH/

To see the sizes and number of files, one can do:
  ./scripts/count_root_files.py -f /net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/

## Interactive test

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

## Batch system

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
                               --corr_dir /net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/corrections/ \
                               --year 2016
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
