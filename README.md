# nano2pico

Utility package for converting NanoAOD to "pico" analysis-ready ntuples.

### Interactive test

Step 1. Produce raw pico ntuple from a nano input file, adding `--isFastsim` and `--isData` if applicable:

~~~~bash
./compile.sh && ./run/process_nano.exe --in_file INFILE --in_dir INDIR --out_dir OUTDIR
~~~~

Note that for interactive jobs, you need to ensure that the output directory subfolders `wgt_sums` and `raw_pico` exist.

:bangbang: Code functionality relies on the input NanoAOD filename! Specifically, `INFILE` is parsed for:

* flag `isData = infile.Contains("Run201") ? true : false;`
* flag `isFastsim = infile.Contains("Fast") ? true : false;`
* variable `year = infile.Contains("RunIISummer16") ? 2016 : (infile.Contains("RunIIFall17") ? 2017 : 2018)`
* output branch `type` is set based on the presence of dataset name substrings (see event_tools.cpp)
* branches related on ISR also depend on the presence of dataset name substrings

Step 2. For each dataset, add up the sums of weights obtained for each file in step 1 and calculate the corrections needed to normalize each individual weight as well as the total weight. Output to `CORR_FILE`. Note that the order of options is fixed with the last argument being the input files in order to allow arbitrary number of input files.

~~~~bash
./compile.sh && ./run/merge_corrections.exe YYYY CORR_FILE WGT_SUMS_FILE1 WGT_SUMS_FILE2 ...
~~~~

Step 3. Using the pico file from step 1 and the corrections file from step 2 as input, we can renormalize the weight branches as follows:

~~~~bash
./compile.sh && ./run/apply_corr.exe --in_file PICO_STEP1 --in_dir INDIR --corr_file CORR_STEP2
~~~~

### Batch system

Setup Jaebak's queue system:

~~~~bash
CMSSW=/net/top/homes/jbkim/analysis/CMSSW
RELEASE=CMSSW_10_2_11_patch1

. /cvmfs/cms.cern.ch/cmsset_default.sh;
cd $CMSSW/$RELEASE/src;
eval `scramv1 runtime -sh`;
cd -;

cd ../
git clone --recurse-submodules https://github.com/richstu/copydataset
cd copydatasets
source set_env.sh
cd ../nano2pico
~~~~

Step 1. Converting Nano to Pico:
Generate a python file that prints the commands to be run in the batch (input for the queue system):

~~~~bash 
./scripts/write_cmd_file.py --in_dir /mnt/hadoop/pico/NanoAODv5/nano/2016/mc/ --production higgsino_angeles --datasets_file datasets/higgsino_2016_mc_dataset_paths.txt --out_dir --out_cmd_file cmds.py
~~~~

Submit the jobs with a script that compares the input and output number of entries. Note the check can be performed later if one needs to detach the session. Or this command can be started in screen:

~~~~bash 
convert_cl_to_jobs_info.py cmds.py higgsino_angeles.json
auto_submit_jobs.py higgsino_angeles.json -c scripts/check_process_nano_job.py
~~~~

To check the jobs later on:

~~~~bash 
check_jobs.py auto_higgsino_angeles.json -c scripts/check_process_nano_job.py
~~~~

The json file here is the most recent trusted record of the status of the jobs. To resubmit failed jobs:

~~~~bash 
select_resubmit_jobs.py checked_auto_higgsino_angeles.json -c scripts/check_process_nano_job.py 
auto_submit_jobs.py resubmit_checked_auto_higgsino_angeles.json -c scripts/check_process_nano_job.py 
~~~~

### Calculating b-tagging efficiencies

Use `parameterize_efficiency.cxx`, giving the directory with all the MC files and the year as arguments. Below is an example run for 2016 MC.

~~~~bash
./compile.sh && ./run/parameterize_efficiency.exe -i /mnt/hadoop/jbkim/2019_09_30/2016/mc/ -y 2016
~~~~
