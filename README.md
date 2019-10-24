# nano2pico

Utility package for converting NanoAOD to "pico" analysis-ready ntuples.

### Interactive test

Step 1. Produce raw pico ntuple from a nano input file, adding `--isFastsim` and `--isData` if applicable:

~~~~bash
./compile.sh && ./run/process_nano.exe --in_file INFILE --wgt_sums_file WGT_SUMS_FILE --out_file OUTFILE
~~~~

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
./compile.sh && ./run/apply_corr.exe --in_file PICO_STEP1 --corr_file CORR_STEP2 --out_file OUTFILE
~~~~

### Calculating b-tagging efficiencies

Use `parameterize_efficiency.cxx`, giving the directory with all the MC files and the year as arguments. Below is an example run for 2016 MC.

~~~~bash
./compile.sh && ./run/parameterize_efficiency.exe -i /mnt/hadoop/jbkim/2019_09_30/2016/mc/ -y 2016
~~~~
