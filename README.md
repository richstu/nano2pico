# nano2pico

Utility package for converting NanoAOD to "pico" analysis-ready ntuples.

### Interactive test

Example given for a file from 2016 MC.

Step 1. Write out raw pico ntuple, adding `--isFastsim` and `--isData` if applicable:

    ./compile.sh && ./run/calc_vars.exe \
        --nano_file /Users/ana/code/sandbox/data/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8_Nano1June2019_102X_upgrade2018_realistic_v19-v1.root \
        --wgt_sums_file out/wgt_sums/wgt_sums_TTJets_TuneCP5_13TeV-madgraphMLM-pythia8_Nano1June2019_102X_upgrade2018_realistic_v19-v1.root \
        --pico_file out/raw_pico/pico_TTJets_TuneCP5_13TeV-madgraphMLM-pythia8_Nano1June2019_102X_upgrade2018_realistic_v19-v1.root
