export INDIR=/net/cms17/cms17r0/pico/NanoAODv9/nano/2017/signal/
export INFILE=GluGluHToZG_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL17NanoAODv9__1000.root
./compile.sh && ./run/process_nano.exe --in_file $INFILE --in_dir $INDIR --out_dir out/
