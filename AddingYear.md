To add an year for data do below steps.
- Add file in txt/datasets 
  - Filename example: txt/datasets/NanoAODv11_htozgamma_2022_data_dataset_paths
  - Content is dataset name. Example: /Muon/Run2022C-ReRecoNanoAODv11-v1/NANOAOD
- Add goldenJson for data in src/in_json.cpp and src/process_nano.cxx
- Add btag_wpts for new year in src/process_nano.cxx
- Add year in src/event_weighter.cxx
- Add year in src/jet_producer.cpp
- Add year in src/photon_producer.cpp
- Add year in src/el_producer.cpp

Example commands for data
source set_env.sh
scons
export OUTDIR=out/NanoAODv11p9/htozgamma_joshuatree_v1/2023/data
export INDIR=/net/cms11/cms11r0/pico/NanoAODv11p9/nano/2023/data/
export INFILE=Muon0__Run2023C__PromptNanoAODv11p9_v1-v1__70000__22570a17-17a0-4f0e-842d-964e8a4c9f95.root
mkdir -p $OUTDIR/wgt_sums
mkdir -p $OUTDIR/raw_pico
./run/process_nano.exe -f $INFILE -i $INDIR -o $OUTDIR --nent 10000
mkdir -p $OUTDIR/skim_llg
./scripts/skim_file.py -k llg -i $OUTDIR/raw_pico/raw_pico_${INFILE} -o $OUTDIR/skim_llg/
./scripts/write_slim_and_merge_cmds.py -f --in_dir $OUTDIR/skim_llg/ --slim_name zgdata --tag ${TAG}_2023_data

Example commands for mc
source set_env.sh
scons
export OUTDIR=out/NanoAODv9/higgsino_apple_v0/2016/mc
export INDIR=/net/cms11/cms11r0/pico/NanoAODv9/nano/2016/mc/
export INFILE=TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8__RunIISummer20UL16NanoAODv9__106X_mcRun2_asymptotic_v17-v1__120000__40245E81-67EF-5944-B335-2BB3F5F348B3.root
mkdir -p $OUTDIR/wgt_sums
mkdir -p $OUTDIR/raw_pico
./run/process_nano.exe -f $INFILE -i $INDIR -o $OUTDIR --nent 10000
mkdir -p $OUTDIR/skim_llg
./scripts/skim_file.py -k llg -i $OUTDIR/raw_pico/raw_pico_${INFILE} -o $OUTDIR/skim_llg/
./scripts/write_slim_and_merge_cmds.py -f --in_dir $OUTDIR/skim_llg/ --slim_name zgdata --tag ${TAG}_2023_data

Production commands
export TAG=htozgamma_joshuatree_v1
./scripts/write_process_nano_cmds.py --in_dir $INDIR --production $TAG --dataset_list txt/datasets/NanoAODv11p9_htozgamma_2023_data_dataset_paths --data --tag ${TAG}_2023_data
./scripts/write_skim_cmds.py --in_dir $OUTDIR/raw_pico/ --skim_name llg --tag ${TAG}_2023_data
./scripts/write_slim_and_merge_cmds.py -f --in_dir $OUTDIR/skim_llg/ --slim_name zgdata --tag ${TAG}_2023_data
./scripts/confirm_slim.py $OUTDIR/skim_llg $OUTDIR/merged_zgdata_llg
