#!/bin/env python
source_directory = "/net/cms11/cms11r0/pico/NanoAODv15/nano/2024/SMS-TChiZH-Hto2G-FullSIM-Central-Run3/"
target_directory = "/net/cms11/cms11r0/pico/NanoAODv15/nano/2024/SMS-TChiZH-Hto2G-FullSIM-Central-Run3_split/"
mass_points = [(150, 0), (200, 0), (200, 50), (250, 0), (250, 50), (250, 100), (300, 0), (300, 50), (300, 100), (300, 150), (350, 0), (350, 50), (350, 100), (350, 150), (350, 200), (400, 0), (400, 50), (400, 100), (400, 150), (400, 200), (400, 250), (450, 0), (450, 50), (450, 100), (450, 150), (450, 200), (450, 250), (450, 300), (500, 0), (500, 50), (500, 100), (500, 150), (500, 200), (500, 250), (500, 300), (500, 350)]
for mass_point in mass_points:
  print("/net/cms37/data1/mhussain/HH-MET/nano2pico_apr2025/nano2pico/scripts/skim_file.py -m "+str(mass_point[0])+" -l "+str(mass_point[1])+" -i '"+source_directory+"SMS-TChiZH-Hto2G_*_TuneCP5_13p6TeV_madgraphMLM-pythia8__RunIII2024Summer24NanoAODv15__150X_mcRun3_2024_realistic_v2-v2*.root' -o "+target_directory)
