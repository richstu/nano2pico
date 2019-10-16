#! /usr/bin/env python

#!/usr/bin/env python

import numpy as np
import uproot
import uproot_methods
import awkward

# histogram creation and manipulation
from coffea import hist

import matplotlib.pyplot as plt

nent = -1

fpico = uproot.open("pico1k.root")
pico_tree = fpico["tree"]
len_pico = len(pico_tree)

fbaby = uproot.open("baby1k.root")
baby_tree = fbaby["tree"]
len_baby = len(baby_tree)

if nent<0: nent = len_pico if len_pico<len_baby else len_baby

filters = [
  'pass_jets',
  'pass_hbhe',
  'pass_hbheiso',
  'pass_goodv',
  'pass_cschalo_tight',
  'pass_eebadsc',
  'pass_ecaldeadcell',
  'pass_fsjets',
  'pass_badpfmu',
  # 'pass_ra2_badmu',
  'pass_badcalib',
  'pass'
]

template = '{:>20} {:>15} {:>15}'
cols = template.format('Variable', 'Pico', 'Baby') + '\n'
cols += template.format('# events', len_pico, len_baby)
print(cols)

event = baby_tree.array(b'event',entrystop=nent)
for filt in filters:
  px = pico_tree.array(filt.encode(),entrystop=nent)
  bx = baby_tree.array(filt.encode(),entrystop=nent)
  cols = template.format(filt, np.count_nonzero(px), np.count_nonzero(bx))
  print(cols)

  if (filt=='pass_jets'):
    print("Events with !pass_jets:")
    b_bad_events = np.where(bx==False)
    p_bad_events = np.where(px==False)

diff = np.setdiff1d(b_bad_events, p_bad_events)

baby_jets = baby_tree.arrays([b'jets_*'])
pico_jets = pico_tree.arrays([b'jet_*'])
for idx in diff:
  # print out the jets pt, eta, id decision
  baby_njets = len(baby_jets[b'jets_pt'][idx])
  pico_njets = len(pico_jets[b'jet_pt'][idx])
  if (baby_njets != pico_njets): print(' Different number of jets saved!')

  for ijet in range(baby_njets):
    print(''.format(baby_jets[b'jets_pt'][idx][ijet], baby_jets[b'jets_pt'][idx][ijet], baby_jets[b'jets_isgood'][idx][ijet]))
  for ijet in range(pico_njets):
    print(''.format(pico_jets[b'jet_pt'][idx][ijet], pico_jets[b'jet_pt'][idx][ijet], pico_jets[b'jet_id'][idx][ijet]))











