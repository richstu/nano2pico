"""Entry counting for the scripts/check_*.py job checkers.
Credit: Claude
count_entries() uses count_entries.exe and only falls back to PyROOT if the executable is
missing or fails.
"""
import os
import subprocess
import sys

_TOP = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')
_KERNEL = os.uname()[2].split('-')[0]
_CANDIDATES = [os.path.join(_TOP, 'kernel', _KERNEL, 'run', 'count_entries.exe'),
               os.path.join(_TOP, 'run', 'count_entries.exe')]
EXE = next((path for path in _CANDIDATES if os.access(path, os.X_OK)), None)
INPROCESS = False

def _count_exe(path, tree, cut):
  cmd = [EXE, '-t', tree] + (['-c', cut] if cut else []) + [path]
  result = subprocess.run(cmd, stdout=subprocess.PIPE, universal_newlines=True)
  if result.returncode != 0:
    raise RuntimeError('exit code '+str(result.returncode))
  for line in result.stdout.splitlines():
    if line.startswith('COUNT '):
      fields = line.split(' ', 3)  # COUNT <entries> <branches> <path>
      return int(fields[1]), int(fields[2])
  raise RuntimeError('no COUNT line in output')


def _count_pyroot(path, tree, cut):
  from ROOT import TChain
  chain = TChain(tree)
  chain.Add(path)
  nent = chain.GetEntries(cut) if cut else chain.GetEntries()
  return nent, chain.GetNbranches()


def count_entries(path, tree='tree', cut=''):
  """Return (number of entries passing cut, number of branches) of TChain(tree) over path."""
  # run_checks_batch.py sets this: ROOT is already imported in that process, so count here
  if INPROCESS:
    return _count_pyroot(path, tree, cut)
  if EXE is not None:
    try:
      return _count_exe(path, tree, cut)
    except (OSError, RuntimeError, ValueError) as err:
      print('[pico_count] count_entries.exe failed ('+str(err)+'), falling back to PyROOT', file=sys.stderr)
  return _count_pyroot(path, tree, cut)

