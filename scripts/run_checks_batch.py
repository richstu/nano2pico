#!/usr/bin/env python3
"""Run many job-check commands in one process, so ROOT is imported once instead of once per job.

Used automatically by queue_system: when this file sits next to the job check script (as in
scripts/), the checks of finished jobs are run through it in groups (batch_check_size in
queue_system.py) instead of one process per job. Delete or rename it to go back to one-by-one checks.

Usage: run_checks_batch.py <commands.json> <results.json>
  commands.json: list of check commands, exactly the strings queue_system would otherwise run one
                 by one, e.g. "scripts/check_skim.py <compressed log> <compressed arguments>"
  results.json:  list with the stdout of each command, in the same order.

Each command is run with runpy as if it were started from the shell (same argv, SystemExit
handled), and pico_count is switched to count entries with PyROOT in this process rather than
starting count_entries.exe for every job. A command that raises or exits non-zero gets the
same "check script crashed" result queue_system gives to a failing command.
"""
import contextlib
import io
import json
import os
import runpy
import shlex
import sys


def run_command(command):
  tokens = shlex.split(command)
  script = tokens[0]
  script_dir = os.path.dirname(os.path.realpath(script))
  crashed = '[For queue_system] fail: check script crashed: '+command
  old_argv = sys.argv
  sys.argv = tokens
  sys.path.insert(0, script_dir) # as when the script is run directly: lets it import its neighbours
  output = io.StringIO()
  try:
    with contextlib.redirect_stdout(output):
      runpy.run_path(script, run_name='__main__')
  except SystemExit as exit_error:
    if exit_error.code not in (None, 0):
      return crashed
  except Exception as error:
    print('[run_checks_batch] '+script+' raised '+repr(error), file=sys.stderr)
    return crashed
  finally:
    sys.argv = old_argv
    sys.path.remove(script_dir)
  return output.getvalue()


if __name__ == '__main__':
  sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
  import pico_count
  pico_count.INPROCESS = True
  with open(sys.argv[1]) as commands_file:
    commands = json.load(commands_file)
  results = [run_command(command) for command in commands]
  with open(sys.argv[2], 'w') as results_file:
    json.dump(results, results_file)

