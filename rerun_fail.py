#!/usr/bin/env python3
import json
import subprocess
import sys

json_file = "checked_auto_process_nano_cmds_higgsino_run3_blanc_v6_2025_data.json"

with open(json_file) as f:
    jobs = json.load(f)

failed = [j for j in jobs if j.get("job_status") != "success"]

print(f"Found {len(failed)} failed/non-success jobs out of {len(jobs)} total\n")

for j in failed:
    print(f"  Status: {j.get('job_status')}")
    print(f"  Command: {j.get('command')}")
    print()

# dry_run = False to actually execute
dry_run = False

for j in failed:
    cmd = j.get("command")
    if not cmd:
        continue
    print(f"Running: {cmd}")
    if not dry_run:
        ret = subprocess.run(cmd, shell=True)
        print(f"  -> exit code {ret.returncode}")

