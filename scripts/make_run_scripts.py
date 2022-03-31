#!/usr/bin/env python3
import sys
import os

# source_target = run/path_to_file.exe
source_target = sys.argv[1]

out_string = '''#!/bin/sh
RUN_KERNEL=$(uname -r | cut -d '-' -f1)
kernel/$RUN_KERNEL/'''+source_target+''' "$@"
'''

with open(source_target, 'w') as source_script:
  source_script.write(out_string)
  #print('Made '+source_script_path)
  os.chmod(source_target, 0o755)
