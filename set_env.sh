#!/bin/bash
# Sets root, python, scons, ssh, and batch enviornment for ucsb servers

# Setup root environment according to kernel of ucsb server
. /cvmfs/cms.cern.ch/cmsset_default.sh
RUN_KERNEL=$(uname -r | cut -d '-' -f1)
if [ "$RUN_KERNEL" == "3.10.0" ]; then
  export SCRAM_ARCH=slc7_amd64_gcc12
  cd /net/cms11/data/pico/cc7/CMSSW_14_2_2/src
elif [ "$RUN_KERNEL" == "2.6.32" ]; then
  cd /net/cms29/cms29r0/pico/CMSSW_10_2_11_patch1/src
fi
eval `scramv1 runtime -sh`
cd -

# path to scons for python 3.9
export PATH=./lib/python3.9/site-packages/bin:$PATH
export PYTHON3PATH=./lib/python3.9/site-packages/:$PYTHON3PATH

# Do multi-core scons
export SCONSFLAGS="-j $(nproc --all)"
export SET_ENV_PATH=set_env.sh # environment to use for build

# Change python to be in unbuffer mode for scripts to run commands
export PYTHONUNBUFFERED=1

# Prevent asking password using x11-gui
unset SSH_ASKPASS

# Setup environment for batch submission
source $(dirname $(readlink -e "$BASH_SOURCE"))/modules/jb_utils/set_env.sh
source $(dirname $(readlink -e "$BASH_SOURCE"))/modules/queue_system/set_env.sh
export JOBBIN=/net/cms2/cms2r0/Job
export JOBS=/net/cms2/cms2r0/${USER}/jobs
export LOG=/net/cms2/cms2r0/${USER}/log
export PATH=$JOBBIN${PATH:+:${PATH}}
alias bsub='JobSubmit.csh'
alias bjobs='JobShow.csh'
alias bkill='JobKill.csh'
unset -f bkillall 
function bkillall {
  cat $JOBS/running.list | awk '{print $1}' | xargs -I {} JobKill.csh {}
  cat $JOBS/queued.list | awk '{print $1}' | xargs -I {} JobKill.csh {}
  cat $JOBS/ready.list | awk '{print $1}' | xargs -I {} JobKill.csh {}
}
unset -f blog
function blog {
  cat logs/out.$1
}
alias bproc="tail -f /net/cms2/cms2r0/${USER}/log/JobProc.log"
