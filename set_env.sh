#!/bin/bash
. /cvmfs/cms.cern.ch/cmsset_default.sh

RUN_KERNEL=$(uname -r | cut -d '-' -f1)
if [ "$RUN_KERNEL" == "3.10.0" ]; then
  export SCRAM_ARCH=slc7_amd64_gcc700
  cd /net/cms29/cms29r0/pico/cc7/CMSSW_10_2_11_patch1/src
elif [ "$RUN_KERNEL" == "2.6.32" ]; then
  cd /net/cms29/cms29r0/pico/CMSSW_10_2_11_patch1/src
fi

eval `scramv1 runtime -sh`
cd -

export SCONSFLAGS="-j $(nproc --all)"

source $(dirname $(readlink -e "$BASH_SOURCE"))/modules/jb_utils/set_env.sh
source $(dirname $(readlink -e "$BASH_SOURCE"))/modules/queue_system/set_env.sh

# Setup batch
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
  cat /net/cms2/cms2r0/jbkim/log/$1.log
}
