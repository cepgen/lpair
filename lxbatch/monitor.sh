#!/bin/sh
export BASE=`pwd`
export SCRAM_ARCH="slc5_amd64_gcc462"
#cd /afs/cern.ch/user/l/lforthom/scratch0/FastSimHLT/CMSSW_6_0_0_pre2
cd /afs/cern.ch/user/l/lforthom/scratch0/FastSimHLT/CMSSW_5_3_5
source /afs/cern.ch/cms/cmsset_default.sh
eval `scramv1 runtime -sh`
cd $BASE
FAILEDJOBS=`python /afs/cern.ch/user/l/lforthom/lxbatch/jobsStats.py $1 1`
#echo $FAILEDJOBS
if [[ ! -z $FAILEDJOBS ]] ; then
    echo "Failed jobs :"
    RESUB=`python /afs/cern.ch/user/l/lforthom/lxbatch/resubmitJobs.py $1 $FAILEDJOBS`
    echo $RESUB
else
    echo "No failed jobs ; exiting..."
fi
