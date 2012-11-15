#!/bin/sh
export BASE=`pwd`
export SCRAM_ARCH="slc5_amd64_gcc462"
#cd /afs/cern.ch/user/l/lforthom/scratch0/FastSimHLT/CMSSW_6_0_0_pre2
cd /afs/cern.ch/user/l/lforthom/scratch0/FastSimHLT/CMSSW_5_3_5
source /afs/cern.ch/cms/cmsset_default.sh
eval `scramv1 runtime -sh`
cd $BASE
python jobsStats.py $1