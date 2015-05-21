#!/bin/bash

#USER INPUTS
CMS3Tag=CMS3_V07-04-01
CMSSW_release=7_4_1_patch1
export SCRAM_ARCH=slc6_amd64_gcc491

#--Here there be dragons----
export CMS_PATH=/cvmfs/cms.cern.ch
source /cvmfs/cms.cern.ch/cmsset_default.sh
scramv1 p -n CMSSW_$CMSSW_release CMSSW CMSSW_$CMSSW_release
cd CMSSW_$CMSSW_release/src
eval `scramv1 runtime -sh`
git clone git@github.com:cmstas/NtupleMaker.git CMS3/NtupleMaker
cd CMS3/NtupleMaker
git checkout $CMS3Tag
source setup/patchesToSource.sh
cd $CMSSW_BASE/src
scram b -j 10
cd ..
