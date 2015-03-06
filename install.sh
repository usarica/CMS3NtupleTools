#!/bin/bash

#USER INPUTS
CMS3Tag=master

export CMS_PATH=/cvmfs/cms.cern.ch
export SCRAM_ARCH=slc6_amd64_gcc481
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/root/5.34.18/bin/thisroot.sh
scramv1 p -n CMSSW_7_2_0 CMSSW CMSSW_7_2_0
cd CMSSW_7_2_0/src
eval `scramv1 runtime -sh`
git clone git@github.com:cmstas/NtupleMaker.git CMS3/NtupleMaker
cd CMS3/NtupleMaker
git checkout $CMS3Tag
source setup/patchesToSource.sh
cd $CMSSW_BASE/src
scram b -j 10
cd ..
