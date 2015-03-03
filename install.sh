#!/bin/bash

#USER INPUTS
CMS3Tag=master

export SCRAM_ARCH=slc6_amd64_gcc481
scramv1 p -n CMSSW_7_2_0 CMSSW CMSSW_7_2_0
cd CMSSW_7_2_0/src
cmsenv
git clone git@github.com:cmstas/NtupleMaker.git CMS3/NtupleMaker
cd CMS3/NtupleMaker
git checkout $CMS3Tag
source setup/patchesToSource.sh
cd $CMSSW_BASE/src
scram b -j 10
cd ..
