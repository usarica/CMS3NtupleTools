#!/bin/bash

#USER INPUTS
CMS3Tag=diet
CMSSW_release=CMSSW_9_2_0_patch2
CMSSW_release_name=    #Leave this blank if you don't know what it is.  It's just a marker in case you have multiple identical directories.  Don't forget the underscore!
export SCRAM_ARCH=slc6_amd64_gcc530

#--Here there be dragons----
export CMS_PATH=/cvmfs/cms.cern.ch
source /cvmfs/cms.cern.ch/cmsset_default.sh
scramv1 p -n ${CMSSW_release}${CMSSW_release_name} CMSSW $CMSSW_release
cd ${CMSSW_release}${CMSSW_release_name}/src
eval `scramv1 runtime -sh`
git cms-init
git clone git@github.com:cmstas/NtupleMaker.git CMS3/NtupleMaker
cd CMS3/NtupleMaker
git checkout $CMS3Tag
source setup/patchesToSource.sh
cd $CMSSW_BASE/src
scram b -j 20
cd ..

# see comment in patchesToSource.sh
rm $CMSSW_BASE/lib/$SCRAM_ARCH/.poisonededmplugincache
