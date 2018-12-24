#!/bin/bash

#USER INPUTS
CMS3Tag=combined
CMSSW_release=CMSSW_10_2_5
CMSSW_release_name=    #Leave this blank if you don't know what it is.  It's just a marker in case you have multiple identical directories.  Don't forget the underscore!
export SCRAM_ARCH=slc6_amd64_gcc700

#--Here there be dragons----
export CMS_PATH=/cvmfs/cms.cern.ch
source /cvmfs/cms.cern.ch/cmsset_default.sh
scramv1 p -n ${CMSSW_release}${CMSSW_release_name} CMSSW $CMSSW_release
cd ${CMSSW_release}${CMSSW_release_name}
eval `scramv1 runtime -sh`

# new upstream-only ignores user's cmssw, but makes cms-init much, much faster
git cms-init --upstream-only

# For MET recipe for 2017 EE noise fix
git cms-merge-topic cms-met:METFixEE2017_949_v2_backport_to_102X

cd src

git clone git@github.com:cmstas/NtupleMaker.git CMS3/NtupleMaker
cd CMS3/NtupleMaker
git checkout $CMS3Tag
source setup/patchesToSource.sh

#######################################
# No CMSSW packages beyond this point #
#######################################

# MELA
git clone git@github.com:cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
(cd ZZMatrixElement; git checkout -b from-v217b1 v2.1.7b1; source setup.sh -j 12;)

# MELA Analytics
git clone git@github.com:usarica/MelaAnalytics.git
(cd MelaAnalytics; git checkout -b from-v11 v1.1)

# Common LHE tools
git clone git@github.com:usarica/CommonLHETools.git
(cd CommonLHETools; git checkout -b from-v122 v1.2.2)

cd $CMSSW_BASE/src
scram b -j 25
cd ..

# see comment in patchesToSource.sh
rm $CMSSW_BASE/lib/$SCRAM_ARCH/.poisonededmplugincache

cd $CMSSW_BASE/src/CMS3/NtupleMaker/test/
