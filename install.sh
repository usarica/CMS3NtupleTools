#!/bin/bash

#USER INPUTS
CMS3Tag=combined
CMSSW_release=CMSSW_10_2_15
CMSSW_release_name=    #Leave this blank if you don't know what it is.  It's just a marker in case you have multiple identical directories.  Don't forget the underscore!
export SCRAM_ARCH=slc6_amd64_gcc700

#--Here there be dragons----
export CMS_PATH=/cvmfs/cms.cern.ch
source /cvmfs/cms.cern.ch/cmsset_default.sh
scramv1 p -n ${CMSSW_release}${CMSSW_release_name} CMSSW $CMSSW_release
cd ${CMSSW_release}${CMSSW_release_name}
eval $(scramv1 runtime -sh)

# new upstream-only ignores user's cmssw, but makes cms-init much, much faster
git cms-init --upstream-only

#Preliminary electron scale and smearing corrections according to https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes#2018_Preliminary_Energy_Correcti
#We need the ElectronTools package to calculate smear and scale uncertainties so just download the ScaleAndSmearing files manualy 
git cms-merge-topic cms-egamma:EgammaPostRecoTools
git cms-merge-topic cms-egamma:PhotonIDValueMapSpeedup1029
git cms-merge-topic cms-egamma:slava77-btvDictFix_10210
git cms-addpkg EgammaAnalysis/ElectronTools
(rm -rf EgammaAnalysis/ElectronTools/data; git clone git@github.com:cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data;)

# For MET recipe for 2017 EE noise fix
# Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_0_for_M
git cms-merge-topic cms-met:METFixEE2017_949_v2_backport_to_102X

## For reading 2 electron iso branches in 102X from 94X sample
## https://github.com/cms-sw/cmssw/issues/25573
## NOTE this should be taken out once merged/backported into CMSSW
## Merged after CMSSW_10_2_11
#git cms-merge-topic Sam-Harper:IORulesForPFClusIso_1025

cd src

git clone git@github.com:usarica/CMS3-NtupleMaker.git CMS3/NtupleMaker
cd CMS3/NtupleMaker
git checkout $CMS3Tag
source setup/patchesToSource.sh

#######################################
# No CMSSW packages beyond this point #
#######################################

# MELA
git clone git@github.com:cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
(cd ZZMatrixElement; source setup.sh -j;)
#(cd ZZMatrixElement; git fetch; git checkout -b from-v222 v2.2.2; source setup.sh -j 12;)

# MELA Analytics
git clone git@github.com:usarica/MelaAnalytics.git
#(cd MelaAnalytics; git fetch; git checkout -b from-v19 v1.9)

# Common LHE tools
git clone git@github.com:usarica/CommonLHETools.git
#(cd CommonLHETools; git fetch; git checkout -b from-v132 v1.3.2)

cd $CMSSW_BASE/src
scram b -j
cd ..

# see comment in patchesToSource.sh
rm $CMSSW_BASE/lib/$SCRAM_ARCH/.poisonededmplugincache

#cd $CMSSW_BASE/src/CMS3/NtupleMaker/test/
