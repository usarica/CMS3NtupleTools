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

# pT-dependent JERs and phi-dependent JECs
git cms-merge-topic 09e9a17cb26020bc507dffec70137bacc6f124c6 # This is the commit id directly from cmssw

## For reading 2 electron iso branches in 102X from 94X sample
## https://github.com/cms-sw/cmssw/issues/25573
## NOTE this should be taken out once merged/backported into CMSSW
## Merged after CMSSW_10_2_11
#git cms-merge-topic Sam-Harper:IORulesForPFClusIso_1025

cd src

git clone git@github.com:usarica/CMS3-NtupleMaker.git CMS3/NtupleMaker
(cd CMS3/NtupleMaker; git checkout ${CMS3Tag}; )

git clone git@github.com:cmstas/Dictionaries CMS3/Dictionaries

mkdir $CMSSW_BASE/bullshit
mv $CMSSW_BASE/src/* $CMSSW_BASE/bullshit/
git cms-addpkg RecoEcal/EgammaClusterProducers
mv $CMSSW_BASE/bullshit/* $CMSSW_BASE/src/
rmdir $CMSSW_BASE/bullshit

#########################
#  DeepAK8 fat jet tagger
# #######################
# check out the package - note, need ssh key in gitlab.cern.ch
# because this is top secret code that needs to be password protected apparently
# and thus, the user must either configure ssh keys or manually type their password.
# the latter ruins the whole "run this install script, get a coffee, use the ntuplemaker" workflow.
# git clone ssh://git@gitlab.cern.ch:7999/DeepAK8/NNKit.git -b ver_2018-03-08_for94X
cp -r /nfs-7/userdata/NtupleModules/NNKit_ver_2018-03-08_for94X NNKit
# setup mxnet library
# cp /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/cmssw/CMSSW_10_2_0/config/toolbox/$SCRAM_ARCH/tools/selected/mxnet-predict.xml $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected
scram setup mxnet-predict
# rm $CMSSW_BASE/external/$SCRAM_ARCH/lib/libmxnet_predict.so
# cp NNKit/misc/lib/libmxnet_predict.so $CMSSW_BASE/external/$SCRAM_ARCH/lib/libmxnet_predict.so
# copy json files to test directory (or wherever you are doing cmsRun)
# cp NNKit/data/ak8/*.{json,params} $CMSSW_BASE/src/CMS3/NtupleMaker/test/
# #######################


#######################################
# No CMSSW packages beyond this point #
#######################################

# MELA
git clone git@github.com:cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
(cd ZZMatrixElement; source setup.sh -j;)
#(cd ZZMatrixElement; git fetch; git checkout -b from-v222 v2.2.2; source setup.sh -j;)

# MELA Analytics
git clone git@github.com:usarica/MelaAnalytics.git
#(cd MelaAnalytics; git fetch; git checkout -b from-v19 v1.9)

# Common LHE tools
git clone git@github.com:usarica/CommonLHETools.git
#(cd CommonLHETools; git fetch; git checkout -b from-v132 v1.3.2)

# CMSDataTools
git clone git@github.com:usarica/CMSDataTools.git
#(cd CMSDataTools; git fetch; git checkout -b from-v111 v1.1.1)

(cd $CMSSW_BASE/src/CMS3/NtupleMaker/data/JECs ; . download.sh)


cd $CMSSW_BASE/src
scram b -j
cd ..

# see comment in patchesToSource.sh
rm $CMSSW_BASE/lib/$SCRAM_ARCH/.poisonededmplugincache

#cd $CMSSW_BASE/src/CMS3/NtupleMaker/test/
