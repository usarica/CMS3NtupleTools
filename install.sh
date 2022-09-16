#!/bin/bash

# USER INPUTS
CMS3Tag=combined
CMSSW_release=CMSSW_10_2_22
CMSSW_release_name=    #Leave this blank if you don't know what it is.  It's just a marker in case you have multiple identical directories. No need for the underscore.
SCRAM_ARCH_name="amd64_gcc700" # Leave slc6/7 out


MACHINESPECS="$(uname -a)"
echo "Machine specifics: ${MACHINESPECS}"
declare -i FOUND_EL6=0
if [[ "${MACHINESPECS}" == *"el6"* ]]; then
  FOUND_EL6=1
else
  for evar in $(env); do
    if [[ "$evar" == *"SINGULARITY_IMAGE_HUMAN"* ]]; then
      # This means you are running in a condor job with a singularity image loaded.
      if [[ "$evar" == *"rhel6"* ]] || [[ "$evar" == *"slc6"* ]]; then
        FOUND_EL6=1
      fi
    fi
  done
fi


if [[ ${FOUND_EL6} -eq 1 ]]; then
  SCRAM_ARCH_name="slc6_${SCRAM_ARCH_name}"
else
  SCRAM_ARCH_name="slc7_${SCRAM_ARCH_name}"
fi

export SCRAM_ARCH=${SCRAM_ARCH_name}

if [[ -z ${CMSSW_release_name+x} ]]; then
  CMSSW_release_name="${CMSSW_release}"
else
  CMSSW_release_name="${CMSSW_release}_${CMSSW_release_name}"
fi

#--Here there be dragons----
export CMS_PATH=/cvmfs/cms.cern.ch
source /cvmfs/cms.cern.ch/cmsset_default.sh
scramv1 p -n ${CMSSW_release_name} CMSSW $CMSSW_release
cd ${CMSSW_release_name}/src
eval $(scramv1 runtime -sh)

# new upstream-only ignores user's cmssw, but makes cms-init much, much faster
git cms-init --upstream-only

# Preliminary electron scale and smearing corrections according to https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes#2018_Preliminary_Energy_Correcti
# We need the ElectronTools package to calculate smear and scale uncertainties so just download the ScaleAndSmearing files manualy 
git cms-merge-topic cms-egamma:EgammaPostRecoTools
git cms-merge-topic cms-egamma:PhotonIDValueMapSpeedup1029
git cms-merge-topic cms-egamma:slava77-btvDictFix_10210
git cms-addpkg EgammaAnalysis/ElectronTools
(rm -rf EgammaAnalysis/ElectronTools/data; git clone git@github.com:cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data;)

# Add HZZ re-trained electron MVA 'id+iso's
git cms-merge-topic mkovac:Electron_XGBoost_MVA_2016_and_2018_CMSSW_10_2_15

# For MET recipe for 2017 EE noise fix
# Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_0_for_M
git cms-merge-topic cms-met:METFixEE2017_949_v2_backport_to_102X

## For reading 2 electron iso branches in 102X from 94X sample
## https://github.com/cms-sw/cmssw/issues/25573
## NOTE this should be taken out once merged/backported into CMSSW
## Merged after CMSSW_10_2_11
#git cms-merge-topic Sam-Harper:IORulesForPFClusIso_1025

# e/gamma cluster production
git cms-addpkg RecoEcal/EgammaClusterProducers

# PU jet id updates
# Twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID#9X_2017_and_10X_2018_recipes
git cms-merge-topic alefisico:PUID_102X # For CMSSW_10_2_X older than CMSSW_10_2_21. Need PUID_106X for CMSSW_10_6_X older than CMSSW_10_6_11
git cms-addpkg RecoJets/JetProducers
(
  git clone -b 94X_weights_DYJets_inc_v2 git@github.com:cms-jet/PUjetID.git RecoJets/JetProducers/data/PUJetIDweights
  cp -r RecoJets/JetProducers/data/PUJetIDweights/weights/pileupJetId_{94,102}X_Eta* RecoJets/JetProducers/data/
  rm -rf RecoJets/JetProducers/data/PUJetIDweights
)


#######################################
# No CMSSW packages beyond this point #
#######################################

# MELA
git clone git@github.com:JHUGen/JHUGenMELA.git
./JHUGenMELA/setup.sh -j
eval $(./JHUGenMELA/setup.sh env)

# MELA Analytics
git clone git@github.com:MELALabs/MelaAnalytics.git
./MelaAnalytics/setup.sh -j
eval $(./MelaAnalytics/setup.sh env)

# Common LHE tools
git clone git@github.com:usarica/CommonLHETools.git

# IvyFramework
git clone git@github.com:IvyFramework/IvyDataTools.git IvyFramework/IvyDataTools
./IvyFramework/IvyDataTools/setup.sh -j
git clone git@github.com:IvyFramework/IvyAutoMELA.git IvyFramework/IvyAutoMELA
./IvyFramework/IvyAutoMELA/setup.sh -j
eval $(./IvyFramework/IvyAutoMELA/setup.sh env)

#########################
#  DeepAK8 fat jet tagger
# #######################
# check out the package - note, need ssh key in gitlab.cern.ch
# because this is top secret code that needs to be password protected apparently
# and thus, the user must either configure ssh keys or manually type their password.
# the latter ruins the whole "run this install script, get a coffee, use the ntuplemaker" workflow.
if [[ -d /nfs-7/userdata/NtupleModules/NNKit_ver_2018-03-08_for94X ]]; then
  cp -r /nfs-7/userdata/NtupleModules/NNKit_ver_2018-03-08_for94X NNKit
elif [[ $(readlink -f ~)/CMS3/NNKit_ver_2018-03-08_for94X ]]; then
  cp -r $(readlink -f ~)/CMS3/NNKit_ver_2018-03-08_for94X NNKit
else
  git clone ssh://git@gitlab.cern.ch:7999/DeepAK8/NNKit.git -b ver_2018-03-08_for94X
fi
# setup mxnet library
# cp /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/cmssw/CMSSW_10_2_0/config/toolbox/$SCRAM_ARCH/tools/selected/mxnet-predict.xml $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected
scram setup mxnet-predict
# rm $CMSSW_BASE/external/$SCRAM_ARCH/lib/libmxnet_predict.so
# cp NNKit/misc/lib/libmxnet_predict.so $CMSSW_BASE/external/$SCRAM_ARCH/lib/libmxnet_predict.so
# #######################

git clone git@github.com:usarica/CMS3NtupleTools.git CMS3
(
  cd CMS3; git checkout ${CMS3Tag}; cd -
  cd $CMSSW_BASE/src/CMS3/NtupleMaker/data/JECs ; . download.sh; cd -
  cd $CMSSW_BASE/src/CMS3/NtupleMaker/data/JERs ; . download.sh; cd -
)
./CMS3/NtupleMaker/setup.sh -j

# copy json files to test directory (or wherever you are doing cmsRun)
# cp NNKit/data/ak8/*.{json,params} $CMSSW_BASE/src/CMS3/NtupleMaker/test/

scram b -j

# see comment in patchesToSource.sh
rm $CMSSW_BASE/lib/$SCRAM_ARCH/.poisonededmplugincache
