#!/bin/bash

#CMSSW Environment
: ${CMSSW_BASE:?"[ERROR] patchesToSource.sh: \$CMSSW_BASE not set. Perhaps you forgot to set your CMSSW environment."} 

#cd
cd $CMSSW_BASE/src/CMS3/NtupleMaker

##############
## MVA JetId #
##############
 
git clone https://github.com/latinos/UserCode-CMG-CMGTools-External $CMSSW_BASE/src/CMGTools/External
pushd $CMSSW_BASE/src/CMGTools/External
git checkout V00-03-01
rm plugins/PileupJetIdProducer.cc
popd


# #######################
# # LCG dictionaries #
# #######################

git clone https://github.com/cmstas/Dictionaries $CMSSW_BASE/src/CMS3/Dictionaries

# ####################
# # jet tool box     #
# ####################

git clone https://github.com/cms-jet/JetToolbox $CMSSW_BASE/src/JMEAnalysis/JetToolbox -b jetToolbox_74X
pushd $CMSSW_BASE/src/JMEAnalysis/JetToolbox
git checkout a80553163684d718e65e44d9b385a4aa11475659
popd

mkdir $CMSSW_BASE/bullshit  
mv $CMSSW_BASE/src/* $CMSSW_BASE/bullshit/
git cms-addpkg   RecoEcal/EgammaClusterProducers
mv $CMSSW_BASE/bullshit/* $CMSSW_BASE/src/
rmdir $CMSSW_BASE/bullshit

####### line needs to be added ###############
inputfile="$CMSSW_BASE/src/RecoEcal/EgammaClusterProducers/src/PFECALSuperClusterProducer.cc"
grep "desc.setAllowAnything();" $inputfile 2>&1 > /dev/null
doNothing=$?
if [ ! $doNothing = "0" ]; 
then
  echo "line does not exist. Adding now."
  sed -i 's/edm::ParameterSetDescription desc;/edm::ParameterSetDescription desc;\n  desc.setAllowAnything();/' $inputfile
else
  echo "line already exists. File $inputfile will be unchanged."
fi

####################

#######################

# run checkdeps
printf "\nchecking deps:\n"
git cms-checkdeps -a

########################
#  EGM MVA ID 80X (must be run after checkdeps to avoid checking out useless packages)
#######################
cd $CMSSW_BASE
mkdir $CMSSW_BASE/bullshit  
mv $CMSSW_BASE/src/* $CMSSW_BASE/bullshit/
git cms-merge-topic ikrav:egm_id_80X_v2
mv $CMSSW_BASE/src/RecoEgamma/ $CMSSW_BASE/bullshit/
rm -rf $CMSSW_BASE/src/*
mv $CMSSW_BASE/bullshit/* $CMSSW_BASE/src/
rmdir $CMSSW_BASE/bullshit

cd $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/data/
git clone https://github.com/ikrav/RecoEgamma-ElectronIdentification.git new
cd new
git checkout egm_id_80X_v1
mv Spring16* ../
cd ../
rm -rf new
cd $CMSSW_BASE/src

### End of EGM MVA ID 80X ###

######################
## # DeepFlavour # ###
######################
#  Btagging DeepFlavour recipe for 80X (from https://twiki.cern.ch/twiki/bin/viewauth/CMS/DeepFlavour)
pushd $CMSSW_BASE/src
git cms-merge-topic -u mverzett:DeepFlavour-from-CMSSW_8_0_21
scram b -j 20
mkdir -p RecoBTag/DeepFlavour/data/
cd RecoBTag/DeepFlavour/data/
wget http://home.fnal.gov/~verzetti//DeepFlavour/training/DeepFlavourNoSL.json
popd
### END of DeepFlavour ###


# compile
cd $CMSSW_BASE/src
scram b -j 20

# If you see an error with the words "poisoned" and "plugin", then this is what you need to do
# I don't know what the hell it does, but you just need to do it
# ...after...every...single...scram b...before you can run
rm $CMSSW_BASE/lib/$SCRAM_ARCH/.poisonededmplugincache

