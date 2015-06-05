#!/bin/bash

#CMSSW Environment
: ${CMSSW_BASE:?"[ERROR] patchesToSource.sh: \$CMSSW_BASE not set. Perhaps you forgot to set your CMSSW environment."} 

#cd
cd $CMSSW_BASE/src/CMS3/NtupleMaker

#Hard-code lepton ID
#git update-index --assume-unchanged setup/GsfEleFull5x5SigmaIEtaIEtaCut72X.cc 
mkdir $CMSSW_BASE/bullshit  
mv $CMSSW_BASE/src/* $CMSSW_BASE/bullshit/
pushd $CMSSW_BASE/src/
#git update-index --assume-unchanged setup/cutBasedElectronID_PHYS14_PU20bx25_V2_miniAOD_cff.py 
git cms-merge-topic ikrav:egm_id_74X_v0
mv $CMSSW_BASE/bullshit/* $CMSSW_BASE/src/
popd
#git cms-addpkg RecoEgamma/ElectronIdentification 
#mv $CMSSW_BASE/bullshit/CMS3/NtupleMaker/setup/GsfEleFull5x5SigmaIEtaIEtaCut72X.cc $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/plugins/cuts/
rmdir $CMSSW_BASE/bullshit


#############
# MVA JetId #
#############
 
git clone https://github.com/latinos/UserCode-CMG-CMGTools-External $CMSSW_BASE/src/CMGTools/External
pushd $CMSSW_BASE/src/CMGTools/External
git checkout V00-03-01
rm plugins/PileupJetIdProducer.cc
popd


#######################
# LCG dictionaries #
#######################

git clone https://github.com/cmstas/Dictionaries $CMSSW_BASE/src/CMS3/Dictionaries

####################
# jet tool box     #
####################

git clone https://github.com/cms-jet/JetToolbox $CMSSW_BASE/src/JMEAnalysis/JetToolbox -b jetToolbox_74X

####### line needs to be added ###############
git cms-addpkg   RecoEcal/EgammaClusterProducers
inputfile="RecoEcal/EgammaClusterProducers/src/PFECALSuperClusterProducer.cc"
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

# compile
scram build -c -j 20
