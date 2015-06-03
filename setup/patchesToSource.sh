#!/bin/bash

#CMSSW Environment
: ${CMSSW_BASE:?"[ERROR] patchesToSource.sh: \$CMSSW_BASE not set. Perhaps you forgot to set your CMSSW environment."} 

#cd
cd $CMSSW_BASE/src/CMS3/NtupleMaker

#Hard-code lepton ID
git update-index --assume-unchanged setup/GsfEleFull5x5SigmaIEtaIEtaCut72X.cc 
git update-index --assume-unchanged setup/cutBasedElectronID_PHYS14_PU20bx25_V0_miniAOD_cff.py 
mkdir $CMSSW_BASE/bullshit  
mv $CMSSW_BASE/src/* $CMSSW_BASE/bullshit/
git cms-addpkg RecoEgamma/ElectronIdentification 
mv $CMSSW_BASE/bullshit/CMS3/NtupleMaker/setup/GsfEleFull5x5SigmaIEtaIEtaCut72X.cc $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/plugins/cuts/
mv $CMSSW_BASE/bullshit/CMS3/NtupleMaker/setup/cutBasedElectronID_PHYS14_PU20bx25_V0_miniAOD_cff.py $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/python/Identification
pushd $CMSSW_BASE/src/
scram b -j 20 
popd
mv $CMSSW_BASE/bullshit/* $CMSSW_BASE/src/
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

git clone https://github.com/cms-jet/JetToolbox JMEAnalysis/JetToolbox -b jetToolbox_74X

####### line needs to be added ###############
git cms-addpkg   RecoEcal/EgammaClusterProducers
inputfile="RecoEcal/EgammaClusterProducers/src/PFECALSuperClusterProducer.cc"
grep "desc.setAllowAnything();" $inputfile 2>&1 > NULL
doNothing=$?
if [ ! $doNothing = "0" ]; then
echo "line does not exist. Adding now."
sed -e '/^ edm::ParameterSetDescription desc;/a\ \ desc.setAllowAnything();' $inputfile > temp_inputfile.txt
cat temp_inputfile.txt > $inputfile
# deletes temp file
if [ -e temp_inputfile.txt ]; then
rm temp_inputfile.txt
fi
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
