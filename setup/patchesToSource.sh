#!/bin/bash


#############
# MVA JetId #
#############
 
git clone https://github.com/latinos/UserCode-CMG-CMGTools-External $CMSSW_BASE/src/CMGTools/External
pushd $CMSSW_BASE/src/CMGTools/External
git checkout V00-03-01
popd


#######################
# LCG dictionaries #
#######################

git clone https://github.com/cmstas/Dictionaries $CMSSW_BASE/src/CMS2/Dictionaries


#######################

# run checkdeps
printf "\nchecking deps:\n"
git cms-checkdeps -a

# compile
scram build -c -j 20


#######################
#   miniAOD specific  #
#######################

cmsenv
pushd .
cd $CMSSW_BASE/src
git init
git config core.sparsecheckout true
echo DataFormats/PatCandidates >> .git/info/sparse-checkout
echo PhysicsTools/PatAlgos >> .git/info/sparse-checkout
echo RecoJets/JetProducers >> .git/info/sparse-checkout
git remote add -f origin https://github.com/gpetruc/cmssw.git
git pull origin micro-from700
popd
