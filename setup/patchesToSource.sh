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
