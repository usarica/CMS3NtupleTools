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
echo CommonTools/CandAlgos >> .git/info/sparse-checkout
echo DataFormats/Candidate >> .git/info/sparse-checkout
echo DataFormats/PatCandidates >> .git/info/sparse-checkout
echo DataFormats/ParticleFlowCandidate >> .git/info/sparse-checkout
echo DataFormats/TrackReco >> .git/info/sparse-checkout
echo JetMETCorrections/Type1MET >> .git/info/sparse-checkout
echo PhysicsTools/HepMCCandAlgos >> .git/info/sparse-checkout
echo PhysicsTools/JetMCAlgos >> .git/info/sparse-checkout
echo PhysicsTools/PatAlgos >> .git/info/sparse-checkout
echo PhysicsTools/PatUtils >> .git/info/sparse-checkout
echo RecoJets/JetProducers >> .git/info/sparse-checkout
echo RecoMET/METAlgorithms >> .git/info/sparse-checkout
echo DataFormats/ParticleFlowCandidate >> .git/info/sparse-checkout
git remote add -f origin https://github.com/gpetruc/cmssw.git
git pull origin micro-from700
git checkout 6ce5fed6d4f11ec52f867bd96187549b6b0f73b5
popd
