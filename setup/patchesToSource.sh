#!/bin/bash


#############
# MVA JetId #
#############
 
git clone https://github.com/latinos/UserCode-CMG-CMGTools-External $CMSSW_BASE/src/CMGTools/External
pushd $CMSSW_BASE/src/CMGTools/External
git checkout V00-03-01
popd


#######################
# MC Jet Flavor Truth #
#######################

git clone https://github.com/cms-cvs-history/PhysicsTools-JetMCAlgos $CMSSW_BASE/src/PhysicsTools/JetMCAlgos
pushd $CMSSW_BASE/src/PhysicsTools/JetMCAlgos
git checkout PhysicsTools-JetMCAlgos-V00-07-04
popd

sed '13 s/ParameterSet/Utilities/' <$CMSSW_BASE/src/PhysicsTools/JetMCAlgos/plugins/CandOneToManyDeltaRMatcher.cc >blah
mv blah $CMSSW_BASE/src/PhysicsTools/JetMCAlgos/plugins/CandOneToManyDeltaRMatcher.cc

sed '12 s/ParameterSet/Utilities/' <$CMSSW_BASE/src/PhysicsTools/JetMCAlgos/plugins/CandOneToOneDeltaRMatcher.cc >blah
mv blah $CMSSW_BASE/src/PhysicsTools/JetMCAlgos/plugins/CandOneToOneDeltaRMatcher.cc

sed '13 s/ParameterSet/Utilities/' <$CMSSW_BASE/src/PhysicsTools/JetMCAlgos/plugins/GenJetBCEnergyRatio.cc >blah
mv blah $CMSSW_BASE/src/PhysicsTools/JetMCAlgos/plugins/GenJetBCEnergyRatio.cc

sed '15 s/ParameterSet/Utilities/' <$CMSSW_BASE/src/PhysicsTools/JetMCAlgos/plugins/JetFlavourIdentifier.cc >blah
mv blah $CMSSW_BASE/src/PhysicsTools/JetMCAlgos/plugins/JetFlavourIdentifier.cc

sed '15 s/ParameterSet/Utilities/' <$CMSSW_BASE/src/PhysicsTools/JetMCAlgos/plugins/JetPartonMatcher.cc >blah
mv blah $CMSSW_BASE/src/PhysicsTools/JetMCAlgos/plugins/JetPartonMatcher.cc

sed '12 s/ParameterSet/Utilities/' <$CMSSW_BASE/src/PhysicsTools/JetMCAlgos/plugins/PartonSelector.cc >blah
mv blah $CMSSW_BASE/src/PhysicsTools/JetMCAlgos/plugins/PartonSelector.cc

sed '32 s/GenParticleCandidate/GenParticle/' <$CMSSW_BASE/src/PhysicsTools/JetMCAlgos/plugins/GenJetBCEnergyRatio.cc >blah
mv blah $CMSSW_BASE/src/PhysicsTools/JetMCAlgos/plugins/GenJetBCEnergyRatio.cc


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
