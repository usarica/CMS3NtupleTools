#!/bin/bash

####################
# MET Filters 2012 #
####################

git cms-cvs-history import V00-00-10 RecoMET/METFilters
touch RecoMET/METFilters/data/dummy.txt
sed '3 s/EDFilter/EDProducer/' <$CMSSW_BASE/src/RecoMET/METFilters/python/EcalDeadCellDeltaRFilter_cfi.py > blah
mv blah $CMSSW_BASE/src/RecoMET/METFilters/python/EcalDeadCellDeltaRFilter_cfi.py

# not working -- tag not migrated (using cvs):
# git cms-cvs-history import V00-00-08 RecoMET/METAnalyzers
cvs co -r V00-00-08 RecoMET/METAnalyzers

git cms-cvs-history import V00-03-23 CommonTools/RecoAlgos  

#############
# MVA JetId #
#############
 
# not working -- tag not migrated: git cms-cvs-history import V00-03-01 -d CMGTools/External UserCode/CMG/CMGTools/External
cvs co -r V00-03-01 -d CMGTools/External UserCode/CMG/CMGTools/External

#############################
# do stuff with PF candidates
#############################

git cms-cvs-history import V15-02-06 RecoParticleFlow/PFProducer

# use new recommended tags (keep the original PFProducer)
git cms-cvs-history import V00-03-15 CommonTools/ParticleFlow
git cms-cvs-history import V08-09-21 PhysicsTools/PatAlgos

####################
# 2012 Electron ID #
####################

git clone https://github.com/latinos/UserCode-EGamma-EGammaAnalysisTools.git $CMSSW_BASE/src/EGamma/EGammaAnalysisTools
pushd $CMSSW_BASE/src/EGamma/EGammaAnalysisTools
git checkout CutBasedId_V00-00-03 
popd

########################
# Type-I Corrected MET #
########################

# not working -- tag not migrated (using cvs): git cms-cvs-history
#git cms-cvs-history import V04-06-09 JetMETCorrections/Type1MET 
cvs co -r V04-06-09 JetMETCorrections/Type1MET 

#######################
# MC Jet Flavor Truth #
#######################

# not working -- tag not migrated (using cvs): git cms-cvs-history
# git cms-cvs-history import V00-07-04 PhysicsTools/JetMCAlgos
cvs co -r V00-07-04 PhysicsTools/JetMCAlgos

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
# TAUs
#######################

git cms-cvs-history import V01-04-23 RecoTauTag/RecoTau #HCP + new discriminants

# not working -- tag not migrated (using cvs): git cms-cvs-history
# git cms-cvs-history import V01-04-10 RecoTauTag/Configuration
cvs co -r V01-04-10 RecoTauTag/Configuration

# not working -- tag not migrated (using cvs): git cms-cvs-history
# git cms-cvs-history import V00-04-00 CondFormats/EgammaObjects
cvs co -r V00-04-00 CondFormats/EgammaObjects

#######################

# run checkdeps
printf "\nchecking deps:\n"
git cms-checkdeps -a

# compile
scram build -c -j 20
