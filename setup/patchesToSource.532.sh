#!/bin/bash

####################
# MET Filters 2012 #
####################

git clone https://github.com/cms-cvs-history/RecoMET-METFilters $CMSSW_BASE/src/RecoMET/METFilters
pushd $CMSSW_BASE/src/RecoMET/METFilters
git checkout RecoMET-METFilters-V00-00-10
popd
touch $CMSSW_BASE/src/RecoMET/METFilters/data/dummy.txt
sed '3 s/EDFilter/EDProducer/' <$CMSSW_BASE/src/RecoMET/METFilters/python/EcalDeadCellDeltaRFilter_cfi.py > blah
mv blah $CMSSW_BASE/src/RecoMET/METFilters/python/EcalDeadCellDeltaRFilter_cfi.py


git clone https://github.com/cms-cvs-dump/RecoMET_METAnalyzers $CMSSW_BASE/src/RecoMET/METAnalyzers
pushd $CMSSW_BASE/src/RecoMET/METAnalyzers
git checkout V00-00-08
popd



git clone https://github.com/cms-cvs-history/CommonTools-RecoAlgos $CMSSW_BASE/src/CommonTools/RecoAlgos
pushd $CMSSW_BASE/src/CommonTools/RecoAlgos
git checkout CommonTools-RecoAlgos-V00-03-23
popd


#############
# MVA JetId #
#############
 
git clone https://github.com/latinos/UserCode-CMG-CMGTools-External $CMSSW_BASE/src/CMGTools/External
pushd $CMSSW_BASE/src/CMGTools/External
git checkout V00-03-01
popd


#############################
# do stuff with PF candidates
#############################

git clone https://github.com/cms-cvs-history/RecoParticleFlow-PFProducer $CMSSW_BASE/src/RecoParticleFlow/PFProducer
pushd $CMSSW_BASE/src/RecoParticleFlow/PFProducer
git checkout RecoParticleFlow-PFProducer-V15-02-06
popd

git clone https://github.com/cms-cvs-history/CommonTools-ParticleFlow $CMSSW_BASE/src/CommonTools/ParticleFlow
pushd $CMSSW_BASE/src/CommonTools/ParticleFlow
git checkout CommonTools-ParticleFlow-V00-03-15
popd

git clone https://github.com/cms-cvs-history/PhysicsTools-PatAlgos $CMSSW_BASE/src/PhysicsTools/PatAlgos
pushd $CMSSW_BASE/src/PhysicsTools/PatAlgos
git checkout PhysicsTools-PatAlgos-V08-09-21
popd


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

git clone https://github.com/cms-cvs-history/JetMETCorrections-Type1MET $CMSSW_BASE/src/JetMETCorrections/Type1MET
pushd $CMSSW_BASE/src/JetMETCorrections/Type1MET
git checkout JetMETCorrections-Type1MET-V04-06-09
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
# TAUs
#######################

git clone https://github.com/cms-cvs-history/RecoTauTag-RecoTau $CMSSW_BASE/src/RecoTauTag/RecoTau
pushd $CMSSW_BASE/src/RecoTauTag/RecoTau
git checkout RecoTauTag-RecoTau-V01-04-23
popd


git clone https://github.com/cms-cvs-history/RecoTauTag-Configuration $CMSSW_BASE/src/RecoTauTag/Configuration
pushd $CMSSW_BASE/src/RecoTauTag/Configuration
git checkout RecoTauTag-Configuration-V01-04-10
popd

git clone https://github.com/cms-cvs-history/CondFormats-EgammaObjects $CMSSW_BASE/src/CondFormats/EgammaObjects
pushd $CMSSW_BASE/src/CondFormats/EgammaObjects
git checkout CondFormats-EgammaObjects-V00-04-00
popd


#######################
#  LCG dictionaries   #
#######################

git clone https://github.com/cmstas/Dictionaries $CMSSW_BASE/src/CMS2/Dictionaries


#######################

# run checkdeps
printf "\nchecking deps:\n"
git cms-checkdeps -a

# compile
scram build -c -j 20
