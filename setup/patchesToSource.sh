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


#############
# MVA JetId #
#############
 
git clone https://github.com/latinos/UserCode-CMG-CMGTools-External $CMSSW_BASE/src/CMGTools/External
pushd $CMSSW_BASE/src/CMGTools/External
git checkout V00-03-01
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
# LCG dictionaries #
#######################

git clone https://github.com/cmstas/Dictionaries $CMSSW_BASE/src/CMS2/Dictionaries


#######################

# run checkdeps
printf "\nchecking deps:\n"
git cms-checkdeps -a

# compile
scram build -c -j 20
