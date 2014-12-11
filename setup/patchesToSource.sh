#!/bin/bash


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


#######################

# run checkdeps
printf "\nchecking deps:\n"
git cms-checkdeps -a

# compile
scram build -c -j 20


#######################
#   miniAOD specific  #
#######################

#cmsenv
#pushd .
#cd $CMSSW_BASE/src
#git init
#git config core.sparsecheckout true
#echo CommonTools/CandAlgos                    >> .git/info/sparse-checkout
#echo Configuration/Applications		      >> .git/info/sparse-checkout
#echo Configuration/EventContent		      >> .git/info/sparse-checkout
#echo Configuration/PyReleaseValidation	      >> .git/info/sparse-checkout
#echo Configuration/StandardSequences	      >> .git/info/sparse-checkout
#echo DataFormats/Candidate		      >> .git/info/sparse-checkout
#echo DataFormats/PatCandidates		      >> .git/info/sparse-checkout
#echo DataFormats/ParticleFlowCandidate	      >> .git/info/sparse-checkout
#echo DataFormats/TrackReco		      >> .git/info/sparse-checkout
#echo DataFormats/EgammaReco		      >> .git/info/sparse-checkout
#echo DataFormats/ParticleFlowCandidate	      >> .git/info/sparse-checkout
#echo DataFormats/EgammaCandidates	      >> .git/info/sparse-checkout
#echo DataFormats/BTauReco		      >> .git/info/sparse-checkout
#echo DataFormats/JetReco		      >> .git/info/sparse-checkout
#echo DQM/SiStripMonitorSummary		      >> .git/info/sparse-checkout
#echo JetMETCorrections/Type1MET		      >> .git/info/sparse-checkout
#echo PhysicsTools/HepMCCandAlgos	      >> .git/info/sparse-checkout
#echo PhysicsTools/JetMCAlgos		      >> .git/info/sparse-checkout
#echo PhysicsTools/JetMCUtils		      >> .git/info/sparse-checkout
#echo PhysicsTools/PatAlgos		      >> .git/info/sparse-checkout
#echo PhysicsTools/PatUtils		      >> .git/info/sparse-checkout
#echo RecoJets/JetProducers		      >> .git/info/sparse-checkout
#echo RecoMET/METAlgorithms		      >> .git/info/sparse-checkout
#echo RecoEcal/EgammaCoreTools		      >> .git/info/sparse-checkout
#echo RecoEgamma/EgammaPhotonProducers	      >> .git/info/sparse-checkout
#echo RecoEgamma/EgammaTools		      >> .git/info/sparse-checkout
#echo RecoBTag/SecondaryVertex		      >> .git/info/sparse-checkout
#echo RecoEcal/EgammaClusterProducers	      >> .git/info/sparse-checkout
#echo RecoVertex/AdaptiveVertexFinder          >> .git/info/sparse-checkout     
#echo SimDataFormats/JetMatching               >> .git/info/sparse-checkout
#git remote add -f origin https://github.com/gpetruc/cmssw.git
#git pull origin miniAOD-from704-part3
#git checkout c64fdc9b0f5fe842291e36e284db8823d108f614
#popd
