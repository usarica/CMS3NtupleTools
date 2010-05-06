import FWCore.ParameterSet.Config as cms

prunedUncorrectedCMS2TrackJets = cms.EDFilter("TrackJetCollectionPruner",
   inputUncorrectedJetCollection = cms.InputTag("ak5TrackJets"),
   uncorrectedJetPtCut           = cms.double(3.0) ##cut on the UNCORRECTED reco jets!!!!!
)
