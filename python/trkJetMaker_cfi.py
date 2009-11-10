import FWCore.ParameterSet.Config as cms

trkJetMaker = cms.EDProducer("TrkJetMaker", 
                             trkJetsInputTag = cms.InputTag('ak5TrackJets'),
                             trkJetPtCut     = cms.double(3.0)
)


