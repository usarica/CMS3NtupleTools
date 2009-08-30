import FWCore.ParameterSet.Config as cms

trkJetMaker = cms.EDProducer("TrkJetMaker", 
                             trkJetsInputTag = cms.InputTag('SISCone5TrkJets')
)

trkjetmaker = cms.Sequence(trkJetMaker)
