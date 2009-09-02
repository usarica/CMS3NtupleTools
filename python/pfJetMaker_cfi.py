import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDFilter("PFJetMaker",
             pfJetsInputTag   = cms.InputTag("sisCone5PFJets")
)


