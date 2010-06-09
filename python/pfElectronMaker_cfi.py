import FWCore.ParameterSet.Config as cms

pfElectronMaker = cms.EDFilter("PFElectronMaker",
                               pfCandidatesTag = cms.InputTag("particleFlow","electrons")
)
