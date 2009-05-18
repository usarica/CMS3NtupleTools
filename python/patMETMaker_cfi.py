import FWCore.ParameterSet.Config as cms

patMETMaker = cms.EDFilter("PATMETMaker",
                           patMETsInputTag = cms.InputTag("layer1METs")
)


