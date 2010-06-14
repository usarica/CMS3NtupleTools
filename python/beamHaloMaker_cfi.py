import FWCore.ParameterSet.Config as cms

beamHaloMaker = cms.EDProducer("BeamHaloMaker",
                             beamHaloInputTag = cms.InputTag("BeamHaloSummary"),
                             cscHaloInputTag  = cms.InputTag("CSCHaloData")
)


