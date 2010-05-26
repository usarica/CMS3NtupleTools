import FWCore.ParameterSet.Config as cms

beamHaloMaker = cms.EDFilter("BeamHaloMaker",
                             beamHaloInputTag = cms.InputTag("BeamHaloSummary"),
                             cscHaloInputTag  = cms.InputTag("CSCHaloData")
)


