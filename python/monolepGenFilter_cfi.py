import FWCore.ParameterSet.Config as cms

monolepGenFilter = cms.EDFilter("LepGenFilter",
        ntupleDaughters       = cms.bool(True),
	nGenLepsRequired = cms.int32(1)
)
