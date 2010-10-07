import FWCore.ParameterSet.Config as cms

monolepGenFilter = cms.EDFilter("LepGenFilter",
	nGenLepsRequired = cms.int32(1)
)
