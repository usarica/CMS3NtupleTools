import FWCore.ParameterSet.Config as cms

monolepGenFilter = cms.EDProducer("LepGenFilter",
	nGenLepsRequired = cms.int32(1)
)
