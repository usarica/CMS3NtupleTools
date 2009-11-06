import FWCore.ParameterSet.Config as cms

dilepGenFilter = cms.EDFilter("LepGenFilter",
	nGenLepsRequired = cms.int32(2)
)
