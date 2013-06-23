import FWCore.ParameterSet.Config as cms

dilepGenFilter = cms.EDFilter("LepGenFilter",
        ntupleDaughters       = cms.bool(True),
	nGenLepsRequired = cms.int32(2)
)
