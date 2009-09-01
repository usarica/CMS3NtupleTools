import FWCore.ParameterSet.Config as cms

hypGenMaker = cms.EDFilter(
	"HypGenMaker",
	candToGenAssTag = cms.InputTag("candToGenAssMaker"),
	hypInputTag = cms.InputTag("hypDilepMaker")
)

